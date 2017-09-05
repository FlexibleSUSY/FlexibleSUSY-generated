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

// File generated at Tue 5 Sep 2017 11:02:33

/**
 * @file MRSSMtower_mass_eigenstates.cpp
 * @brief implementation of the MRSSMtower model class
 *
 * Contains the definition of the MRSSMtower model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 5 Sep 2017 11:02:33 with FlexibleSUSY
 * 1.7.5 (git commit: c98e024e1e74ea3309b68f7006d5f91f8df6c678) and SARAH 4.12.0 .
 */

#include "MRSSMtower_mass_eigenstates.hpp"
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

#define CLASSNAME MRSSMtower_mass_eigenstates

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

CLASSNAME::MRSSMtower_mass_eigenstates(const MRSSMtower_input_parameters& input_)
   : MRSSMtower_soft_parameters(input_)
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
   , problems(MRSSMtower_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSRdp(0), MSRum(0)
      , MsigmaO(0), MphiO(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(
      Eigen::Array<double,3,1>::Zero()), MSu(Eigen::Array<double,6,1>::Zero()),
      MSe(Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,4,1>::Zero())
      , MAh(Eigen::Array<double,4,1>::Zero()), MRh(Eigen::Array<double,2,1>::Zero(
      )), MHpm(Eigen::Array<double,4,1>::Zero()), MChi(Eigen::Array<double,4,1>
      ::Zero()), MCha1(Eigen::Array<double,2,1>::Zero()), MCha2(Eigen::Array<
      double,2,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MFd(
      Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,3,1>::Zero()),
      MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,4,4>::Zero()), ZA(Eigen::Matrix<double,4,
      4>::Zero()), ZHR(Eigen::Matrix<double,2,2>::Zero()), ZP(Eigen::Matrix<double
      ,4,4>::Zero()), ZN1(Eigen::Matrix<std::complex<double>,4,4>::Zero()), ZN2(
      Eigen::Matrix<std::complex<double>,4,4>::Zero()), UM1(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), UP1(Eigen::Matrix<std::complex<double>,2
      ,2>::Zero()), UM2(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP2(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())


{
}

CLASSNAME::~MRSSMtower_mass_eigenstates()
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

const MRSSMtower_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MRSSMtower_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MRSSMtower_physical& physical_)
{
   physical = physical_;
}

const Problems<MRSSMtower_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<MRSSMtower_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
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
   MRSSMtower_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mHd2(gsl_vector_get(x, 0));
   model->set_mHu2(gsl_vector_get(x, 1));
   model->set_vS(gsl_vector_get(x, 2));
   model->set_vT(gsl_vector_get(x, 3));


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

   mHd2 = solver->get_solution(0);
   mHu2 = solver->get_solution(1);
   vS = solver->get_solution(2);
   vT = solver->get_solution(3);


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

   const double old_vT = vT;
   const double old_vS = vS;
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   vT = Re((0.2*(-20*g2*MDWBT*Power(vd,4)*AbsSqr(LamSD) + 20*g2*MDWBT*Power(vu,
      4)*AbsSqr(LamSU) - 5.477225575051661*g1*LamTD*MDBS*Power(vd,4)*Conj(LamSD) -
      5.477225575051661*g1*LamTU*MDBS*Power(vu,4)*Conj(LamSU) - 5.477225575051661
      *g1*LamSD*MDBS*Power(vd,4)*Conj(LamTD) - 10*MuD*Power(vd,4)*AbsSqr(LamSD)*
      Conj(LamTD) - 5.477225575051661*g1*LamSU*MDBS*Power(vu,4)*Conj(LamTU) + 10*
      MuU*Power(vu,4)*AbsSqr(LamSU)*Conj(LamTU) - 5.477225575051661*g1*LamTD*Power
      (vd,4)*Conj(LamSD)*Conj(MDBS) - 5.477225575051661*g1*LamTU*Power(vu,4)*Conj(
      LamSU)*Conj(MDBS) - 5.477225575051661*g1*LamSD*Power(vd,4)*Conj(LamTD)*Conj(
      MDBS) - 5.477225575051661*g1*LamSU*Power(vu,4)*Conj(LamTU)*Conj(MDBS) - 20*
      g2*Power(vd,4)*AbsSqr(LamSD)*Conj(MDWBT) + 20*g2*Power(vu,4)*AbsSqr(LamSU)*
      Conj(MDWBT) - 10*LamTD*Power(vd,4)*AbsSqr(LamSD)*Conj(MuD) + 10*LamTU*Power(
      vu,4)*AbsSqr(LamSU)*Conj(MuU) + 10*Power(vd,4)*Conj(LamTD)*Conj(MuD)*Sqr(
      LamSD) - 10*Power(vu,4)*Conj(LamTU)*Conj(MuU)*Sqr(LamSU) - 40*g2*MDWBT*mS2*
      Sqr(vd) - 40*mS2*MuD*Conj(LamTD)*Sqr(vd) - 40*g2*mS2*Conj(MDWBT)*Sqr(vd) -
      40*LamTD*mS2*Conj(MuD)*Sqr(vd) - 160*g2*MDWBT*Sqr(MDBS)*Sqr(vd) - 160*MuD*
      Conj(LamTD)*Sqr(MDBS)*Sqr(vd) - 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vd) - 160*
      LamTD*Conj(MuD)*Sqr(MDBS)*Sqr(vd) + 40*g2*MDWBT*mS2*Sqr(vu) + 40*mS2*MuU*
      Conj(LamTU)*Sqr(vu) + 40*g2*mS2*Conj(MDWBT)*Sqr(vu) + 40*LamTU*mS2*Conj(MuU)
      *Sqr(vu) + 160*g2*MDWBT*Sqr(MDBS)*Sqr(vu) + 160*MuU*Conj(LamTU)*Sqr(MDBS)*
      Sqr(vu) + 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vu) + 160*LamTU*Conj(MuU)*Sqr(
      MDBS)*Sqr(vu) + 20*g2*MDWBT*AbsSqr(LamSD)*Sqr(vd)*Sqr(vu) - 20*g2*MDWBT*
      AbsSqr(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTD*MDBS*Conj(LamSD)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*MDBS*Conj(LamSU)*Sqr(vd)*Sqr(vu
      ) - 10*LamTU*MuD*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuU*Conj
      (LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSD*MDBS*Conj(
      LamTD)*Sqr(vd)*Sqr(vu) - 20*MuD*AbsSqr(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) +
      10*LamSD*MuU*Conj(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*
      LamSU*MDBS*Conj(LamTU)*Sqr(vd)*Sqr(vu) + 20*MuU*AbsSqr(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) - 10*LamSU*MuD*Conj(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTD*Conj(LamSD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) + 20*g2*
      AbsSqr(LamSD)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*g2*AbsSqr(LamSU)*Conj(MDWBT)*
      Sqr(vd)*Sqr(vu) - 20*LamTD*AbsSqr(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*
      LamSD*LamTU*Conj(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*LamSD*LamSU*Conj(
      LamTU)*Conj(MuD)*Sqr(vd)*Sqr(vu) + 20*LamTU*AbsSqr(LamSD)*Conj(MuU)*Sqr(vd)*
      Sqr(vu) + 10*LamSU*LamTD*Conj(LamSD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamSD*
      LamSU*Conj(LamTD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuD*Power(vd,4)*Sqr(
      Conj(LamSD)) - 10*LamTU*MuU*Power(vu,4)*Sqr(Conj(LamSU))))/(32*mS2*mT2 + 2*
      Power(vd,4)*AbsSqr(LamSD)*AbsSqr(LamTD) + 2*Power(vu,4)*AbsSqr(LamSU)*AbsSqr
      (LamTU) + 128*mT2*Sqr(MDBS) + 128*mS2*Sqr(MDWBT) + 512*Sqr(MDBS)*Sqr(MDWBT)
      + 16*mT2*AbsSqr(LamSD)*Sqr(vd) + 8*mS2*AbsSqr(LamTD)*Sqr(vd) + 32*AbsSqr(
      LamTD)*Sqr(MDBS)*Sqr(vd) + 64*AbsSqr(LamSD)*Sqr(MDWBT)*Sqr(vd) + 16*mT2*
      AbsSqr(LamSU)*Sqr(vu) + 8*mS2*AbsSqr(LamTU)*Sqr(vu) + 32*AbsSqr(LamTU)*Sqr(
      MDBS)*Sqr(vu) + 64*AbsSqr(LamSU)*Sqr(MDWBT)*Sqr(vu) + 4*AbsSqr(LamSU)*AbsSqr
      (LamTD)*Sqr(vd)*Sqr(vu) + 4*AbsSqr(LamSD)*AbsSqr(LamTU)*Sqr(vd)*Sqr(vu) + 2*
      LamTD*LamTU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 2*LamSD*LamTU*Conj(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 2*LamSU*LamTD*Conj(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) + 2*LamSD*LamSU*Conj(LamTD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) -
      Power(vd,4)*Sqr(LamTD)*Sqr(Conj(LamSD)) - Power(vu,4)*Sqr(LamTU)*Sqr(Conj(
      LamSU)) - Power(vd,4)*Sqr(LamSD)*Sqr(Conj(LamTD)) - Power(vu,4)*Sqr(LamSU)*
      Sqr(Conj(LamTU))));
   vS = Re((-1.4142135623730951*(4*mT2*vT + 16*vT*Sqr(MDWBT) + g2*MDWBT*Sqr(vd)
      + vT*AbsSqr(LamTD)*Sqr(vd) + MuD*Conj(LamTD)*Sqr(vd) + g2*Conj(MDWBT)*Sqr(
      vd) + LamTD*Conj(MuD)*Sqr(vd) - g2*MDWBT*Sqr(vu) + vT*AbsSqr(LamTU)*Sqr(vu)
      - MuU*Conj(LamTU)*Sqr(vu) - g2*Conj(MDWBT)*Sqr(vu) - LamTU*Conj(MuU)*Sqr(vu)
      ))/(LamTD*Conj(LamSD)*Sqr(vd) + LamSD*Conj(LamTD)*Sqr(vd) - LamTU*Conj(LamSU
      )*Sqr(vu) - LamSU*Conj(LamTU)*Sqr(vu)));
   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*
      vd*AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS
      *Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(
      LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*
      Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40
      *vu*AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*
      vu*Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*
      vu*Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr
      (vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu*AbsSqr(
      LamTU)*Sqr(vT)))/vu);

   const bool is_finite = IsFinite(vT) && IsFinite(vS) && IsFinite(mHd2) &&
      IsFinite(mHu2);

   if (!is_finite) {
      vT = old_vT;
      vS = old_vS;
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = 0;

   const double old_vT = vT;
   const double old_vS = vS;
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   vT = Re((0.2*(-20*g2*MDWBT*Power(vd,4)*AbsSqr(LamSD) + 20*g2*MDWBT*Power(vu,
      4)*AbsSqr(LamSU) - 5.477225575051661*g1*LamTD*MDBS*Power(vd,4)*Conj(LamSD) -
      5.477225575051661*g1*LamTU*MDBS*Power(vu,4)*Conj(LamSU) - 5.477225575051661
      *g1*LamSD*MDBS*Power(vd,4)*Conj(LamTD) - 10*MuD*Power(vd,4)*AbsSqr(LamSD)*
      Conj(LamTD) - 5.477225575051661*g1*LamSU*MDBS*Power(vu,4)*Conj(LamTU) + 10*
      MuU*Power(vu,4)*AbsSqr(LamSU)*Conj(LamTU) - 5.477225575051661*g1*LamTD*Power
      (vd,4)*Conj(LamSD)*Conj(MDBS) - 5.477225575051661*g1*LamTU*Power(vu,4)*Conj(
      LamSU)*Conj(MDBS) - 5.477225575051661*g1*LamSD*Power(vd,4)*Conj(LamTD)*Conj(
      MDBS) - 5.477225575051661*g1*LamSU*Power(vu,4)*Conj(LamTU)*Conj(MDBS) - 20*
      g2*Power(vd,4)*AbsSqr(LamSD)*Conj(MDWBT) + 20*g2*Power(vu,4)*AbsSqr(LamSU)*
      Conj(MDWBT) - 10*LamTD*Power(vd,4)*AbsSqr(LamSD)*Conj(MuD) + 10*LamTU*Power(
      vu,4)*AbsSqr(LamSU)*Conj(MuU) + 10*Power(vd,4)*Conj(LamTD)*Conj(MuD)*Sqr(
      LamSD) - 10*Power(vu,4)*Conj(LamTU)*Conj(MuU)*Sqr(LamSU) - 40*g2*MDWBT*mS2*
      Sqr(vd) - 40*mS2*MuD*Conj(LamTD)*Sqr(vd) - 40*g2*mS2*Conj(MDWBT)*Sqr(vd) -
      40*LamTD*mS2*Conj(MuD)*Sqr(vd) - 160*g2*MDWBT*Sqr(MDBS)*Sqr(vd) - 160*MuD*
      Conj(LamTD)*Sqr(MDBS)*Sqr(vd) - 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vd) - 160*
      LamTD*Conj(MuD)*Sqr(MDBS)*Sqr(vd) + 40*g2*MDWBT*mS2*Sqr(vu) + 40*mS2*MuU*
      Conj(LamTU)*Sqr(vu) + 40*g2*mS2*Conj(MDWBT)*Sqr(vu) + 40*LamTU*mS2*Conj(MuU)
      *Sqr(vu) + 160*g2*MDWBT*Sqr(MDBS)*Sqr(vu) + 160*MuU*Conj(LamTU)*Sqr(MDBS)*
      Sqr(vu) + 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vu) + 160*LamTU*Conj(MuU)*Sqr(
      MDBS)*Sqr(vu) + 20*g2*MDWBT*AbsSqr(LamSD)*Sqr(vd)*Sqr(vu) - 20*g2*MDWBT*
      AbsSqr(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTD*MDBS*Conj(LamSD)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*MDBS*Conj(LamSU)*Sqr(vd)*Sqr(vu
      ) - 10*LamTU*MuD*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuU*Conj
      (LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSD*MDBS*Conj(
      LamTD)*Sqr(vd)*Sqr(vu) - 20*MuD*AbsSqr(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) +
      10*LamSD*MuU*Conj(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*
      LamSU*MDBS*Conj(LamTU)*Sqr(vd)*Sqr(vu) + 20*MuU*AbsSqr(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) - 10*LamSU*MuD*Conj(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTD*Conj(LamSD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) + 20*g2*
      AbsSqr(LamSD)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*g2*AbsSqr(LamSU)*Conj(MDWBT)*
      Sqr(vd)*Sqr(vu) - 20*LamTD*AbsSqr(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*
      LamSD*LamTU*Conj(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*LamSD*LamSU*Conj(
      LamTU)*Conj(MuD)*Sqr(vd)*Sqr(vu) + 20*LamTU*AbsSqr(LamSD)*Conj(MuU)*Sqr(vd)*
      Sqr(vu) + 10*LamSU*LamTD*Conj(LamSD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamSD*
      LamSU*Conj(LamTD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuD*Power(vd,4)*Sqr(
      Conj(LamSD)) - 10*LamTU*MuU*Power(vu,4)*Sqr(Conj(LamSU))))/(32*mS2*mT2 + 2*
      Power(vd,4)*AbsSqr(LamSD)*AbsSqr(LamTD) + 2*Power(vu,4)*AbsSqr(LamSU)*AbsSqr
      (LamTU) + 128*mT2*Sqr(MDBS) + 128*mS2*Sqr(MDWBT) + 512*Sqr(MDBS)*Sqr(MDWBT)
      + 16*mT2*AbsSqr(LamSD)*Sqr(vd) + 8*mS2*AbsSqr(LamTD)*Sqr(vd) + 32*AbsSqr(
      LamTD)*Sqr(MDBS)*Sqr(vd) + 64*AbsSqr(LamSD)*Sqr(MDWBT)*Sqr(vd) + 16*mT2*
      AbsSqr(LamSU)*Sqr(vu) + 8*mS2*AbsSqr(LamTU)*Sqr(vu) + 32*AbsSqr(LamTU)*Sqr(
      MDBS)*Sqr(vu) + 64*AbsSqr(LamSU)*Sqr(MDWBT)*Sqr(vu) + 4*AbsSqr(LamSU)*AbsSqr
      (LamTD)*Sqr(vd)*Sqr(vu) + 4*AbsSqr(LamSD)*AbsSqr(LamTU)*Sqr(vd)*Sqr(vu) + 2*
      LamTD*LamTU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 2*LamSD*LamTU*Conj(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 2*LamSU*LamTD*Conj(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) + 2*LamSD*LamSU*Conj(LamTD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) -
      Power(vd,4)*Sqr(LamTD)*Sqr(Conj(LamSD)) - Power(vu,4)*Sqr(LamTU)*Sqr(Conj(
      LamSU)) - Power(vd,4)*Sqr(LamSD)*Sqr(Conj(LamTD)) - Power(vu,4)*Sqr(LamSU)*
      Sqr(Conj(LamTU))));
   vS = Re((-1.4142135623730951*(4*mT2*vT + 16*vT*Sqr(MDWBT) + g2*MDWBT*Sqr(vd)
      + vT*AbsSqr(LamTD)*Sqr(vd) + MuD*Conj(LamTD)*Sqr(vd) + g2*Conj(MDWBT)*Sqr(
      vd) + LamTD*Conj(MuD)*Sqr(vd) - g2*MDWBT*Sqr(vu) + vT*AbsSqr(LamTU)*Sqr(vu)
      - MuU*Conj(LamTU)*Sqr(vu) - g2*Conj(MDWBT)*Sqr(vu) - LamTU*Conj(MuU)*Sqr(vu)
      ))/(LamTD*Conj(LamSD)*Sqr(vd) + LamSD*Conj(LamTD)*Sqr(vd) - LamTU*Conj(LamSU
      )*Sqr(vu) - LamSU*Conj(LamTU)*Sqr(vu)));
   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*
      vd*AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS
      *Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(
      LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*
      Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40
      *vu*AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*
      vu*Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*
      vu*Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr
      (vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu*AbsSqr(
      LamTU)*Sqr(vT)))/vu);

   const bool is_finite = IsFinite(vT) && IsFinite(vS) && IsFinite(mHd2) &&
      IsFinite(mHu2);

   if (!is_finite) {
      vT = old_vT;
      vS = old_vS;
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      error = 1;
   }


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

   x_init[0] = mHd2;
   x_init[1] = mHu2;
   x_init[2] = vS;
   x_init[3] = vT;


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
      tadpole[2] += Re(tadpole_hh(3));
      tadpole[3] += Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {

      }
   }

   double vT;
   double vS;
   double mHd2;
   double mHu2;

   vT = Re((0.2*(-20*g2*MDWBT*Power(vd,4)*AbsSqr(LamSD) + 20*g2*MDWBT*Power(vu,
      4)*AbsSqr(LamSU) - 5.477225575051661*g1*LamTD*MDBS*Power(vd,4)*Conj(LamSD) -
      5.477225575051661*g1*LamTU*MDBS*Power(vu,4)*Conj(LamSU) - 5.477225575051661
      *g1*LamSD*MDBS*Power(vd,4)*Conj(LamTD) - 10*MuD*Power(vd,4)*AbsSqr(LamSD)*
      Conj(LamTD) - 5.477225575051661*g1*LamSU*MDBS*Power(vu,4)*Conj(LamTU) + 10*
      MuU*Power(vu,4)*AbsSqr(LamSU)*Conj(LamTU) - 5.477225575051661*g1*LamTD*Power
      (vd,4)*Conj(LamSD)*Conj(MDBS) - 5.477225575051661*g1*LamTU*Power(vu,4)*Conj(
      LamSU)*Conj(MDBS) - 5.477225575051661*g1*LamSD*Power(vd,4)*Conj(LamTD)*Conj(
      MDBS) - 5.477225575051661*g1*LamSU*Power(vu,4)*Conj(LamTU)*Conj(MDBS) - 20*
      g2*Power(vd,4)*AbsSqr(LamSD)*Conj(MDWBT) + 20*g2*Power(vu,4)*AbsSqr(LamSU)*
      Conj(MDWBT) - 10*LamTD*Power(vd,4)*AbsSqr(LamSD)*Conj(MuD) + 10*LamTU*Power(
      vu,4)*AbsSqr(LamSU)*Conj(MuU) + 160*mS2*tadpole[2] + 10*Power(vd,4)*Conj(
      LamTD)*Conj(MuD)*Sqr(LamSD) - 10*Power(vu,4)*Conj(LamTU)*Conj(MuU)*Sqr(LamSU
      ) + 640*tadpole[2]*Sqr(MDBS) - 40*g2*MDWBT*mS2*Sqr(vd) - 40*mS2*MuD*Conj(
      LamTD)*Sqr(vd) - 40*g2*mS2*Conj(MDWBT)*Sqr(vd) - 40*LamTD*mS2*Conj(MuD)*Sqr(
      vd) + 80*AbsSqr(LamSD)*tadpole[2]*Sqr(vd) - 28.284271247461902*LamTD*Conj(
      LamSD)*tadpole[3]*Sqr(vd) - 28.284271247461902*LamSD*Conj(LamTD)*tadpole[3]*
      Sqr(vd) - 160*g2*MDWBT*Sqr(MDBS)*Sqr(vd) - 160*MuD*Conj(LamTD)*Sqr(MDBS)*Sqr
      (vd) - 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vd) - 160*LamTD*Conj(MuD)*Sqr(MDBS)*
      Sqr(vd) + 40*g2*MDWBT*mS2*Sqr(vu) + 40*mS2*MuU*Conj(LamTU)*Sqr(vu) + 40*g2*
      mS2*Conj(MDWBT)*Sqr(vu) + 40*LamTU*mS2*Conj(MuU)*Sqr(vu) + 80*AbsSqr(LamSU)*
      tadpole[2]*Sqr(vu) + 28.284271247461902*LamTU*Conj(LamSU)*tadpole[3]*Sqr(vu)
      + 28.284271247461902*LamSU*Conj(LamTU)*tadpole[3]*Sqr(vu) + 160*g2*MDWBT*
      Sqr(MDBS)*Sqr(vu) + 160*MuU*Conj(LamTU)*Sqr(MDBS)*Sqr(vu) + 160*g2*Conj(
      MDWBT)*Sqr(MDBS)*Sqr(vu) + 160*LamTU*Conj(MuU)*Sqr(MDBS)*Sqr(vu) + 20*g2*
      MDWBT*AbsSqr(LamSD)*Sqr(vd)*Sqr(vu) - 20*g2*MDWBT*AbsSqr(LamSU)*Sqr(vd)*Sqr(
      vu) + 5.477225575051661*g1*LamTD*MDBS*Conj(LamSD)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTU*MDBS*Conj(LamSU)*Sqr(vd)*Sqr(vu) - 10*LamTU*MuD*
      Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuU*Conj(LamSD)*Conj(
      LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSD*MDBS*Conj(LamTD)*Sqr(vd)
      *Sqr(vu) - 20*MuD*AbsSqr(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 10*LamSD*MuU*
      Conj(LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSU*MDBS*
      Conj(LamTU)*Sqr(vd)*Sqr(vu) + 20*MuU*AbsSqr(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(
      vu) - 10*LamSU*MuD*Conj(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTD*Conj(LamSD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Sqr(vd)*Sqr(vu) + 20*g2*
      AbsSqr(LamSD)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*g2*AbsSqr(LamSU)*Conj(MDWBT)*
      Sqr(vd)*Sqr(vu) - 20*LamTD*AbsSqr(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*
      LamSD*LamTU*Conj(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*LamSD*LamSU*Conj(
      LamTU)*Conj(MuD)*Sqr(vd)*Sqr(vu) + 20*LamTU*AbsSqr(LamSD)*Conj(MuU)*Sqr(vd)*
      Sqr(vu) + 10*LamSU*LamTD*Conj(LamSD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamSD*
      LamSU*Conj(LamTD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuD*Power(vd,4)*Sqr(
      Conj(LamSD)) - 10*LamTU*MuU*Power(vu,4)*Sqr(Conj(LamSU))))/(32*mS2*mT2 + 2*
      Power(vd,4)*AbsSqr(LamSD)*AbsSqr(LamTD) + 2*Power(vu,4)*AbsSqr(LamSU)*AbsSqr
      (LamTU) + 128*mT2*Sqr(MDBS) + 128*mS2*Sqr(MDWBT) + 512*Sqr(MDBS)*Sqr(MDWBT)
      + 16*mT2*AbsSqr(LamSD)*Sqr(vd) + 8*mS2*AbsSqr(LamTD)*Sqr(vd) + 32*AbsSqr(
      LamTD)*Sqr(MDBS)*Sqr(vd) + 64*AbsSqr(LamSD)*Sqr(MDWBT)*Sqr(vd) + 16*mT2*
      AbsSqr(LamSU)*Sqr(vu) + 8*mS2*AbsSqr(LamTU)*Sqr(vu) + 32*AbsSqr(LamTU)*Sqr(
      MDBS)*Sqr(vu) + 64*AbsSqr(LamSU)*Sqr(MDWBT)*Sqr(vu) + 4*AbsSqr(LamSU)*AbsSqr
      (LamTD)*Sqr(vd)*Sqr(vu) + 4*AbsSqr(LamSD)*AbsSqr(LamTU)*Sqr(vd)*Sqr(vu) + 2*
      LamTD*LamTU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 2*LamSD*LamTU*Conj(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 2*LamSU*LamTD*Conj(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) + 2*LamSD*LamSU*Conj(LamTD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) -
      Power(vd,4)*Sqr(LamTD)*Sqr(Conj(LamSD)) - Power(vu,4)*Sqr(LamTU)*Sqr(Conj(
      LamSU)) - Power(vd,4)*Sqr(LamSD)*Sqr(Conj(LamTD)) - Power(vu,4)*Sqr(LamSU)*
      Sqr(Conj(LamTU))));
   vS = Re((-1.4142135623730951*(4*mT2*vT - 4*tadpole[2] + 16*vT*Sqr(MDWBT) +
      g2*MDWBT*Sqr(vd) + vT*AbsSqr(LamTD)*Sqr(vd) + MuD*Conj(LamTD)*Sqr(vd) + g2*
      Conj(MDWBT)*Sqr(vd) + LamTD*Conj(MuD)*Sqr(vd) - g2*MDWBT*Sqr(vu) + vT*AbsSqr
      (LamTU)*Sqr(vu) - MuU*Conj(LamTU)*Sqr(vu) - g2*Conj(MDWBT)*Sqr(vu) - LamTU*
      Conj(MuU)*Sqr(vu)))/(LamTD*Conj(LamSD)*Sqr(vd) + LamSD*Conj(LamTD)*Sqr(vd) -
      LamTU*Conj(LamSU)*Sqr(vu) - LamSU*Conj(LamTU)*Sqr(vu)));
   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*
      vd*AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS
      *Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) + 40*tadpole[0] - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) -
      20*vd*AbsSqr(LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr
      (vu) + 5*vd*Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40
      *vu*AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*
      vu*Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*
      vu*Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) + 40*tadpole[1] - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) +
      3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) -
      10*vu*AbsSqr(LamTU)*Sqr(vT)))/vu);

   const bool is_finite = IsFinite(vT) && IsFinite(vS) && IsFinite(mHd2) &&
      IsFinite(mHu2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mHd2;
   ewsb_parameters[1] = mHu2;
   ewsb_parameters[2] = vS;
   ewsb_parameters[3] = vT;


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
   MRSSMtower_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mHd2 = gsl_vector_get(x, 0);
   const double mHu2 = gsl_vector_get(x, 1);
   const double vS = gsl_vector_get(x, 2);
   const double vT = gsl_vector_get(x, 3);

   model->set_mHd2(mHd2);
   model->set_mHu2(mHu2);
   model->set_vS(vS);
   model->set_vT(vT);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = mHd2;
   ewsb_parameters[1] = mHu2;
   ewsb_parameters[2] = vS;
   ewsb_parameters[3] = vT;


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
           "MRSSMtower\n"
           "========================================\n";
   MRSSMtower_soft_parameters::print(ostr);
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
   const auto old_vS = vS;
   const auto old_vT = vT;

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

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;
   vS = old_vS;
   vT = old_vT;

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
   std::future<void> fut_MGlu;
   std::future<void> fut_MFv;
   std::future<void> fut_MSRdp;
   std::future<void> fut_MSRum;
   std::future<void> fut_MsigmaO;
   std::future<void> fut_MphiO;
   std::future<void> fut_MVP;
   std::future<void> fut_MVZ;
   std::future<void> fut_MSd;
   std::future<void> fut_MSv;
   std::future<void> fut_MSu;
   std::future<void> fut_MSe;
   std::future<void> fut_Mhh;
   std::future<void> fut_MAh;
   std::future<void> fut_MRh;
   std::future<void> fut_MHpm;
   std::future<void> fut_MChi;
   std::future<void> fut_MCha1;
   std::future<void> fut_MCha2;
   std::future<void> fut_MFe;
   std::future<void> fut_MFd;
   std::future<void> fut_MFu;
   std::future<void> fut_MVWm;

   if (calculate_bsm_pole_masses) {
      fut_MAh = run_async([obj_ptr] () { obj_ptr->calculate_MAh_pole(); });
      fut_MCha1 = run_async([obj_ptr] () { obj_ptr->calculate_MCha1_pole(); });
      fut_MCha2 = run_async([obj_ptr] () { obj_ptr->calculate_MCha2_pole(); });
      fut_MChi = run_async([obj_ptr] () { obj_ptr->calculate_MChi_pole(); });
      fut_MGlu = run_async([obj_ptr] () { obj_ptr->calculate_MGlu_pole(); });
      fut_Mhh = run_async([obj_ptr] () { obj_ptr->calculate_Mhh_pole(); });
      fut_MHpm = run_async([obj_ptr] () { obj_ptr->calculate_MHpm_pole(); });
      fut_MphiO = run_async([obj_ptr] () { obj_ptr->calculate_MphiO_pole(); });
      fut_MRh = run_async([obj_ptr] () { obj_ptr->calculate_MRh_pole(); });
      fut_MSd = run_async([obj_ptr] () { obj_ptr->calculate_MSd_pole(); });
      fut_MSe = run_async([obj_ptr] () { obj_ptr->calculate_MSe_pole(); });
      fut_MsigmaO = run_async([obj_ptr] () { obj_ptr->calculate_MsigmaO_pole(); });
      fut_MSRdp = run_async([obj_ptr] () { obj_ptr->calculate_MSRdp_pole(); });
      fut_MSRum = run_async([obj_ptr] () { obj_ptr->calculate_MSRum_pole(); });
      fut_MSu = run_async([obj_ptr] () { obj_ptr->calculate_MSu_pole(); });
      fut_MSv = run_async([obj_ptr] () { obj_ptr->calculate_MSv_pole(); });
   }

   if (calculate_sm_pole_masses) {
      fut_MVG = run_async([obj_ptr] () { obj_ptr->calculate_MVG_pole(); });
      fut_MFv = run_async([obj_ptr] () { obj_ptr->calculate_MFv_pole(); });
      fut_MVP = run_async([obj_ptr] () { obj_ptr->calculate_MVP_pole(); });
      fut_MVZ = run_async([obj_ptr] () { obj_ptr->calculate_MVZ_pole(); });
      fut_MFe = run_async([obj_ptr] () { obj_ptr->calculate_MFe_pole(); });
      fut_MFd = run_async([obj_ptr] () { obj_ptr->calculate_MFd_pole(); });
      fut_MFu = run_async([obj_ptr] () { obj_ptr->calculate_MFu_pole(); });
      fut_MVWm = run_async([obj_ptr] () { obj_ptr->calculate_MVWm_pole(); });
   }

   if (fut_MAh.valid()) fut_MAh.get();
   if (fut_MCha1.valid()) fut_MCha1.get();
   if (fut_MCha2.valid()) fut_MCha2.get();
   if (fut_MChi.valid()) fut_MChi.get();
   if (fut_MGlu.valid()) fut_MGlu.get();
   if (fut_Mhh.valid()) fut_Mhh.get();
   if (fut_MHpm.valid()) fut_MHpm.get();
   if (fut_MphiO.valid()) fut_MphiO.get();
   if (fut_MRh.valid()) fut_MRh.get();
   if (fut_MSd.valid()) fut_MSd.get();
   if (fut_MSe.valid()) fut_MSe.get();
   if (fut_MsigmaO.valid()) fut_MsigmaO.get();
   if (fut_MSRdp.valid()) fut_MSRdp.get();
   if (fut_MSRum.valid()) fut_MSRum.get();
   if (fut_MSu.valid()) fut_MSu.get();
   if (fut_MSv.valid()) fut_MSv.get();
   if (fut_MVG.valid()) fut_MVG.get();
   if (fut_MFv.valid()) fut_MFv.get();
   if (fut_MVP.valid()) fut_MVP.get();
   if (fut_MVZ.valid()) fut_MVZ.get();
   if (fut_MFe.valid()) fut_MFe.get();
   if (fut_MFd.valid()) fut_MFd.get();
   if (fut_MFu.valid()) fut_MFu.get();
   if (fut_MVWm.valid()) fut_MVWm.get();

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
   if (PHYSICAL(MSRdp) < 0.) problems.flag_tachyon(MRSSMtower_info::SRdp);
   if (PHYSICAL(MSRum) < 0.) problems.flag_tachyon(MRSSMtower_info::SRum);
   if (PHYSICAL(MsigmaO) < 0.) problems.flag_tachyon(MRSSMtower_info::sigmaO);
   if (PHYSICAL(MphiO) < 0.) problems.flag_tachyon(MRSSMtower_info::phiO);
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Sd);
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Sv);
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Su);
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Se);
   if (PHYSICAL(Mhh).tail<4>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::hh);
   if (PHYSICAL(MAh).tail<3>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Ah);
   if (PHYSICAL(MRh).tail<2>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Rh);
   if (PHYSICAL(MHpm).tail<3>().minCoeff() < 0.) problems.flag_tachyon(MRSSMtower_info::Hpm);

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
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   MRSSMtower_soft_parameters::clear();
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

std::string CLASSNAME::name() const
{
   return "MRSSMtower";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   MRSSMtower_soft_parameters::run_to(scale, eps);
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
      problems.flag_tachyon(MRSSMtower_info::SRdp);
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
      problems.flag_tachyon(MRSSMtower_info::SRum);
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
      problems.flag_tachyon(MRSSMtower_info::sigmaO);
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
      problems.flag_tachyon(MRSSMtower_info::phiO);
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
   problems.flag_bad_mass(MRSSMtower_info::Sd, eigenvalue_error >
      precision * Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif


   if (MSd.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Sd);
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
   problems.flag_bad_mass(MRSSMtower_info::Sv, eigenvalue_error >
      precision * Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif


   if (MSv.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Sv);
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
   problems.flag_bad_mass(MRSSMtower_info::Su, eigenvalue_error >
      precision * Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif


   if (MSu.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Su);
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
   problems.flag_bad_mass(MRSSMtower_info::Se, eigenvalue_error >
      precision * Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif


   if (MSe.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Se);
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
   problems.flag_bad_mass(MRSSMtower_info::hh, eigenvalue_error >
      precision * Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif


   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::hh);
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
   problems.flag_bad_mass(MRSSMtower_info::Ah, eigenvalue_error >
      precision * Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif


   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Ah);
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
   problems.flag_bad_mass(MRSSMtower_info::Rh, eigenvalue_error >
      precision * Abs(MRh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Rh, MRh, ZHR);
#endif


   if (MRh.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Rh);
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
   problems.flag_bad_mass(MRSSMtower_info::Hpm, eigenvalue_error >
      precision * Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif


   if (MHpm.minCoeff() < 0.) {
      problems.flag_tachyon(MRSSMtower_info::Hpm);
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
   problems.flag_bad_mass(MRSSMtower_info::Chi, eigenvalue_error >
      precision * Abs(MChi(0)));
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
   problems.flag_bad_mass(MRSSMtower_info::Cha1, eigenvalue_error >
      precision * Abs(MCha1(0)));
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
   problems.flag_bad_mass(MRSSMtower_info::Cha2, eigenvalue_error >
      precision * Abs(MCha2(0)));
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
   problems.flag_bad_mass(MRSSMtower_info::Fe, eigenvalue_error >
      precision * Abs(MFe(0)));
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
   problems.flag_bad_mass(MRSSMtower_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
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
   problems.flag_bad_mass(MRSSMtower_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
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
      problems.flag_tachyon(MRSSMtower_info::VWm);
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
   double result = Re(mHd2*vd - 0.3872983346207417*g1*MDBS*vd*vS + 0.5*g2*MDWBT
      *vd*vT + vd*AbsSqr(MuD) + vd*AbsSqr(Mu) - 0.5*vu*BMu + 0.7071067811865475*
      MuD*vd*vS*Conj(LamSD) + 0.35355339059327373*LamTD*vd*vS*vT*Conj(LamSD) + 0.5
      *MuD*vd*vT*Conj(LamTD) + 0.35355339059327373*LamSD*vd*vS*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*vS*Conj(MDBS) + 0.5*g2*vd*vT*Conj(MDWBT) +
      0.7071067811865475*LamSD*vd*vS*Conj(MuD) + 0.5*LamTD*vd*vT*Conj(MuD) - 0.5*
      vu*Conj(BMu) + 0.075*Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) + 0.5*
      vd*AbsSqr(LamSD)*Sqr(vS) + 0.25*vd*AbsSqr(LamTD)*Sqr(vT) - 0.075*vd*Sqr(g1)*
      Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu));

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
      vd*Conj(BMu) + 0.075*Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) - 0.075
      *vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(LamSU)*Sqr(vS
      ) + 0.25*vu*AbsSqr(LamTU)*Sqr(vT));

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

std::complex<double> CLASSNAME::CpUSdconjUSdconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2;
   std::complex<double> tmp_3;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_2 += tmp_3;
   result += (-0.1*Sqr(g1)) * tmp_2;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_4 += tmp_5;
   result += (0.1*Sqr(g1)) * tmp_4;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
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

std::complex<double> CLASSNAME::CpUSdconjUSdconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_6;
   std::complex<double> tmp_7;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_6 += tmp_7;
   result += (-0.1*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0)) * tmp_6;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   std::complex<double> tmp_8;
   std::complex<double> tmp_9;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_9 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_8 += tmp_9;
   result += (0.1*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1)) * tmp_8;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1);
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

std::complex<double> CLASSNAME::CpconjUSdSdphiO(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_62;
   std::complex<double> tmp_63;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_63 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_62 += tmp_63;
   result += (g3*MDGoc) * tmp_62;
   std::complex<double> tmp_64;
   std::complex<double> tmp_65;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_65 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_64 += tmp_65;
   result += (g3*Conj(MDGoc)) * tmp_64;
   if (gO2 < 3) {
      result += -(g3*MDGoc*Conj(ZD(gI1,gO2)));
   }
   if (gO2 < 3) {
      result += -(g3*Conj(MDGoc)*Conj(ZD(gI1,gO2)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuSRum(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_66;
      std::complex<double> tmp_67;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_67 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_66 += tmp_67;
      result += (-MuU) * tmp_66;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_69 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_68 += tmp_69;
      result += (-0.7071067811865475*LamSU*vS) * tmp_68;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_71 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_70 += tmp_71;
      result += (-0.5*LamTU*vT) * tmp_70;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_72;
   std::complex<double> tmp_74;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_74 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_72 += tmp_74;
   std::complex<double> tmp_73;
   std::complex<double> tmp_75;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_75 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_73 += tmp_75;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_72 * tmp_73;
   std::complex<double> tmp_76;
   std::complex<double> tmp_78;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_78 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_76 += tmp_78;
   std::complex<double> tmp_77;
   std::complex<double> tmp_79;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_79 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_77 += tmp_79;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_76 * tmp_77;
   std::complex<double> tmp_80;
   std::complex<double> tmp_82;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_82 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_80 += tmp_82;
   std::complex<double> tmp_81;
   std::complex<double> tmp_83;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_83 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_81 += tmp_83;
   result += (-0.05*Sqr(g1)) * tmp_80 * tmp_81;
   std::complex<double> tmp_84;
   std::complex<double> tmp_86;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_86 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_84 += tmp_86;
   std::complex<double> tmp_85;
   std::complex<double> tmp_87;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_87 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_85 += tmp_87;
   result += (-0.1*Sqr(g1)) * tmp_84 * tmp_85;
   std::complex<double> tmp_88;
   std::complex<double> tmp_90;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_90 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_88 += tmp_90;
   std::complex<double> tmp_89;
   std::complex<double> tmp_91;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_91 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_89 += tmp_91;
   result += (-0.05*Sqr(g1)) * tmp_88 * tmp_89;
   std::complex<double> tmp_92;
   std::complex<double> tmp_94;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_94 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_92 += tmp_94;
   std::complex<double> tmp_93;
   std::complex<double> tmp_95;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_95 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_93 += tmp_95;
   result += (-0.1*Sqr(g1)) * tmp_92 * tmp_93;
   std::complex<double> tmp_96;
   std::complex<double> tmp_98;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_98 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_96 += tmp_98;
   std::complex<double> tmp_97;
   std::complex<double> tmp_99;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_99 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_97 += tmp_99;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_96 * tmp_97;
   std::complex<double> tmp_100;
   std::complex<double> tmp_102;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_102 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_100 += tmp_102;
   std::complex<double> tmp_101;
   std::complex<double> tmp_103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_103 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_101 += tmp_103;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_100 * tmp_101;
   std::complex<double> tmp_104;
   std::complex<double> tmp_106;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_107;
      std::complex<double> tmp_108;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_108 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_107 += tmp_108;
      tmp_106 += (Conj(ZD(gI2,j2))) * tmp_107;
   }
   tmp_104 += tmp_106;
   std::complex<double> tmp_105;
   std::complex<double> tmp_109;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_111 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_110 += tmp_111;
      tmp_109 += (ZD(gI1,j4)) * tmp_110;
   }
   tmp_105 += tmp_109;
   result += (-1) * tmp_104 * tmp_105;
   if (gO1 < 3) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_114;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_114 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_112 += tmp_114;
      std::complex<double> tmp_113;
      std::complex<double> tmp_115;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_116;
         std::complex<double> tmp_117;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_117 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_116 += tmp_117;
         tmp_115 += (ZD(gI1,j4)) * tmp_116;
      }
      tmp_113 += tmp_115;
      result += (-3) * tmp_112 * tmp_113;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_118;
      std::complex<double> tmp_119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_119 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_118 += tmp_119;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_118;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_121 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_120 += tmp_121;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_120;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_123;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_123 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_122 += tmp_123;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_122;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_125;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_125 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_124 += tmp_125;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_124;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_126;
      std::complex<double> tmp_127;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_127 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_126 += tmp_127;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_126;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_129 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_128 += tmp_129;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_128;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_130;
      std::complex<double> tmp_131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_131 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_130 += tmp_131;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_130;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_133 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_132 += tmp_133;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_132;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_135 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_134 += tmp_135;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_134;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_136;
      std::complex<double> tmp_137;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_137 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_136 += tmp_137;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_136;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_138;
      std::complex<double> tmp_140;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_141;
         std::complex<double> tmp_142;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_142 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_141 += tmp_142;
         tmp_140 += (Conj(ZD(gI2,j2))) * tmp_141;
      }
      tmp_138 += tmp_140;
      std::complex<double> tmp_139;
      std::complex<double> tmp_143;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_143 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_139 += tmp_143;
      result += (-3) * tmp_138 * tmp_139;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_144;
      std::complex<double> tmp_146;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_146 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_144 += tmp_146;
      std::complex<double> tmp_145;
      std::complex<double> tmp_147;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_147 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_145 += tmp_147;
      result += (-1) * tmp_144 * tmp_145;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_148;
      std::complex<double> tmp_149;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_149 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_148 += tmp_149;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_148;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_150;
      std::complex<double> tmp_151;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_151 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_150 += tmp_151;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_150;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_153 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_152 += tmp_153;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_152;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_154;
      std::complex<double> tmp_155;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_155 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_154 += tmp_155;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_154;
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

   std::complex<double> tmp_156;
   std::complex<double> tmp_158;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_158 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_156 += tmp_158;
   std::complex<double> tmp_157;
   std::complex<double> tmp_159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_159 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_157 += tmp_159;
   result += (0.05*Sqr(g1)) * tmp_156 * tmp_157;
   std::complex<double> tmp_160;
   std::complex<double> tmp_162;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_162 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_160 += tmp_162;
   std::complex<double> tmp_161;
   std::complex<double> tmp_163;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_163 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_161 += tmp_163;
   result += (-0.1*Sqr(g1)) * tmp_160 * tmp_161;
   std::complex<double> tmp_164;
   std::complex<double> tmp_166;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_166 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_164 += tmp_166;
   std::complex<double> tmp_165;
   std::complex<double> tmp_167;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_167 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_165 += tmp_167;
   result += (0.05*Sqr(g1)) * tmp_164 * tmp_165;
   std::complex<double> tmp_168;
   std::complex<double> tmp_170;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_170 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_168 += tmp_170;
   std::complex<double> tmp_169;
   std::complex<double> tmp_171;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_171 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_169 += tmp_171;
   result += (-0.1*Sqr(g1)) * tmp_168 * tmp_169;
   if (gO1 < 3) {
      std::complex<double> tmp_172;
      std::complex<double> tmp_174;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_174 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_172 += tmp_174;
      std::complex<double> tmp_173;
      std::complex<double> tmp_175;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_176;
         std::complex<double> tmp_177;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_177 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_176 += tmp_177;
         tmp_175 += (ZE(gI1,j4)) * tmp_176;
      }
      tmp_173 += tmp_175;
      result += (-1) * tmp_172 * tmp_173;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_178;
      std::complex<double> tmp_179;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_179 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_178 += tmp_179;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_181;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_181 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_180 += tmp_181;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_180;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_182;
      std::complex<double> tmp_183;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_183 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_182 += tmp_183;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_182;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_184;
      std::complex<double> tmp_185;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_185 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_184 += tmp_185;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_184;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_186;
      std::complex<double> tmp_187;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_187 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_186 += tmp_187;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_186;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_188;
      std::complex<double> tmp_189;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_189 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_188 += tmp_189;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_188;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_190;
      std::complex<double> tmp_192;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_193;
         std::complex<double> tmp_194;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_194 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_193 += tmp_194;
         tmp_192 += (Conj(ZE(gI2,j2))) * tmp_193;
      }
      tmp_190 += tmp_192;
      std::complex<double> tmp_191;
      std::complex<double> tmp_195;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_195 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_191 += tmp_195;
      result += (-1) * tmp_190 * tmp_191;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_196;
   std::complex<double> tmp_198;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_198 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_196 += tmp_198;
   std::complex<double> tmp_197;
   std::complex<double> tmp_199;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_199 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_197 += tmp_199;
   result += (-0.05*Sqr(g1)) * tmp_196 * tmp_197;
   std::complex<double> tmp_200;
   std::complex<double> tmp_202;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_202 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_200 += tmp_202;
   std::complex<double> tmp_201;
   std::complex<double> tmp_203;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_203 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_201 += tmp_203;
   result += (0.2*Sqr(g1)) * tmp_200 * tmp_201;
   std::complex<double> tmp_204;
   std::complex<double> tmp_206;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_206 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_204 += tmp_206;
   std::complex<double> tmp_205;
   std::complex<double> tmp_207;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_207 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_205 += tmp_207;
   result += (-0.05*Sqr(g1)) * tmp_204 * tmp_205;
   std::complex<double> tmp_208;
   std::complex<double> tmp_210;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_210 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_208 += tmp_210;
   std::complex<double> tmp_209;
   std::complex<double> tmp_211;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_211 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_209 += tmp_211;
   result += (0.2*Sqr(g1)) * tmp_208 * tmp_209;
   std::complex<double> tmp_212;
   std::complex<double> tmp_214;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_215;
      std::complex<double> tmp_216;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_216 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_215 += tmp_216;
      tmp_214 += (Conj(ZU(gI2,j2))) * tmp_215;
   }
   tmp_212 += tmp_214;
   std::complex<double> tmp_213;
   std::complex<double> tmp_217;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_218;
      std::complex<double> tmp_219;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_219 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_218 += tmp_219;
      tmp_217 += (ZU(gI1,j4)) * tmp_218;
   }
   tmp_213 += tmp_217;
   result += (-1) * tmp_212 * tmp_213;
   if (gO1 < 3) {
      std::complex<double> tmp_220;
      std::complex<double> tmp_221;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_221 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_220 += tmp_221;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_220;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_222;
      std::complex<double> tmp_223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_223 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_222 += tmp_223;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_222;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_224;
      std::complex<double> tmp_225;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_225 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_224 += tmp_225;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_224;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_226;
      std::complex<double> tmp_227;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_227 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_226 += tmp_227;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_226;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_228;
      std::complex<double> tmp_229;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_229 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_228 += tmp_229;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_228;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_230;
      std::complex<double> tmp_231;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_231 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_230 += tmp_231;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_230;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_232;
      std::complex<double> tmp_234;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_234 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_232 += tmp_234;
      std::complex<double> tmp_233;
      std::complex<double> tmp_235;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_235 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_233 += tmp_235;
      result += (-1) * tmp_232 * tmp_233;
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
      std::complex<double> tmp_236;
      std::complex<double> tmp_237;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_237 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_236 += tmp_237;
      result += (MuD*ZHR(gI2,0)) * tmp_236;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_238;
      std::complex<double> tmp_239;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_239 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_238 += tmp_239;
      result += (0.7071067811865475*LamSD*vS*ZHR(gI2,0)) * tmp_238;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_240;
      std::complex<double> tmp_241;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_241 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_240 += tmp_241;
      result += (0.5*LamTD*vT*ZHR(gI2,0)) * tmp_240;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_242;
   std::complex<double> tmp_243;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_244;
      std::complex<double> tmp_245;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_245 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_244 += tmp_245;
      tmp_243 += (Conj(ZD(gI1,j2))) * tmp_244;
   }
   tmp_242 += tmp_243;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1))
      * tmp_242;
   if (gO2 < 3) {
      std::complex<double> tmp_246;
      std::complex<double> tmp_247;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_247 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_246 += tmp_247;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_246;
   }
   std::complex<double> tmp_248;
   std::complex<double> tmp_249;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_249 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_248 += tmp_249;
   result += (std::complex<double>(0.,-0.2581988897471611)*g1*MDBS*ZA(gI2,2)) *
      tmp_248;
   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_251 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_250 += tmp_251;
   result += (std::complex<double>(0.,0.2581988897471611)*g1*Conj(MDBS)*ZA(gI2,
      2)) * tmp_250;
   if (gO2 < 3) {
      result += std::complex<double>(0.,-0.12909944487358055)*g1*MDBS*Conj(
         ZD(gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0.,0.12909944487358055)*g1*Conj(MDBS)*
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

   std::complex<double> tmp_252;
   std::complex<double> tmp_253;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_253 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_252 += tmp_253;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_252;
   std::complex<double> tmp_254;
   std::complex<double> tmp_255;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_258;
         std::complex<double> tmp_259;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_259 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_258 += tmp_259;
         tmp_257 += (KroneckerDelta(gO2,3 + j2)) * tmp_258;
      }
      tmp_256 += tmp_257;
      tmp_255 += (Conj(ZD(gI1,3 + j3))) * tmp_256;
   }
   tmp_254 += tmp_255;
   result += (-(vd*ZH(gI2,0))) * tmp_254;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_260;
      std::complex<double> tmp_261;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_262;
         std::complex<double> tmp_263;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_263 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_262 += tmp_263;
         tmp_261 += (Conj(ZD(gI1,j2))) * tmp_262;
      }
      tmp_260 += tmp_261;
      result += (-(vd*ZH(gI2,0))) * tmp_260;
   }
   std::complex<double> tmp_264;
   std::complex<double> tmp_265;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_265 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_264 += tmp_265;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_264;
   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_269 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_268 += tmp_269;
      tmp_267 += (Conj(ZD(gI1,j2))) * tmp_268;
   }
   tmp_266 += tmp_267;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_266;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_270;
      std::complex<double> tmp_271;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_271 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_270 += tmp_271;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_270;
   }
   std::complex<double> tmp_272;
   std::complex<double> tmp_273;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_273 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_272 += tmp_273;
   result += (-0.2581988897471611*g1*MDBS*ZH(gI2,2)) * tmp_272;
   std::complex<double> tmp_274;
   std::complex<double> tmp_275;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_275 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_274 += tmp_275;
   result += (-0.2581988897471611*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_274;
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

   std::complex<double> tmp_276;
   std::complex<double> tmp_277;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_278;
      std::complex<double> tmp_279;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_280;
         std::complex<double> tmp_281;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_281 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_280 += tmp_281;
         tmp_279 += (KroneckerDelta(gO2,3 + j2)) * tmp_280;
      }
      tmp_278 += tmp_279;
      tmp_277 += (Conj(ZU(gI1,3 + j3))) * tmp_278;
   }
   tmp_276 += tmp_277;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_276;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_282;
      std::complex<double> tmp_283;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_283 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_282 += tmp_283;
      result += (Mu*ZP(gI2,0)) * tmp_282;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_284;
      std::complex<double> tmp_285;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_286;
         std::complex<double> tmp_287;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_287 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_286 += tmp_287;
         tmp_285 += (Conj(ZU(gI1,j2))) * tmp_286;
      }
      tmp_284 += tmp_285;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_284;
   }
   std::complex<double> tmp_288;
   std::complex<double> tmp_289;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_290;
      std::complex<double> tmp_291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_291 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_290 += tmp_291;
      tmp_289 += (Conj(ZU(gI1,j2))) * tmp_290;
   }
   tmp_288 += tmp_289;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_288;
   std::complex<double> tmp_292;
   std::complex<double> tmp_293;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_294;
      std::complex<double> tmp_295;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_296;
         std::complex<double> tmp_297;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_297 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_296 += tmp_297;
         tmp_295 += (KroneckerDelta(gO2,3 + j2)) * tmp_296;
      }
      tmp_294 += tmp_295;
      tmp_293 += (Conj(ZU(gI1,3 + j3))) * tmp_294;
   }
   tmp_292 += tmp_293;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_292;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_298;
      std::complex<double> tmp_299;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_300;
         std::complex<double> tmp_301;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_301 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_300 += tmp_301;
         tmp_299 += (Conj(ZU(gI1,j2))) * tmp_300;
      }
      tmp_298 += tmp_299;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_298;
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

   std::complex<double> tmp_302;
   std::complex<double> tmp_303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_303 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_302 += tmp_303;
   result += (1.4142135623730951*g3) * tmp_302;

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

std::complex<double> CLASSNAME::CpconjUSdsigmaOSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_305 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_304 += tmp_305;
   result += (std::complex<double>(0,1)*g3*MDGoc) * tmp_304;
   std::complex<double> tmp_306;
   std::complex<double> tmp_307;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_307 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_306 += tmp_307;
   result += (std::complex<double>(0,-1)*g3*Conj(MDGoc)) * tmp_306;
   if (gO2 < 3) {
      result += std::complex<double>(0,-1)*g3*MDGoc*Conj(ZD(gI2,gO2));
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,1)*g3*Conj(MDGoc)*Conj(ZD(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdconjSRdpSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_308;
   std::complex<double> tmp_309;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_310;
      std::complex<double> tmp_311;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_311 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_310 += tmp_311;
      tmp_309 += (Conj(ZU(gI2,j2))) * tmp_310;
   }
   tmp_308 += tmp_309;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(
      MuD))) * tmp_308;

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

   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_313 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_312 += tmp_313;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_312;
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

   std::complex<double> tmp_314;
   std::complex<double> tmp_315;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_315 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_314 += tmp_315;
   result += (0.2581988897471611*g1*Sin(ThetaW())) * tmp_314;
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

std::complex<double> CLASSNAME::CpUSvconjUSvconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*KroneckerDelta(gO1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2));

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
      std::complex<double> tmp_316;
      std::complex<double> tmp_317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_317 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_316 += tmp_317;
      result += (UM1(gI2,1)) * tmp_316;
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
      result += std::complex<double>(0.,0.3872983346207417)*g1*MDBS*Conj(ZV(
         gI1,gO2))*ZA(gI2,2);
   }
   if (gI1 < 3) {
      result += std::complex<double>(0.,-0.3872983346207417)*g1*Conj(MDBS)*
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
      std::complex<double> tmp_318;
      std::complex<double> tmp_319;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_319 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_318 += tmp_319;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_318;
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
      std::complex<double> tmp_320;
      std::complex<double> tmp_321;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_322;
         std::complex<double> tmp_323;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_323 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_322 += tmp_323;
         tmp_321 += (Conj(ZE(gI2,j2))) * tmp_322;
      }
      tmp_320 += tmp_321;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_320;
   }
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_324 += tmp_325;
      result += (Mu*ZP(gI1,1)) * tmp_324;
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

   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   std::complex<double> tmp_328;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_328 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_327 += tmp_328;
   tmp_326 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_327;
   std::complex<double> tmp_329;
   std::complex<double> tmp_330;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_330 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_329 += tmp_330;
   tmp_326 += (std::complex<double>(0,0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_329;
   std::complex<double> tmp_331;
   std::complex<double> tmp_332;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_332 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_331 += tmp_332;
   tmp_326 += (std::complex<double>(0,0.1)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_331;
   result += (std::complex<double>(0,-1)) * tmp_326;

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_333;
   std::complex<double> tmp_334;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_334 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_333 += tmp_334;
   result += (-0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_333;
   std::complex<double> tmp_335;
   std::complex<double> tmp_336;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_336 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_335 += tmp_336;
   result += (0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_335;
   std::complex<double> tmp_337;
   std::complex<double> tmp_338;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_338 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_337 += tmp_338;
   result += (0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_337;
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_339;
      std::complex<double> tmp_341;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_341 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_339 += tmp_341;
      std::complex<double> tmp_340;
      std::complex<double> tmp_342;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_342 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_340 += tmp_342;
      result += (-1) * tmp_339 * tmp_340;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_343;
   std::complex<double> tmp_344;
   std::complex<double> tmp_345;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_345 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_344 += tmp_345;
   tmp_343 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_344;
   std::complex<double> tmp_346;
   std::complex<double> tmp_347;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_347 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_346 += tmp_347;
   tmp_343 += (std::complex<double>(0,-0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_346;
   std::complex<double> tmp_348;
   std::complex<double> tmp_349;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_349 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_348 += tmp_349;
   tmp_343 += (std::complex<double>(0,-0.2)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_348;
   result += (std::complex<double>(0,-1)) * tmp_343;

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

std::complex<double> CLASSNAME::CpconjUSvSRdpSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_350;
      std::complex<double> tmp_351;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_351 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_350 += tmp_351;
      result += (MuD) * tmp_350;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_352;
      std::complex<double> tmp_353;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_353 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_352 += tmp_353;
      result += (0.7071067811865475*LamSD*vS) * tmp_352;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_354;
      std::complex<double> tmp_355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_355 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_354 += tmp_355;
      result += (-0.5*LamTD*vT) * tmp_354;
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

   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_357 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_356 += tmp_357;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_356;
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

std::complex<double> CLASSNAME::CpUSuconjUSuconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_359 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_358 += tmp_359;
   result += (0.2*Sqr(g1)) * tmp_358;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_361 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_360 += tmp_361;
   result += (-0.2*Sqr(g1)) * tmp_360;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
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

std::complex<double> CLASSNAME::CpUSuconjUSuconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_362;
   std::complex<double> tmp_363;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_363 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_362 += tmp_363;
   result += (0.2*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0)) * tmp_362;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   std::complex<double> tmp_364;
   std::complex<double> tmp_365;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_365 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_364 += tmp_365;
   result += (-0.2*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1)) * tmp_364;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1);
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

   std::complex<double> tmp_366;
   std::complex<double> tmp_367;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_368;
      std::complex<double> tmp_369;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_369 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_368 += tmp_369;
      tmp_367 += (Conj(ZDL(gI2,j2))) * tmp_368;
   }
   tmp_366 += tmp_367;
   result += (Conj(UP2(gI1,1))) * tmp_366;

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_370;
   std::complex<double> tmp_371;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_371 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_370 += tmp_371;
   result += (-0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_370;
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
      std::complex<double> tmp_372;
      std::complex<double> tmp_373;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_373 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_372 += tmp_373;
      result += (UM1(gI2,1)) * tmp_372;
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
      std::complex<double> tmp_374;
      std::complex<double> tmp_375;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_375 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_374 += tmp_375;
      result += (-ZN2(gI2,3)) * tmp_374;
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

   std::complex<double> tmp_376;
   std::complex<double> tmp_377;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_377 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_376 += tmp_377;
   result += (-0.2*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_376;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   std::complex<double> tmp_378;
   std::complex<double> tmp_379;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_379 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_378 += tmp_379;
   result += (0.2*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_378;
   std::complex<double> tmp_380;
   std::complex<double> tmp_381;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_382;
      std::complex<double> tmp_383;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_384;
         std::complex<double> tmp_385;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_385 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_384 += tmp_385;
         tmp_383 += (KroneckerDelta(gO2,3 + j2)) * tmp_384;
      }
      tmp_382 += tmp_383;
      tmp_381 += (KroneckerDelta(gO1,3 + j3)) * tmp_382;
   }
   tmp_380 += tmp_381;
   result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_380;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_386;
      std::complex<double> tmp_387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_387 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_386 += tmp_387;
      result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_386;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_389 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_388 += tmp_389;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_388;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_390;
      std::complex<double> tmp_391;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_391 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_390 += tmp_391;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_390;
   }
   std::complex<double> tmp_392;
   std::complex<double> tmp_393;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_393 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_392 += tmp_393;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_392;
   std::complex<double> tmp_394;
   std::complex<double> tmp_395;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_396;
      std::complex<double> tmp_397;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_398;
         std::complex<double> tmp_399;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_399 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_398 += tmp_399;
         tmp_397 += (KroneckerDelta(gO2,3 + j2)) * tmp_398;
      }
      tmp_396 += tmp_397;
      tmp_395 += (KroneckerDelta(gO1,3 + j3)) * tmp_396;
   }
   tmp_394 += tmp_395;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_394;
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

   std::complex<double> tmp_400;
   std::complex<double> tmp_401;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_401 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_400 += tmp_401;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_400;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   std::complex<double> tmp_402;
   std::complex<double> tmp_403;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_403 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_402 += tmp_403;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_402;
   std::complex<double> tmp_404;
   std::complex<double> tmp_405;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_406;
      std::complex<double> tmp_407;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_408;
         std::complex<double> tmp_409;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_409 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_408 += tmp_409;
         tmp_407 += (KroneckerDelta(gO2,3 + j2)) * tmp_408;
      }
      tmp_406 += tmp_407;
      tmp_405 += (KroneckerDelta(gO1,3 + j3)) * tmp_406;
   }
   tmp_404 += tmp_405;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_404;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_410;
      std::complex<double> tmp_411;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_411 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_410 += tmp_411;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_410;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChiFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_413 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_412 += tmp_413;
   result += (0.7302967433402214*g1*ZN1(gI1,0)) * tmp_412;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChiFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_416;
      std::complex<double> tmp_417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_417 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_416 += tmp_417;
      tmp_415 += (Conj(ZUL(gI2,j2))) * tmp_416;
   }
   tmp_414 += tmp_415;
   result += (-Conj(ZN2(gI1,3))) * tmp_414;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_420;
      std::complex<double> tmp_421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_421 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_420 += tmp_421;
      tmp_419 += (Conj(ZD(gI2,j2))) * tmp_420;
   }
   tmp_418 += tmp_419;
   result += (Conj(Mu)*ZP(gI1,0)) * tmp_418;
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_424;
      std::complex<double> tmp_425;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_426;
         std::complex<double> tmp_427;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_427 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_426 += tmp_427;
         tmp_425 += (KroneckerDelta(gO2,3 + j2)) * tmp_426;
      }
      tmp_424 += tmp_425;
      tmp_423 += (Conj(ZD(gI2,3 + j3))) * tmp_424;
   }
   tmp_422 += tmp_423;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_422;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_429;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_430;
         std::complex<double> tmp_431;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_431 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_430 += tmp_431;
         tmp_429 += (Conj(ZD(gI2,j2))) * tmp_430;
      }
      tmp_428 += tmp_429;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_428;
   }
   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_434;
      std::complex<double> tmp_435;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_436;
         std::complex<double> tmp_437;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_437 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_436 += tmp_437;
         tmp_435 += (KroneckerDelta(gO2,3 + j2)) * tmp_436;
      }
      tmp_434 += tmp_435;
      tmp_433 += (Conj(ZD(gI2,3 + j3))) * tmp_434;
   }
   tmp_432 += tmp_433;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_432;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_438;
      std::complex<double> tmp_439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_439 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_438 += tmp_439;
      result += (Mu*ZP(gI1,1)) * tmp_438;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_440;
      std::complex<double> tmp_441;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_442;
         std::complex<double> tmp_443;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_443 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_442 += tmp_443;
         tmp_441 += (Conj(ZD(gI2,j2))) * tmp_442;
      }
      tmp_440 += tmp_441;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_440;
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

std::complex<double> CLASSNAME::CpconjUSuSuphiO(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_444;
   std::complex<double> tmp_445;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_445 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_444 += tmp_445;
   result += (g3*MDGoc) * tmp_444;
   std::complex<double> tmp_446;
   std::complex<double> tmp_447;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_447 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_446 += tmp_447;
   result += (g3*Conj(MDGoc)) * tmp_446;
   if (gO2 < 3) {
      result += -(g3*MDGoc*Conj(ZU(gI1,gO2)));
   }
   if (gO2 < 3) {
      result += -(g3*Conj(MDGoc)*Conj(ZU(gI1,gO2)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSusigmaO(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_448;
   std::complex<double> tmp_449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_449 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_448 += tmp_449;
   result += (std::complex<double>(0,1)*g3*MDGoc) * tmp_448;
   std::complex<double> tmp_450;
   std::complex<double> tmp_451;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_451 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_450 += tmp_451;
   result += (std::complex<double>(0,-1)*g3*Conj(MDGoc)) * tmp_450;
   if (gO2 < 3) {
      result += std::complex<double>(0,-1)*g3*MDGoc*Conj(ZU(gI1,gO2));
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,1)*g3*Conj(MDGoc)*Conj(ZU(gI1,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_452;
   std::complex<double> tmp_454;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_454 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_452 += tmp_454;
   std::complex<double> tmp_453;
   std::complex<double> tmp_455;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_455 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_453 += tmp_455;
   result += (0.1*Sqr(g1)) * tmp_452 * tmp_453;
   std::complex<double> tmp_456;
   std::complex<double> tmp_458;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_458 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_456 += tmp_458;
   std::complex<double> tmp_457;
   std::complex<double> tmp_459;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_459 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_457 += tmp_459;
   result += (0.2*Sqr(g1)) * tmp_456 * tmp_457;
   std::complex<double> tmp_460;
   std::complex<double> tmp_462;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_462 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_460 += tmp_462;
   std::complex<double> tmp_461;
   std::complex<double> tmp_463;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_463 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_461 += tmp_463;
   result += (0.1*Sqr(g1)) * tmp_460 * tmp_461;
   std::complex<double> tmp_464;
   std::complex<double> tmp_466;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_466 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_464 += tmp_466;
   std::complex<double> tmp_465;
   std::complex<double> tmp_467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_467 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_465 += tmp_467;
   result += (0.2*Sqr(g1)) * tmp_464 * tmp_465;
   std::complex<double> tmp_468;
   std::complex<double> tmp_470;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_471;
      std::complex<double> tmp_472;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_472 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_471 += tmp_472;
      tmp_470 += (Conj(ZD(gI2,j2))) * tmp_471;
   }
   tmp_468 += tmp_470;
   std::complex<double> tmp_469;
   std::complex<double> tmp_473;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_474;
      std::complex<double> tmp_475;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_475 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_474 += tmp_475;
      tmp_473 += (ZD(gI1,j4)) * tmp_474;
   }
   tmp_469 += tmp_473;
   result += (-1) * tmp_468 * tmp_469;
   if (gO1 < 3) {
      std::complex<double> tmp_476;
      std::complex<double> tmp_477;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_477 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_476 += tmp_477;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_476;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_478;
      std::complex<double> tmp_479;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_479 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_478 += tmp_479;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_478;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_480;
      std::complex<double> tmp_481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_481 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_480 += tmp_481;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_480;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_482;
      std::complex<double> tmp_483;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_483 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_482 += tmp_483;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_482;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_484;
      std::complex<double> tmp_485;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_485 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_484 += tmp_485;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_484;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_486;
      std::complex<double> tmp_487;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_487 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_486 += tmp_487;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_486;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_488;
      std::complex<double> tmp_490;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_490 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_488 += tmp_490;
      std::complex<double> tmp_489;
      std::complex<double> tmp_491;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_491 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_489 += tmp_491;
      result += (-1) * tmp_488 * tmp_489;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_492;
   std::complex<double> tmp_494;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_494 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_492 += tmp_494;
   std::complex<double> tmp_493;
   std::complex<double> tmp_495;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_495 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_493 += tmp_495;
   result += (-0.1*Sqr(g1)) * tmp_492 * tmp_493;
   std::complex<double> tmp_496;
   std::complex<double> tmp_498;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_498 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_496 += tmp_498;
   std::complex<double> tmp_497;
   std::complex<double> tmp_499;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_499 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_497 += tmp_499;
   result += (0.2*Sqr(g1)) * tmp_496 * tmp_497;
   std::complex<double> tmp_500;
   std::complex<double> tmp_502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_502 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_500 += tmp_502;
   std::complex<double> tmp_501;
   std::complex<double> tmp_503;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_503 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_501 += tmp_503;
   result += (-0.1*Sqr(g1)) * tmp_500 * tmp_501;
   std::complex<double> tmp_504;
   std::complex<double> tmp_506;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_506 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_504 += tmp_506;
   std::complex<double> tmp_505;
   std::complex<double> tmp_507;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_507 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_505 += tmp_507;
   result += (0.2*Sqr(g1)) * tmp_504 * tmp_505;
   if (gO1 < 3) {
      std::complex<double> tmp_508;
      std::complex<double> tmp_509;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_509 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_508 += tmp_509;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_508;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_510;
      std::complex<double> tmp_511;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_511 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_510 += tmp_511;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_510;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_512;
      std::complex<double> tmp_513;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_513 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_512 += tmp_513;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_512;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_514;
      std::complex<double> tmp_515;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_515 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_514 += tmp_515;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_514;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_516;
      std::complex<double> tmp_517;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_517 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_516 += tmp_517;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_516;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_518;
      std::complex<double> tmp_519;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_519 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_518 += tmp_519;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_518;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_520;
   std::complex<double> tmp_522;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_522 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_520 += tmp_522;
   std::complex<double> tmp_521;
   std::complex<double> tmp_523;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_523 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_521 += tmp_523;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_520 * tmp_521;
   std::complex<double> tmp_524;
   std::complex<double> tmp_526;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_526 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_524 += tmp_526;
   std::complex<double> tmp_525;
   std::complex<double> tmp_527;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_527 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_525 += tmp_527;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_524 * tmp_525;
   std::complex<double> tmp_528;
   std::complex<double> tmp_530;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_530 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_528 += tmp_530;
   std::complex<double> tmp_529;
   std::complex<double> tmp_531;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_531 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_529 += tmp_531;
   result += (0.1*Sqr(g1)) * tmp_528 * tmp_529;
   std::complex<double> tmp_532;
   std::complex<double> tmp_534;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_534 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_532 += tmp_534;
   std::complex<double> tmp_533;
   std::complex<double> tmp_535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_535 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_533 += tmp_535;
   result += (-0.4*Sqr(g1)) * tmp_532 * tmp_533;
   std::complex<double> tmp_536;
   std::complex<double> tmp_538;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_538 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_536 += tmp_538;
   std::complex<double> tmp_537;
   std::complex<double> tmp_539;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_539 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_537 += tmp_539;
   result += (0.1*Sqr(g1)) * tmp_536 * tmp_537;
   std::complex<double> tmp_540;
   std::complex<double> tmp_542;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_542 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_540 += tmp_542;
   std::complex<double> tmp_541;
   std::complex<double> tmp_543;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_543 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_541 += tmp_543;
   result += (-0.4*Sqr(g1)) * tmp_540 * tmp_541;
   std::complex<double> tmp_544;
   std::complex<double> tmp_546;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_546 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_544 += tmp_546;
   std::complex<double> tmp_545;
   std::complex<double> tmp_547;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_547 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_545 += tmp_547;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_544 * tmp_545;
   std::complex<double> tmp_548;
   std::complex<double> tmp_550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_550 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_548 += tmp_550;
   std::complex<double> tmp_549;
   std::complex<double> tmp_551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_551 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_549 += tmp_551;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_548 * tmp_549;
   std::complex<double> tmp_552;
   std::complex<double> tmp_554;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_555;
      std::complex<double> tmp_556;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_556 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_555 += tmp_556;
      tmp_554 += (Conj(ZU(gI2,j2))) * tmp_555;
   }
   tmp_552 += tmp_554;
   std::complex<double> tmp_553;
   std::complex<double> tmp_557;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_558;
      std::complex<double> tmp_559;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_559 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_558 += tmp_559;
      tmp_557 += (ZU(gI1,j4)) * tmp_558;
   }
   tmp_553 += tmp_557;
   result += (-1) * tmp_552 * tmp_553;
   if (gO1 < 3) {
      std::complex<double> tmp_560;
      std::complex<double> tmp_562;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_562 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_560 += tmp_562;
      std::complex<double> tmp_561;
      std::complex<double> tmp_563;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_564;
         std::complex<double> tmp_565;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_565 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_564 += tmp_565;
         tmp_563 += (ZU(gI1,j4)) * tmp_564;
      }
      tmp_561 += tmp_563;
      result += (-3) * tmp_560 * tmp_561;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_566;
      std::complex<double> tmp_567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_567 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_566 += tmp_567;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_566;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_568;
      std::complex<double> tmp_569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_569 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_568 += tmp_569;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_568;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_570;
      std::complex<double> tmp_571;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_571 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_570 += tmp_571;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_570;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_572;
      std::complex<double> tmp_573;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_573 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_572 += tmp_573;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_572;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_574;
      std::complex<double> tmp_575;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_575 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_574 += tmp_575;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_574;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_576;
      std::complex<double> tmp_577;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_577 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_576 += tmp_577;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_576;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_578;
      std::complex<double> tmp_579;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_579 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_578 += tmp_579;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_578;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_580;
      std::complex<double> tmp_581;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_581 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_580 += tmp_581;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_580;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_582;
      std::complex<double> tmp_583;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_583 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_582 += tmp_583;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_582;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_584;
      std::complex<double> tmp_585;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_585 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_584 += tmp_585;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_584;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_586;
      std::complex<double> tmp_588;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_589;
         std::complex<double> tmp_590;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_590 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_589 += tmp_590;
         tmp_588 += (Conj(ZU(gI2,j2))) * tmp_589;
      }
      tmp_586 += tmp_588;
      std::complex<double> tmp_587;
      std::complex<double> tmp_591;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_591 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_587 += tmp_591;
      result += (-3) * tmp_586 * tmp_587;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_592;
      std::complex<double> tmp_594;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_594 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_592 += tmp_594;
      std::complex<double> tmp_593;
      std::complex<double> tmp_595;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_595 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_593 += tmp_595;
      result += (-1) * tmp_592 * tmp_593;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_596;
      std::complex<double> tmp_597;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_597 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_596 += tmp_597;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_596;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_598;
      std::complex<double> tmp_599;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_599 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_598 += tmp_599;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_598;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_600;
      std::complex<double> tmp_601;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_601 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_600 += tmp_601;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_600;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_602;
      std::complex<double> tmp_603;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_603 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_602 += tmp_603;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_602;
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
      std::complex<double> tmp_604;
      std::complex<double> tmp_605;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_605 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_604 += tmp_605;
      result += (-(MuU*ZHR(gI2,1))) * tmp_604;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_606;
      std::complex<double> tmp_607;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_607 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_606 += tmp_607;
      result += (-0.7071067811865475*LamSU*vS*ZHR(gI2,1)) * tmp_606;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_608;
      std::complex<double> tmp_609;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_609 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_608 += tmp_609;
      result += (0.5*LamTU*vT*ZHR(gI2,1)) * tmp_608;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_610;
   std::complex<double> tmp_611;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_612;
      std::complex<double> tmp_613;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_613 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_612 += tmp_613;
      tmp_611 += (Conj(ZU(gI1,j2))) * tmp_612;
   }
   tmp_610 += tmp_611;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,0))
      * tmp_610;
   if (gO2 < 3) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_615 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_614 += tmp_615;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,0)) *
         tmp_614;
   }
   std::complex<double> tmp_616;
   std::complex<double> tmp_617;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_617 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_616 += tmp_617;
   result += (std::complex<double>(0.,0.5163977794943222)*g1*MDBS*ZA(gI2,2)) *
      tmp_616;
   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_619 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_618 += tmp_619;
   result += (std::complex<double>(0.,-0.5163977794943222)*g1*Conj(MDBS)*ZA(gI2
      ,2)) * tmp_618;
   if (gO2 < 3) {
      result += std::complex<double>(0.,-0.12909944487358055)*g1*MDBS*Conj(
         ZU(gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0.,0.12909944487358055)*g1*Conj(MDBS)*
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

   std::complex<double> tmp_620;
   std::complex<double> tmp_621;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_621 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_620 += tmp_621;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_620;
   std::complex<double> tmp_622;
   std::complex<double> tmp_623;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_624;
      std::complex<double> tmp_625;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_625 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_624 += tmp_625;
      tmp_623 += (Conj(ZU(gI1,j2))) * tmp_624;
   }
   tmp_622 += tmp_623;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,0)) * tmp_622;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_626 += tmp_627;
      result += (0.7071067811865475*Mu*ZH(gI2,0)) * tmp_626;
   }
   std::complex<double> tmp_628;
   std::complex<double> tmp_629;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_629 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_628 += tmp_629;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_628;
   std::complex<double> tmp_630;
   std::complex<double> tmp_631;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_632;
      std::complex<double> tmp_633;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_634;
         std::complex<double> tmp_635;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_635 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_634 += tmp_635;
         tmp_633 += (KroneckerDelta(gO2,3 + j2)) * tmp_634;
      }
      tmp_632 += tmp_633;
      tmp_631 += (Conj(ZU(gI1,3 + j3))) * tmp_632;
   }
   tmp_630 += tmp_631;
   result += (-(vu*ZH(gI2,1))) * tmp_630;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_636;
      std::complex<double> tmp_637;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_638;
         std::complex<double> tmp_639;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_639 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_638 += tmp_639;
         tmp_637 += (Conj(ZU(gI1,j2))) * tmp_638;
      }
      tmp_636 += tmp_637;
      result += (-(vu*ZH(gI2,1))) * tmp_636;
   }
   std::complex<double> tmp_640;
   std::complex<double> tmp_641;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_641 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_640 += tmp_641;
   result += (0.5163977794943222*g1*MDBS*ZH(gI2,2)) * tmp_640;
   std::complex<double> tmp_642;
   std::complex<double> tmp_643;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_643 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_642 += tmp_643;
   result += (0.5163977794943222*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_642;
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

   std::complex<double> tmp_644;
   std::complex<double> tmp_645;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_645 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_644 += tmp_645;
   result += (1.4142135623730951*g3) * tmp_644;

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

std::complex<double> CLASSNAME::CpconjUSuSRdpSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_646;
      std::complex<double> tmp_647;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_647 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_646 += tmp_647;
      result += (MuD) * tmp_646;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_648;
      std::complex<double> tmp_649;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_649 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_648 += tmp_649;
      result += (0.7071067811865475*LamSD*vS) * tmp_648;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_650;
      std::complex<double> tmp_651;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_651 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_650 += tmp_651;
      result += (-0.5*LamTD*vT) * tmp_650;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjSRumSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_652;
   std::complex<double> tmp_653;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_654;
      std::complex<double> tmp_655;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_655 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_654 += tmp_655;
      tmp_653 += (Conj(ZD(gI2,j2))) * tmp_654;
   }
   tmp_652 += tmp_653;
   result += (0.5*(-1.4142135623730951*vS*Conj(LamSU) - vT*Conj(LamTU) - 2*Conj
      (MuU))) * tmp_652;

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

   std::complex<double> tmp_656;
   std::complex<double> tmp_657;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_657 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_656 += tmp_657;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_656;
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

   std::complex<double> tmp_658;
   std::complex<double> tmp_659;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_659 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_658 += tmp_659;
   result += (-0.5163977794943222*g1*Sin(ThetaW())) * tmp_658;
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

   std::complex<double> tmp_660;
   std::complex<double> tmp_661;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_661 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_660 += tmp_661;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_660;
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

std::complex<double> CLASSNAME::CpUSeconjUSeconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_662;
   std::complex<double> tmp_663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_663 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_662 += tmp_663;
   result += (-0.3*Sqr(g1)) * tmp_662;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_664;
   std::complex<double> tmp_665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_665 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_664 += tmp_665;
   result += (0.3*Sqr(g1)) * tmp_664;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2);
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

double CLASSNAME::CpconjUSebarCha1FvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSebarCha1FvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_670;
   std::complex<double> tmp_671;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_671 += KroneckerDelta(gO1,3 + j1)*Ye(j1,gI2);
   }
   tmp_670 += tmp_671;
   result += (Conj(UM1(gI1,1))) * tmp_670;

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_672;
   std::complex<double> tmp_674;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_675;
      std::complex<double> tmp_676;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_676 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_675 += tmp_676;
      tmp_674 += (Conj(ZV(gI2,j2))) * tmp_675;
   }
   tmp_672 += tmp_674;
   std::complex<double> tmp_673;
   std::complex<double> tmp_677;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_678;
      std::complex<double> tmp_679;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_679 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_678 += tmp_679;
      tmp_677 += (ZV(gI1,j4)) * tmp_678;
   }
   tmp_673 += tmp_677;
   result += (-1) * tmp_672 * tmp_673;
   std::complex<double> tmp_680;
   std::complex<double> tmp_681;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_681 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_680 += tmp_681;
   result += (0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_680;
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
      std::complex<double> tmp_682;
      std::complex<double> tmp_683;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_684;
         std::complex<double> tmp_685;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_685 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_684 += tmp_685;
         tmp_683 += (Conj(ZV(gI1,j2))) * tmp_684;
      }
      tmp_682 += tmp_683;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_682;
   }
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
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_686;
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
      std::complex<double> tmp_690;
      std::complex<double> tmp_691;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_691 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_690 += tmp_691;
      result += (-ZN2(gI2,2)) * tmp_690;
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

   std::complex<double> tmp_692;
   std::complex<double> tmp_693;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_693 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_692 += tmp_693;
   result += (0.3*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_692;
   std::complex<double> tmp_694;
   std::complex<double> tmp_695;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_696;
      std::complex<double> tmp_697;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_698;
         std::complex<double> tmp_699;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_699 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_698 += tmp_699;
         tmp_697 += (KroneckerDelta(gO2,3 + j2)) * tmp_698;
      }
      tmp_696 += tmp_697;
      tmp_695 += (KroneckerDelta(gO1,3 + j3)) * tmp_696;
   }
   tmp_694 += tmp_695;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_694;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_700;
      std::complex<double> tmp_701;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_701 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_700 += tmp_701;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_700;
   }
   std::complex<double> tmp_702;
   std::complex<double> tmp_703;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_703 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_702 += tmp_703;
   result += (-0.3*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_702;
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

   std::complex<double> tmp_704;
   std::complex<double> tmp_705;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_705 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_704 += tmp_705;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_704;
   std::complex<double> tmp_706;
   std::complex<double> tmp_707;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_708;
      std::complex<double> tmp_709;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_710;
         std::complex<double> tmp_711;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_711 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_710 += tmp_711;
         tmp_709 += (KroneckerDelta(gO2,3 + j2)) * tmp_710;
      }
      tmp_708 += tmp_709;
      tmp_707 += (KroneckerDelta(gO1,3 + j3)) * tmp_708;
   }
   tmp_706 += tmp_707;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_706;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_712;
   std::complex<double> tmp_713;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_713 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_712 += tmp_713;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_712;
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

   std::complex<double> tmp_714;
   std::complex<double> tmp_715;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_715 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_714 += tmp_715;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_714;
   std::complex<double> tmp_716;
   std::complex<double> tmp_717;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_718;
      std::complex<double> tmp_719;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_720;
         std::complex<double> tmp_721;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_721 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_720 += tmp_721;
         tmp_719 += (KroneckerDelta(gO2,3 + j2)) * tmp_720;
      }
      tmp_718 += tmp_719;
      tmp_717 += (KroneckerDelta(gO1,3 + j3)) * tmp_718;
   }
   tmp_716 += tmp_717;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_716;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_722;
      std::complex<double> tmp_723;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_723 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_722 += tmp_723;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_722;
   }
   std::complex<double> tmp_724;
   std::complex<double> tmp_725;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_725 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_724 += tmp_725;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_724;
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

   std::complex<double> tmp_726;
   std::complex<double> tmp_727;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_727 += KroneckerDelta(gO2,3 + j1)*ZER(gI2,j1);
   }
   tmp_726 += tmp_727;
   result += (-1.0954451150103321*g1*ZN1(gI1,0)) * tmp_726;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSebarChiFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_728;
   std::complex<double> tmp_729;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_730;
      std::complex<double> tmp_731;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_731 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_730 += tmp_731;
      tmp_729 += (Conj(ZEL(gI2,j2))) * tmp_730;
   }
   tmp_728 += tmp_729;
   result += (-Conj(ZN2(gI1,2))) * tmp_728;

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_732;
   std::complex<double> tmp_734;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_734 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_732 += tmp_734;
   std::complex<double> tmp_733;
   std::complex<double> tmp_735;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_735 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_733 += tmp_735;
   result += (-0.05*Sqr(g1)) * tmp_732 * tmp_733;
   std::complex<double> tmp_736;
   std::complex<double> tmp_738;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_738 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_736 += tmp_738;
   std::complex<double> tmp_737;
   std::complex<double> tmp_739;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_739 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_737 += tmp_739;
   result += (-0.1*Sqr(g1)) * tmp_736 * tmp_737;
   std::complex<double> tmp_740;
   std::complex<double> tmp_742;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_742 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_740 += tmp_742;
   std::complex<double> tmp_741;
   std::complex<double> tmp_743;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_743 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_741 += tmp_743;
   result += (-0.05*Sqr(g1)) * tmp_740 * tmp_741;
   std::complex<double> tmp_744;
   std::complex<double> tmp_746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_746 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_744 += tmp_746;
   std::complex<double> tmp_745;
   std::complex<double> tmp_747;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_747 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_745 += tmp_747;
   result += (-0.1*Sqr(g1)) * tmp_744 * tmp_745;
   if (gO1 < 3) {
      std::complex<double> tmp_748;
      std::complex<double> tmp_750;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_750 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_748 += tmp_750;
      std::complex<double> tmp_749;
      std::complex<double> tmp_751;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_752;
         std::complex<double> tmp_753;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_753 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_752 += tmp_753;
         tmp_751 += (ZD(gI1,j4)) * tmp_752;
      }
      tmp_749 += tmp_751;
      result += (-1) * tmp_748 * tmp_749;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_754;
      std::complex<double> tmp_755;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_755 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_754 += tmp_755;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_754;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_756;
      std::complex<double> tmp_757;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_757 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_756 += tmp_757;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_756;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_758;
      std::complex<double> tmp_759;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_759 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_758 += tmp_759;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_758;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_760;
      std::complex<double> tmp_761;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_761 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_760 += tmp_761;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_760;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_762;
      std::complex<double> tmp_763;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_763 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_762 += tmp_763;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_762;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_764;
      std::complex<double> tmp_765;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_765 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_764 += tmp_765;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_764;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_766;
      std::complex<double> tmp_768;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_769;
         std::complex<double> tmp_770;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_770 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_769 += tmp_770;
         tmp_768 += (Conj(ZD(gI2,j2))) * tmp_769;
      }
      tmp_766 += tmp_768;
      std::complex<double> tmp_767;
      std::complex<double> tmp_771;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_771 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_767 += tmp_771;
      result += (-1) * tmp_766 * tmp_767;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_772;
   std::complex<double> tmp_774;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_774 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_772 += tmp_774;
   std::complex<double> tmp_773;
   std::complex<double> tmp_775;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_775 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_773 += tmp_775;
   result += (-0.3*Sqr(g1)) * tmp_772 * tmp_773;
   std::complex<double> tmp_776;
   std::complex<double> tmp_778;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_778 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_776 += tmp_778;
   std::complex<double> tmp_777;
   std::complex<double> tmp_779;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_779 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_777 += tmp_779;
   result += (0.15*Sqr(g1)) * tmp_776 * tmp_777;
   std::complex<double> tmp_780;
   std::complex<double> tmp_782;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_782 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_780 += tmp_782;
   std::complex<double> tmp_781;
   std::complex<double> tmp_783;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_783 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_781 += tmp_783;
   result += (-0.3*Sqr(g1)) * tmp_780 * tmp_781;
   std::complex<double> tmp_784;
   std::complex<double> tmp_786;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_786 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_784 += tmp_786;
   std::complex<double> tmp_785;
   std::complex<double> tmp_787;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_787 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_785 += tmp_787;
   result += (0.15*Sqr(g1)) * tmp_784 * tmp_785;
   std::complex<double> tmp_788;
   std::complex<double> tmp_790;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_790 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_788 += tmp_790;
   std::complex<double> tmp_789;
   std::complex<double> tmp_791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_791 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_789 += tmp_791;
   result += (-0.3*Sqr(g1)) * tmp_788 * tmp_789;
   std::complex<double> tmp_792;
   std::complex<double> tmp_794;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_794 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_792 += tmp_794;
   std::complex<double> tmp_793;
   std::complex<double> tmp_795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_795 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_793 += tmp_795;
   result += (-0.3*Sqr(g1)) * tmp_792 * tmp_793;
   std::complex<double> tmp_796;
   std::complex<double> tmp_798;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_799;
      std::complex<double> tmp_800;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_800 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_799 += tmp_800;
      tmp_798 += (Conj(ZE(gI2,j2))) * tmp_799;
   }
   tmp_796 += tmp_798;
   std::complex<double> tmp_797;
   std::complex<double> tmp_801;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_802;
      std::complex<double> tmp_803;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_803 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_802 += tmp_803;
      tmp_801 += (ZE(gI1,j4)) * tmp_802;
   }
   tmp_797 += tmp_801;
   result += (-1) * tmp_796 * tmp_797;
   if (gO1 < 3) {
      std::complex<double> tmp_804;
      std::complex<double> tmp_806;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_806 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_804 += tmp_806;
      std::complex<double> tmp_805;
      std::complex<double> tmp_807;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_808;
         std::complex<double> tmp_809;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_809 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_808 += tmp_809;
         tmp_807 += (ZE(gI1,j4)) * tmp_808;
      }
      tmp_805 += tmp_807;
      result += (-1) * tmp_804 * tmp_805;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_810;
      std::complex<double> tmp_811;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_811 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_810 += tmp_811;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_810;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_812;
      std::complex<double> tmp_813;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_813 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_812 += tmp_813;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_812;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_814;
      std::complex<double> tmp_815;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_815 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_814 += tmp_815;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_814;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_816;
      std::complex<double> tmp_817;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_817 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_816 += tmp_817;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_816;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_818;
      std::complex<double> tmp_819;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_819 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_818 += tmp_819;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_818;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_820;
      std::complex<double> tmp_821;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_821 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_820 += tmp_821;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_820;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_823;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_823 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_822 += tmp_823;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_822;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_824;
      std::complex<double> tmp_825;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_825 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_824 += tmp_825;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_824;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_826;
      std::complex<double> tmp_828;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_829;
         std::complex<double> tmp_830;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_830 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_829 += tmp_830;
         tmp_828 += (Conj(ZE(gI2,j2))) * tmp_829;
      }
      tmp_826 += tmp_828;
      std::complex<double> tmp_827;
      std::complex<double> tmp_831;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_831 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_827 += tmp_831;
      result += (-1) * tmp_826 * tmp_827;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_832;
      std::complex<double> tmp_834;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_834 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_832 += tmp_834;
      std::complex<double> tmp_833;
      std::complex<double> tmp_835;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_835 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_833 += tmp_835;
      result += (-1) * tmp_832 * tmp_833;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_836;
      std::complex<double> tmp_837;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_837 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_836 += tmp_837;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_836;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_838;
      std::complex<double> tmp_839;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_839 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_838 += tmp_839;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_838;
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

   std::complex<double> tmp_840;
   std::complex<double> tmp_842;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_842 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_840 += tmp_842;
   std::complex<double> tmp_841;
   std::complex<double> tmp_843;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_843 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_841 += tmp_843;
   result += (-0.05*Sqr(g1)) * tmp_840 * tmp_841;
   std::complex<double> tmp_844;
   std::complex<double> tmp_846;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_846 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_844 += tmp_846;
   std::complex<double> tmp_845;
   std::complex<double> tmp_847;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_847 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_845 += tmp_847;
   result += (0.2*Sqr(g1)) * tmp_844 * tmp_845;
   std::complex<double> tmp_848;
   std::complex<double> tmp_850;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_850 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_848 += tmp_850;
   std::complex<double> tmp_849;
   std::complex<double> tmp_851;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_851 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_849 += tmp_851;
   result += (-0.05*Sqr(g1)) * tmp_848 * tmp_849;
   std::complex<double> tmp_852;
   std::complex<double> tmp_854;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_854 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_852 += tmp_854;
   std::complex<double> tmp_853;
   std::complex<double> tmp_855;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_855 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_853 += tmp_855;
   result += (0.2*Sqr(g1)) * tmp_852 * tmp_853;
   if (gO1 < 3) {
      std::complex<double> tmp_856;
      std::complex<double> tmp_857;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_857 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_856 += tmp_857;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_856;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_858;
      std::complex<double> tmp_859;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_859 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_858 += tmp_859;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_858;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_860;
      std::complex<double> tmp_861;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_861 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_860 += tmp_861;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_860;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_862;
      std::complex<double> tmp_863;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_863 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_862 += tmp_863;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_862;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_864;
      std::complex<double> tmp_865;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_865 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_864 += tmp_865;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_864;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_866;
      std::complex<double> tmp_867;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_867 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_866 += tmp_867;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_866;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_868;
      std::complex<double> tmp_869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_869 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_868 += tmp_869;
      result += (MuD*ZHR(gI2,0)) * tmp_868;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_870;
      std::complex<double> tmp_871;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_871 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_870 += tmp_871;
      result += (0.7071067811865475*LamSD*vS*ZHR(gI2,0)) * tmp_870;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_872;
      std::complex<double> tmp_873;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_873 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_872 += tmp_873;
      result += (0.5*LamTD*vT*ZHR(gI2,0)) * tmp_872;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_874;
   std::complex<double> tmp_875;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_876;
      std::complex<double> tmp_877;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_877 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_876 += tmp_877;
      tmp_875 += (Conj(ZE(gI1,j2))) * tmp_876;
   }
   tmp_874 += tmp_875;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1))
      * tmp_874;
   if (gO2 < 3) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_879 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_878 += tmp_879;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_878;
   }
   std::complex<double> tmp_880;
   std::complex<double> tmp_881;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_881 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_880 += tmp_881;
   result += (std::complex<double>(0.,-0.7745966692414834)*g1*MDBS*ZA(gI2,2)) *
      tmp_880;
   std::complex<double> tmp_882;
   std::complex<double> tmp_883;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_883 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_882 += tmp_883;
   result += (std::complex<double>(0.,0.7745966692414834)*g1*Conj(MDBS)*ZA(gI2,
      2)) * tmp_882;
   if (gO2 < 3) {
      result += std::complex<double>(0.,0.3872983346207417)*g1*MDBS*Conj(ZE(
         gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0.,-0.3872983346207417)*g1*Conj(MDBS)*
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

   std::complex<double> tmp_884;
   std::complex<double> tmp_885;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_885 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_884 += tmp_885;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_884;
   std::complex<double> tmp_886;
   std::complex<double> tmp_887;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_888;
      std::complex<double> tmp_889;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_890;
         std::complex<double> tmp_891;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_891 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_890 += tmp_891;
         tmp_889 += (KroneckerDelta(gO2,3 + j2)) * tmp_890;
      }
      tmp_888 += tmp_889;
      tmp_887 += (Conj(ZE(gI1,3 + j3))) * tmp_888;
   }
   tmp_886 += tmp_887;
   result += (-(vd*ZH(gI2,0))) * tmp_886;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_892;
      std::complex<double> tmp_893;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_894;
         std::complex<double> tmp_895;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_895 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_894 += tmp_895;
         tmp_893 += (Conj(ZE(gI1,j2))) * tmp_894;
      }
      tmp_892 += tmp_893;
      result += (-(vd*ZH(gI2,0))) * tmp_892;
   }
   std::complex<double> tmp_896;
   std::complex<double> tmp_897;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_897 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_896 += tmp_897;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_896;
   std::complex<double> tmp_898;
   std::complex<double> tmp_899;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_900;
      std::complex<double> tmp_901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_901 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_900 += tmp_901;
      tmp_899 += (Conj(ZE(gI1,j2))) * tmp_900;
   }
   tmp_898 += tmp_899;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_898;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_902;
      std::complex<double> tmp_903;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_903 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_902 += tmp_903;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_902;
   }
   std::complex<double> tmp_904;
   std::complex<double> tmp_905;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_905 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_904 += tmp_905;
   result += (-0.7745966692414834*g1*MDBS*ZH(gI2,2)) * tmp_904;
   std::complex<double> tmp_906;
   std::complex<double> tmp_907;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_907 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_906 += tmp_907;
   result += (-0.7745966692414834*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_906;
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

std::complex<double> CLASSNAME::CpconjUSeconjSRdpSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_908;
   std::complex<double> tmp_909;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_910;
      std::complex<double> tmp_911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_911 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_910 += tmp_911;
      tmp_909 += (Conj(ZV(gI2,j2))) * tmp_910;
   }
   tmp_908 += tmp_909;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(
      MuD))) * tmp_908;

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

   std::complex<double> tmp_912;
   std::complex<double> tmp_913;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_913 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_912 += tmp_913;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_912;
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

   std::complex<double> tmp_914;
   std::complex<double> tmp_915;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_915 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_914 += tmp_915;
   result += (0.7745966692414834*g1*Sin(ThetaW())) * tmp_914;
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

std::complex<double> CLASSNAME::CpUhhconjSRdpSRdp(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(-7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2) - 20*vS*
      AbsSqr(LamSD)*KroneckerDelta(2,gO2) - 14.142135623730951*MuD*Conj(LamSD)*
      KroneckerDelta(2,gO2) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*
      KroneckerDelta(2,gO2) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*
      KroneckerDelta(2,gO2) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2
      ) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2) - 10*g2*MDWBT*
      KroneckerDelta(3,gO2) - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2) +
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2) + 10*MuD*Conj(
      LamTD)*KroneckerDelta(3,gO2) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*
      KroneckerDelta(3,gO2) - 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2) + 10*LamTD*
      Conj(MuD)*KroneckerDelta(3,gO2) + vd*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD
      ) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(1,gO2)*(-3*vu*Sqr(g1) + 5*vu*Sqr
      (g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSRumSRum(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2) - 20*vS*
      AbsSqr(LamSU)*KroneckerDelta(2,gO2) - 14.142135623730951*MuU*Conj(LamSU)*
      KroneckerDelta(2,gO2) - 7.0710678118654755*LamTU*vT*Conj(LamSU)*
      KroneckerDelta(2,gO2) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*
      KroneckerDelta(2,gO2) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2
      ) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2) + 10*g2*MDWBT*
      KroneckerDelta(3,gO2) - 10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2) -
      7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2) - 10*MuU*Conj(
      LamTU)*KroneckerDelta(3,gO2) - 7.0710678118654755*LamSU*vS*Conj(LamTU)*
      KroneckerDelta(3,gO2) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2) - 10*LamTU*
      Conj(MuU)*KroneckerDelta(3,gO2) + vu*KroneckerDelta(1,gO2)*(-20*AbsSqr(LamTU
      ) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(0,gO2)*(-3*vd*Sqr(g1) + 5*vd*Sqr
      (g2)));

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

std::complex<double> CLASSNAME::CpUhhUhhconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,gO2
      )*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1) -
      2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)
      *(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1)
      + 2*LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)) + Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2)))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)
      *(-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

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

   std::complex<double> tmp_916;
   std::complex<double> tmp_917;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_918;
      std::complex<double> tmp_919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_919 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_918 += tmp_919;
      tmp_917 += (ZDL(gI1,j2)) * tmp_918;
   }
   tmp_916 += tmp_917;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_916;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_920;
   std::complex<double> tmp_921;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_922;
      std::complex<double> tmp_923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_923 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_922 += tmp_923;
      tmp_921 += (Conj(ZDL(gI2,j2))) * tmp_922;
   }
   tmp_920 += tmp_921;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_920;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_924;
   std::complex<double> tmp_925;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_926;
      std::complex<double> tmp_927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_927 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_926 += tmp_927;
      tmp_925 += (ZEL(gI1,j2)) * tmp_926;
   }
   tmp_924 += tmp_925;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_924;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_928;
   std::complex<double> tmp_929;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_930;
      std::complex<double> tmp_931;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_931 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_930 += tmp_931;
      tmp_929 += (Conj(ZEL(gI2,j2))) * tmp_930;
   }
   tmp_928 += tmp_929;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_928;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_932;
   std::complex<double> tmp_933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_934;
      std::complex<double> tmp_935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_935 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_934 += tmp_935;
      tmp_933 += (ZUL(gI1,j2)) * tmp_934;
   }
   tmp_932 += tmp_933;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_932;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_936;
   std::complex<double> tmp_937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_938;
      std::complex<double> tmp_939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_939 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_938 += tmp_939;
      tmp_937 += (Conj(ZUL(gI2,j2))) * tmp_938;
   }
   tmp_936 += tmp_937;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_936;

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

   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   std::complex<double> tmp_942;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_942 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_941 += tmp_942;
   tmp_940 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_941;
   std::complex<double> tmp_943;
   std::complex<double> tmp_944;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_944 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_943 += tmp_944;
   tmp_940 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_943;
   std::complex<double> tmp_945;
   std::complex<double> tmp_946;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_946 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_945 += tmp_946;
   tmp_940 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*Sqr(g1)) * tmp_945;
   std::complex<double> tmp_947;
   std::complex<double> tmp_948;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_949;
      std::complex<double> tmp_950;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_951;
         std::complex<double> tmp_952;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_952 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_951 += tmp_952;
         tmp_950 += (ZD(gI1,3 + j2)) * tmp_951;
      }
      tmp_949 += tmp_950;
      tmp_948 += (Conj(ZD(gI2,3 + j3))) * tmp_949;
   }
   tmp_947 += tmp_948;
   tmp_940 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_947;
   std::complex<double> tmp_953;
   std::complex<double> tmp_954;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_955;
      std::complex<double> tmp_956;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_957;
         std::complex<double> tmp_958;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_958 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_957 += tmp_958;
         tmp_956 += (Conj(ZD(gI2,j2))) * tmp_957;
      }
      tmp_955 += tmp_956;
      tmp_954 += (ZD(gI1,j3)) * tmp_955;
   }
   tmp_953 += tmp_954;
   tmp_940 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_953;
   std::complex<double> tmp_959;
   std::complex<double> tmp_960;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_960 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_959 += tmp_960;
   tmp_940 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_959;
   std::complex<double> tmp_961;
   std::complex<double> tmp_962;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_962 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_961 += tmp_962;
   tmp_940 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_961;
   std::complex<double> tmp_963;
   std::complex<double> tmp_964;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_964 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_963 += tmp_964;
   tmp_940 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_963;
   result += (std::complex<double>(0,-1)) * tmp_940;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_965;
   std::complex<double> tmp_966;
   std::complex<double> tmp_967;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_967 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_966 += tmp_967;
   tmp_965 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_966;
   std::complex<double> tmp_968;
   std::complex<double> tmp_969;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_969 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_968 += tmp_969;
   tmp_965 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_968;
   std::complex<double> tmp_970;
   std::complex<double> tmp_971;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_971 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_970 += tmp_971;
   tmp_965 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*Sqr(g1)) * tmp_970;
   std::complex<double> tmp_972;
   std::complex<double> tmp_973;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_974;
      std::complex<double> tmp_975;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_976;
         std::complex<double> tmp_977;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_977 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_976 += tmp_977;
         tmp_975 += (ZE(gI1,3 + j2)) * tmp_976;
      }
      tmp_974 += tmp_975;
      tmp_973 += (Conj(ZE(gI2,3 + j3))) * tmp_974;
   }
   tmp_972 += tmp_973;
   tmp_965 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_972;
   std::complex<double> tmp_978;
   std::complex<double> tmp_979;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_980;
      std::complex<double> tmp_981;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_982;
         std::complex<double> tmp_983;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_983 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_982 += tmp_983;
         tmp_981 += (Conj(ZE(gI2,j2))) * tmp_982;
      }
      tmp_980 += tmp_981;
      tmp_979 += (ZE(gI1,j3)) * tmp_980;
   }
   tmp_978 += tmp_979;
   tmp_965 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_978;
   std::complex<double> tmp_984;
   std::complex<double> tmp_985;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_985 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_984 += tmp_985;
   tmp_965 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_984;
   std::complex<double> tmp_986;
   std::complex<double> tmp_987;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_987 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_986 += tmp_987;
   tmp_965 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_986;
   std::complex<double> tmp_988;
   std::complex<double> tmp_989;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_989 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_988 += tmp_989;
   tmp_965 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_988;
   result += (std::complex<double>(0,-1)) * tmp_965;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_990;
   std::complex<double> tmp_991;
   std::complex<double> tmp_992;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_992 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_991 += tmp_992;
   tmp_990 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_991;
   std::complex<double> tmp_993;
   std::complex<double> tmp_994;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_994 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_993 += tmp_994;
   tmp_990 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_993;
   std::complex<double> tmp_995;
   std::complex<double> tmp_996;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_996 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_995 += tmp_996;
   tmp_990 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_995;
   std::complex<double> tmp_997;
   std::complex<double> tmp_998;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_998 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_997 += tmp_998;
   tmp_990 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_997;
   std::complex<double> tmp_999;
   std::complex<double> tmp_1000;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1000 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_999 += tmp_1000;
   tmp_990 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_999;
   std::complex<double> tmp_1001;
   std::complex<double> tmp_1002;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1002 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1001 += tmp_1002;
   tmp_990 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)*Sqr(g1)) * tmp_1001;
   std::complex<double> tmp_1003;
   std::complex<double> tmp_1004;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1005;
      std::complex<double> tmp_1006;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1007;
         std::complex<double> tmp_1008;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1008 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1007 += tmp_1008;
         tmp_1006 += (ZU(gI1,3 + j2)) * tmp_1007;
      }
      tmp_1005 += tmp_1006;
      tmp_1004 += (Conj(ZU(gI2,3 + j3))) * tmp_1005;
   }
   tmp_1003 += tmp_1004;
   tmp_990 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)) * tmp_1003;
   std::complex<double> tmp_1009;
   std::complex<double> tmp_1010;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1011;
      std::complex<double> tmp_1012;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1013;
         std::complex<double> tmp_1014;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1014 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1013 += tmp_1014;
         tmp_1012 += (Conj(ZU(gI2,j2))) * tmp_1013;
      }
      tmp_1011 += tmp_1012;
      tmp_1010 += (ZU(gI1,j3)) * tmp_1011;
   }
   tmp_1009 += tmp_1010;
   tmp_990 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)) * tmp_1009;
   result += (std::complex<double>(0,-1)) * tmp_990;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1015;
   std::complex<double> tmp_1016;
   std::complex<double> tmp_1017;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1017 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1016 += tmp_1017;
   tmp_1015 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1016;
   std::complex<double> tmp_1018;
   std::complex<double> tmp_1019;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1019 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1018 += tmp_1019;
   tmp_1015 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1018;
   std::complex<double> tmp_1020;
   std::complex<double> tmp_1021;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1021 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1020 += tmp_1021;
   tmp_1015 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1020;
   std::complex<double> tmp_1022;
   std::complex<double> tmp_1023;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1024;
      std::complex<double> tmp_1025;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1026;
         std::complex<double> tmp_1027;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1027 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1026 += tmp_1027;
         tmp_1025 += (ZD(gI1,3 + j2)) * tmp_1026;
      }
      tmp_1024 += tmp_1025;
      tmp_1023 += (Conj(ZD(gI2,3 + j3))) * tmp_1024;
   }
   tmp_1022 += tmp_1023;
   tmp_1015 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1022
      ;
   std::complex<double> tmp_1028;
   std::complex<double> tmp_1029;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1030;
      std::complex<double> tmp_1031;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1032;
         std::complex<double> tmp_1033;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1033 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1032 += tmp_1033;
         tmp_1031 += (Conj(ZD(gI2,j2))) * tmp_1032;
      }
      tmp_1030 += tmp_1031;
      tmp_1029 += (ZD(gI1,j3)) * tmp_1030;
   }
   tmp_1028 += tmp_1029;
   tmp_1015 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1028
      ;
   std::complex<double> tmp_1034;
   std::complex<double> tmp_1035;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1035 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1034 += tmp_1035;
   tmp_1015 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1034;
   std::complex<double> tmp_1036;
   std::complex<double> tmp_1037;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1037 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1036 += tmp_1037;
   tmp_1015 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1036;
   std::complex<double> tmp_1038;
   std::complex<double> tmp_1039;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1039 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1038 += tmp_1039;
   tmp_1015 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1038;
   std::complex<double> tmp_1040;
   std::complex<double> tmp_1041;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1042;
      std::complex<double> tmp_1043;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1043 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1042 += tmp_1043;
      tmp_1041 += (Conj(ZD(gI2,j2))) * tmp_1042;
   }
   tmp_1040 += tmp_1041;
   tmp_1015 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1040;
   std::complex<double> tmp_1044;
   std::complex<double> tmp_1045;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1046;
      std::complex<double> tmp_1047;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1047 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1046 += tmp_1047;
      tmp_1045 += (ZD(gI1,j2)) * tmp_1046;
   }
   tmp_1044 += tmp_1045;
   tmp_1015 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,
      gO2)*Mu) * tmp_1044;
   std::complex<double> tmp_1048;
   std::complex<double> tmp_1049;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1049 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1048 += tmp_1049;
   tmp_1015 += (std::complex<double>(0.,-0.12909944487358055)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1048;
   std::complex<double> tmp_1050;
   std::complex<double> tmp_1051;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1051 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1050 += tmp_1051;
   tmp_1015 += (std::complex<double>(0.,-0.12909944487358055)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1050;
   std::complex<double> tmp_1052;
   std::complex<double> tmp_1053;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1053 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1052 += tmp_1053;
   tmp_1015 += (std::complex<double>(0.,-0.2581988897471611)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1052;
   std::complex<double> tmp_1054;
   std::complex<double> tmp_1055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1055 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1054 += tmp_1055;
   tmp_1015 += (std::complex<double>(0.,-0.2581988897471611)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1054;
   std::complex<double> tmp_1056;
   std::complex<double> tmp_1057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1057 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1056 += tmp_1057;
   tmp_1015 += (std::complex<double>(0,0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1056;
   std::complex<double> tmp_1058;
   std::complex<double> tmp_1059;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1059 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1058 += tmp_1059;
   tmp_1015 += (std::complex<double>(0,0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,gO2
      )) * tmp_1058;
   result += (std::complex<double>(0,-1)) * tmp_1015;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1060;
   std::complex<double> tmp_1061;
   std::complex<double> tmp_1062;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1062 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1061 += tmp_1062;
   tmp_1060 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1061;
   std::complex<double> tmp_1063;
   std::complex<double> tmp_1064;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1064 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1063 += tmp_1064;
   tmp_1060 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1063;
   std::complex<double> tmp_1065;
   std::complex<double> tmp_1066;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1066 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1065 += tmp_1066;
   tmp_1060 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1065;
   std::complex<double> tmp_1067;
   std::complex<double> tmp_1068;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1069;
      std::complex<double> tmp_1070;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1071;
         std::complex<double> tmp_1072;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1072 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1071 += tmp_1072;
         tmp_1070 += (ZE(gI1,3 + j2)) * tmp_1071;
      }
      tmp_1069 += tmp_1070;
      tmp_1068 += (Conj(ZE(gI2,3 + j3))) * tmp_1069;
   }
   tmp_1067 += tmp_1068;
   tmp_1060 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1067
      ;
   std::complex<double> tmp_1073;
   std::complex<double> tmp_1074;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1075;
      std::complex<double> tmp_1076;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1077;
         std::complex<double> tmp_1078;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1078 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1077 += tmp_1078;
         tmp_1076 += (Conj(ZE(gI2,j2))) * tmp_1077;
      }
      tmp_1075 += tmp_1076;
      tmp_1074 += (ZE(gI1,j3)) * tmp_1075;
   }
   tmp_1073 += tmp_1074;
   tmp_1060 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1073
      ;
   std::complex<double> tmp_1079;
   std::complex<double> tmp_1080;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1080 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1079 += tmp_1080;
   tmp_1060 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1079;
   std::complex<double> tmp_1081;
   std::complex<double> tmp_1082;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1082 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1081 += tmp_1082;
   tmp_1060 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1081;
   std::complex<double> tmp_1083;
   std::complex<double> tmp_1084;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1084 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1083 += tmp_1084;
   tmp_1060 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1083;
   std::complex<double> tmp_1085;
   std::complex<double> tmp_1086;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1087;
      std::complex<double> tmp_1088;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1088 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1087 += tmp_1088;
      tmp_1086 += (Conj(ZE(gI2,j2))) * tmp_1087;
   }
   tmp_1085 += tmp_1086;
   tmp_1060 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1085;
   std::complex<double> tmp_1089;
   std::complex<double> tmp_1090;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1091;
      std::complex<double> tmp_1092;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1092 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1091 += tmp_1092;
      tmp_1090 += (ZE(gI1,j2)) * tmp_1091;
   }
   tmp_1089 += tmp_1090;
   tmp_1060 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,
      gO2)*Mu) * tmp_1089;
   std::complex<double> tmp_1093;
   std::complex<double> tmp_1094;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1094 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1093 += tmp_1094;
   tmp_1060 += (std::complex<double>(0.,0.3872983346207417)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1093;
   std::complex<double> tmp_1095;
   std::complex<double> tmp_1096;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1096 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1095 += tmp_1096;
   tmp_1060 += (std::complex<double>(0.,0.3872983346207417)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1095;
   std::complex<double> tmp_1097;
   std::complex<double> tmp_1098;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1098 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1097 += tmp_1098;
   tmp_1060 += (std::complex<double>(0.,-0.7745966692414834)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1097;
   std::complex<double> tmp_1099;
   std::complex<double> tmp_1100;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1100 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1099 += tmp_1100;
   tmp_1060 += (std::complex<double>(0.,-0.7745966692414834)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1099;
   std::complex<double> tmp_1101;
   std::complex<double> tmp_1102;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1102 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1101 += tmp_1102;
   tmp_1060 += (std::complex<double>(0,0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1101;
   std::complex<double> tmp_1103;
   std::complex<double> tmp_1104;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1104 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1103 += tmp_1104;
   tmp_1060 += (std::complex<double>(0,0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,gO2
      )) * tmp_1103;
   result += (std::complex<double>(0,-1)) * tmp_1060;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1105;
   std::complex<double> tmp_1106;
   std::complex<double> tmp_1107;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1107 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1106 += tmp_1107;
   tmp_1105 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1106;
   std::complex<double> tmp_1108;
   std::complex<double> tmp_1109;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1109 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1108 += tmp_1109;
   tmp_1105 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1108;
   std::complex<double> tmp_1110;
   std::complex<double> tmp_1111;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1111 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1110 += tmp_1111;
   tmp_1105 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1110;
   std::complex<double> tmp_1112;
   std::complex<double> tmp_1113;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1114;
      std::complex<double> tmp_1115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1115 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1114 += tmp_1115;
      tmp_1113 += (Conj(ZU(gI2,j2))) * tmp_1114;
   }
   tmp_1112 += tmp_1113;
   tmp_1105 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(0,gO2)) * tmp_1112;
   std::complex<double> tmp_1116;
   std::complex<double> tmp_1117;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1118;
      std::complex<double> tmp_1119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1119 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1118 += tmp_1119;
      tmp_1117 += (ZU(gI1,j2)) * tmp_1118;
   }
   tmp_1116 += tmp_1117;
   tmp_1105 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,
      gO2)*Mu) * tmp_1116;
   std::complex<double> tmp_1120;
   std::complex<double> tmp_1121;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1121 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1120 += tmp_1121;
   tmp_1105 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1120;
   std::complex<double> tmp_1122;
   std::complex<double> tmp_1123;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1123 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1122 += tmp_1123;
   tmp_1105 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1122;
   std::complex<double> tmp_1124;
   std::complex<double> tmp_1125;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1125 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1124 += tmp_1125;
   tmp_1105 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1124;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1127;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1128;
      std::complex<double> tmp_1129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1130;
         std::complex<double> tmp_1131;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1131 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1130 += tmp_1131;
         tmp_1129 += (ZU(gI1,3 + j2)) * tmp_1130;
      }
      tmp_1128 += tmp_1129;
      tmp_1127 += (Conj(ZU(gI2,3 + j3))) * tmp_1128;
   }
   tmp_1126 += tmp_1127;
   tmp_1105 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1126
      ;
   std::complex<double> tmp_1132;
   std::complex<double> tmp_1133;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1134;
      std::complex<double> tmp_1135;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1136;
         std::complex<double> tmp_1137;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1137 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1136 += tmp_1137;
         tmp_1135 += (Conj(ZU(gI2,j2))) * tmp_1136;
      }
      tmp_1134 += tmp_1135;
      tmp_1133 += (ZU(gI1,j3)) * tmp_1134;
   }
   tmp_1132 += tmp_1133;
   tmp_1105 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1132
      ;
   std::complex<double> tmp_1138;
   std::complex<double> tmp_1139;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1139 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1138 += tmp_1139;
   tmp_1105 += (std::complex<double>(0.,-0.12909944487358055)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1138;
   std::complex<double> tmp_1140;
   std::complex<double> tmp_1141;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1141 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1140 += tmp_1141;
   tmp_1105 += (std::complex<double>(0.,-0.12909944487358055)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1140;
   std::complex<double> tmp_1142;
   std::complex<double> tmp_1143;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1143 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1142 += tmp_1143;
   tmp_1105 += (std::complex<double>(0.,0.5163977794943222)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1142;
   std::complex<double> tmp_1144;
   std::complex<double> tmp_1145;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1145 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1144 += tmp_1145;
   tmp_1105 += (std::complex<double>(0.,0.5163977794943222)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1144;
   std::complex<double> tmp_1146;
   std::complex<double> tmp_1147;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1147 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1146 += tmp_1147;
   tmp_1105 += (std::complex<double>(0,-0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1146;
   std::complex<double> tmp_1148;
   std::complex<double> tmp_1149;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1149 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1148 += tmp_1149;
   tmp_1105 += (std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,
      gO2)) * tmp_1148;
   result += (std::complex<double>(0,-1)) * tmp_1105;

   return result;
}

std::complex<double> CLASSNAME::CpUhhSRdpHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Conj(Mu)*(-1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*ZP(
      gI2,1) + LamTD*KroneckerDelta(3,gO2)*ZP(gI2,1) + 1.4142135623730951*LamTD*
      KroneckerDelta(1,gO2)*ZP(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSRumHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Mu*(1.4142135623730951*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI2
      ,0) + Conj(LamTU)*(KroneckerDelta(3,gO2)*ZP(gI2,0) - 1.4142135623730951*
      KroneckerDelta(0,gO2)*ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpUAhconjSRdpSRdp(unsigned gO2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.05)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(-2*MuD + LamTD*vT)*Conj(LamSD) - 7.0710678118654755*
      LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) + 14.142135623730951*
      LamSD*Conj(MuD))*KroneckerDelta(2,gO2) + 5*(2*g2*MDWBT - 1.4142135623730951*
      LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) + 1.4142135623730951*LamSD*vS*Conj(
      LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*KroneckerDelta(3,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSRumSRum(unsigned gO2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU) - 7.0710678118654755*LamSU
      *vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSU
      *Conj(MuU))*KroneckerDelta(2,gO2) + 5*(2*g2*MDWBT - 1.4142135623730951*LamTU
      *vS*Conj(LamSU) + 2*MuU*Conj(LamTU) + 1.4142135623730951*LamSU*vS*Conj(LamTU
      ) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(MuU))*KroneckerDelta(3,gO2));

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

std::complex<double> CLASSNAME::CpUAhUAhconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,gO2
      )*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1) -
      2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)
      *(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1)
      + 2*LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)) + Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2)))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)
      *(-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

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

   std::complex<double> tmp_1150;
   std::complex<double> tmp_1151;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1152;
      std::complex<double> tmp_1153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1153 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1152 += tmp_1153;
      tmp_1151 += (ZDL(gI1,j2)) * tmp_1152;
   }
   tmp_1150 += tmp_1151;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1150;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1154;
   std::complex<double> tmp_1155;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1156;
      std::complex<double> tmp_1157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1157 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1156 += tmp_1157;
      tmp_1155 += (Conj(ZDL(gI2,j2))) * tmp_1156;
   }
   tmp_1154 += tmp_1155;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1154;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1158;
   std::complex<double> tmp_1159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1160;
      std::complex<double> tmp_1161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1161 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1160 += tmp_1161;
      tmp_1159 += (ZEL(gI1,j2)) * tmp_1160;
   }
   tmp_1158 += tmp_1159;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1158;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1162;
   std::complex<double> tmp_1163;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1164;
      std::complex<double> tmp_1165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1165 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1164 += tmp_1165;
      tmp_1163 += (Conj(ZEL(gI2,j2))) * tmp_1164;
   }
   tmp_1162 += tmp_1163;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1162;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1166;
   std::complex<double> tmp_1167;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1168;
      std::complex<double> tmp_1169;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1169 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1168 += tmp_1169;
      tmp_1167 += (ZUL(gI1,j2)) * tmp_1168;
   }
   tmp_1166 += tmp_1167;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_1166;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1170;
   std::complex<double> tmp_1171;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1172;
      std::complex<double> tmp_1173;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1173 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1172 += tmp_1173;
      tmp_1171 += (Conj(ZUL(gI2,j2))) * tmp_1172;
   }
   tmp_1170 += tmp_1171;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_1170;

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

   std::complex<double> tmp_1174;
   std::complex<double> tmp_1175;
   std::complex<double> tmp_1176;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1176 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1175 += tmp_1176;
   tmp_1174 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1175;
   std::complex<double> tmp_1177;
   std::complex<double> tmp_1178;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1178 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1177 += tmp_1178;
   tmp_1174 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1177;
   std::complex<double> tmp_1179;
   std::complex<double> tmp_1180;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1180 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1179 += tmp_1180;
   tmp_1174 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1179;
   std::complex<double> tmp_1181;
   std::complex<double> tmp_1182;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1183;
      std::complex<double> tmp_1184;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1185;
         std::complex<double> tmp_1186;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1186 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1185 += tmp_1186;
         tmp_1184 += (ZD(gI1,3 + j2)) * tmp_1185;
      }
      tmp_1183 += tmp_1184;
      tmp_1182 += (Conj(ZD(gI2,3 + j3))) * tmp_1183;
   }
   tmp_1181 += tmp_1182;
   tmp_1174 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1181;
   std::complex<double> tmp_1187;
   std::complex<double> tmp_1188;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1189;
      std::complex<double> tmp_1190;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1191;
         std::complex<double> tmp_1192;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1192 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1191 += tmp_1192;
         tmp_1190 += (Conj(ZD(gI2,j2))) * tmp_1191;
      }
      tmp_1189 += tmp_1190;
      tmp_1188 += (ZD(gI1,j3)) * tmp_1189;
   }
   tmp_1187 += tmp_1188;
   tmp_1174 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1187;
   std::complex<double> tmp_1193;
   std::complex<double> tmp_1194;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1194 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1193 += tmp_1194;
   tmp_1174 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1193;
   std::complex<double> tmp_1195;
   std::complex<double> tmp_1196;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1196 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1195 += tmp_1196;
   tmp_1174 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1195;
   std::complex<double> tmp_1197;
   std::complex<double> tmp_1198;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1198 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1197 += tmp_1198;
   tmp_1174 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1197;
   result += (std::complex<double>(0,-1)) * tmp_1174;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1199;
   std::complex<double> tmp_1200;
   std::complex<double> tmp_1201;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1201 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1200 += tmp_1201;
   tmp_1199 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1200;
   std::complex<double> tmp_1202;
   std::complex<double> tmp_1203;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1203 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1202 += tmp_1203;
   tmp_1199 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1202;
   std::complex<double> tmp_1204;
   std::complex<double> tmp_1205;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1205 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1204 += tmp_1205;
   tmp_1199 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1204;
   std::complex<double> tmp_1206;
   std::complex<double> tmp_1207;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1208;
      std::complex<double> tmp_1209;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1210;
         std::complex<double> tmp_1211;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1211 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1210 += tmp_1211;
         tmp_1209 += (ZE(gI1,3 + j2)) * tmp_1210;
      }
      tmp_1208 += tmp_1209;
      tmp_1207 += (Conj(ZE(gI2,3 + j3))) * tmp_1208;
   }
   tmp_1206 += tmp_1207;
   tmp_1199 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1206;
   std::complex<double> tmp_1212;
   std::complex<double> tmp_1213;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1214;
      std::complex<double> tmp_1215;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1216;
         std::complex<double> tmp_1217;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1217 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1216 += tmp_1217;
         tmp_1215 += (Conj(ZE(gI2,j2))) * tmp_1216;
      }
      tmp_1214 += tmp_1215;
      tmp_1213 += (ZE(gI1,j3)) * tmp_1214;
   }
   tmp_1212 += tmp_1213;
   tmp_1199 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1212;
   std::complex<double> tmp_1218;
   std::complex<double> tmp_1219;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1219 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1218 += tmp_1219;
   tmp_1199 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1218;
   std::complex<double> tmp_1220;
   std::complex<double> tmp_1221;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1221 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1220 += tmp_1221;
   tmp_1199 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1220;
   std::complex<double> tmp_1222;
   std::complex<double> tmp_1223;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1223 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1222 += tmp_1223;
   tmp_1199 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1222;
   result += (std::complex<double>(0,-1)) * tmp_1199;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1224;
   std::complex<double> tmp_1225;
   std::complex<double> tmp_1226;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1226 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1225 += tmp_1226;
   tmp_1224 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1225;
   std::complex<double> tmp_1227;
   std::complex<double> tmp_1228;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1228 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1227 += tmp_1228;
   tmp_1224 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1227;
   std::complex<double> tmp_1229;
   std::complex<double> tmp_1230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1230 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1229 += tmp_1230;
   tmp_1224 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1229;
   std::complex<double> tmp_1231;
   std::complex<double> tmp_1232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1232 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1231 += tmp_1232;
   tmp_1224 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1231;
   std::complex<double> tmp_1233;
   std::complex<double> tmp_1234;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1234 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1233 += tmp_1234;
   tmp_1224 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1233;
   std::complex<double> tmp_1235;
   std::complex<double> tmp_1236;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1236 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1235 += tmp_1236;
   tmp_1224 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1235;
   std::complex<double> tmp_1237;
   std::complex<double> tmp_1238;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1239;
      std::complex<double> tmp_1240;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1241;
         std::complex<double> tmp_1242;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1242 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1241 += tmp_1242;
         tmp_1240 += (ZU(gI1,3 + j2)) * tmp_1241;
      }
      tmp_1239 += tmp_1240;
      tmp_1238 += (Conj(ZU(gI2,3 + j3))) * tmp_1239;
   }
   tmp_1237 += tmp_1238;
   tmp_1224 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1237;
   std::complex<double> tmp_1243;
   std::complex<double> tmp_1244;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1245;
      std::complex<double> tmp_1246;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1247;
         std::complex<double> tmp_1248;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1248 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1247 += tmp_1248;
         tmp_1246 += (Conj(ZU(gI2,j2))) * tmp_1247;
      }
      tmp_1245 += tmp_1246;
      tmp_1244 += (ZU(gI1,j3)) * tmp_1245;
   }
   tmp_1243 += tmp_1244;
   tmp_1224 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1243;
   result += (std::complex<double>(0,-1)) * tmp_1224;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1249;
   std::complex<double> tmp_1250;
   std::complex<double> tmp_1251;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1252;
      std::complex<double> tmp_1253;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1253 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1252 += tmp_1253;
      tmp_1251 += (Conj(ZD(gI2,j2))) * tmp_1252;
   }
   tmp_1250 += tmp_1251;
   tmp_1249 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1250;
   std::complex<double> tmp_1254;
   std::complex<double> tmp_1255;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1256;
      std::complex<double> tmp_1257;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1257 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1256 += tmp_1257;
      tmp_1255 += (ZD(gI1,j2)) * tmp_1256;
   }
   tmp_1254 += tmp_1255;
   tmp_1249 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1254;
   std::complex<double> tmp_1258;
   std::complex<double> tmp_1259;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1259 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1258 += tmp_1259;
   tmp_1249 += (0.12909944487358055*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1258;
   std::complex<double> tmp_1260;
   std::complex<double> tmp_1261;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1261 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1260 += tmp_1261;
   tmp_1249 += (-0.12909944487358055*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1260;
   std::complex<double> tmp_1262;
   std::complex<double> tmp_1263;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1263 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1262 += tmp_1263;
   tmp_1249 += (0.2581988897471611*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1262;
   std::complex<double> tmp_1264;
   std::complex<double> tmp_1265;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1265 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1264 += tmp_1265;
   tmp_1249 += (-0.2581988897471611*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1264;
   std::complex<double> tmp_1266;
   std::complex<double> tmp_1267;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1267 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1266 += tmp_1267;
   tmp_1249 += (-0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1266;
   std::complex<double> tmp_1268;
   std::complex<double> tmp_1269;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1269 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1268 += tmp_1269;
   tmp_1249 += (0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1268;
   result += (std::complex<double>(0,-1)) * tmp_1249;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1270;
   std::complex<double> tmp_1271;
   std::complex<double> tmp_1272;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1273;
      std::complex<double> tmp_1274;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1274 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1273 += tmp_1274;
      tmp_1272 += (Conj(ZE(gI2,j2))) * tmp_1273;
   }
   tmp_1271 += tmp_1272;
   tmp_1270 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1271;
   std::complex<double> tmp_1275;
   std::complex<double> tmp_1276;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1277;
      std::complex<double> tmp_1278;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1278 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1277 += tmp_1278;
      tmp_1276 += (ZE(gI1,j2)) * tmp_1277;
   }
   tmp_1275 += tmp_1276;
   tmp_1270 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1275;
   std::complex<double> tmp_1279;
   std::complex<double> tmp_1280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1280 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1279 += tmp_1280;
   tmp_1270 += (-0.3872983346207417*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1279;
   std::complex<double> tmp_1281;
   std::complex<double> tmp_1282;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1282 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1281 += tmp_1282;
   tmp_1270 += (0.3872983346207417*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1281;
   std::complex<double> tmp_1283;
   std::complex<double> tmp_1284;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1284 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1283 += tmp_1284;
   tmp_1270 += (0.7745966692414834*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1283;
   std::complex<double> tmp_1285;
   std::complex<double> tmp_1286;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1286 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1285 += tmp_1286;
   tmp_1270 += (-0.7745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1285;
   std::complex<double> tmp_1287;
   std::complex<double> tmp_1288;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1288 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1287 += tmp_1288;
   tmp_1270 += (-0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1287;
   std::complex<double> tmp_1289;
   std::complex<double> tmp_1290;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1290 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1289 += tmp_1290;
   tmp_1270 += (0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1289;
   result += (std::complex<double>(0,-1)) * tmp_1270;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1291;
   std::complex<double> tmp_1292;
   std::complex<double> tmp_1293;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1294;
      std::complex<double> tmp_1295;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1295 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1294 += tmp_1295;
      tmp_1293 += (Conj(ZU(gI2,j2))) * tmp_1294;
   }
   tmp_1292 += tmp_1293;
   tmp_1291 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(0,gO2)) * tmp_1292;
   std::complex<double> tmp_1296;
   std::complex<double> tmp_1297;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1298;
      std::complex<double> tmp_1299;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1299 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1298 += tmp_1299;
      tmp_1297 += (ZU(gI1,j2)) * tmp_1298;
   }
   tmp_1296 += tmp_1297;
   tmp_1291 += (-0.7071067811865475*KroneckerDelta(0,gO2)*Mu) * tmp_1296;
   std::complex<double> tmp_1300;
   std::complex<double> tmp_1301;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1301 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1300 += tmp_1301;
   tmp_1291 += (0.12909944487358055*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1300;
   std::complex<double> tmp_1302;
   std::complex<double> tmp_1303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1303 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1302 += tmp_1303;
   tmp_1291 += (-0.12909944487358055*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1302;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1305 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1304 += tmp_1305;
   tmp_1291 += (-0.5163977794943222*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1304;
   std::complex<double> tmp_1306;
   std::complex<double> tmp_1307;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1307 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1306 += tmp_1307;
   tmp_1291 += (0.5163977794943222*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1306;
   std::complex<double> tmp_1308;
   std::complex<double> tmp_1309;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1309 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1308 += tmp_1309;
   tmp_1291 += (0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1308;
   std::complex<double> tmp_1310;
   std::complex<double> tmp_1311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1311 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1310 += tmp_1311;
   tmp_1291 += (-0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1310;
   result += (std::complex<double>(0,-1)) * tmp_1291;

   return result;
}

std::complex<double> CLASSNAME::CpUAhSRdpHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*Conj(Mu)*(1.4142135623730951*LamSD*
      KroneckerDelta(2,gO2)*ZP(gI2,1) - LamTD*KroneckerDelta(3,gO2)*ZP(gI2,1) +
      1.4142135623730951*LamTD*KroneckerDelta(1,gO2)*ZP(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSRumHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*Mu*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(2,gO2)*ZP(gI2,0) + Conj(LamTU)*(KroneckerDelta(3,gO2)*ZP(gI2,
      0) + 1.4142135623730951*KroneckerDelta(0,gO2)*ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpURhconjURhconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) - 5*
      Sqr(g2)) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2
      )));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*
      Sqr(g2)) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2
      )));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjURhVWmSRdp(unsigned gO2) const
{
   double result = 0.0;

   result = -0.7071067811865475*g2*KroneckerDelta(0,gO2);

   return result;
}

double CLASSNAME::CpconjURhconjVWmSRum(unsigned gO2) const
{
   double result = 0.0;

   result = -0.7071067811865475*g2*KroneckerDelta(1,gO2);

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

std::complex<double> CLASSNAME::CpconjURhconjHpmSRum(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = 0.25*(-2.8284271247461903*LamSU*vd*Conj(LamSD)*KroneckerDelta(0,gO2
      )*ZP(gI1,1) - 1.4142135623730951*LamTU*Conj(LamTD)*KroneckerDelta(0,gO2)*(2*
      vu*ZP(gI1,0) + vd*ZP(gI1,1)) - KroneckerDelta(1,gO2)*(1.4142135623730951*vd*
      Sqr(g2)*ZP(gI1,0) + 1.4142135623730951*vu*(-2*AbsSqr(LamSU) + AbsSqr(LamTU)
      + Sqr(g2))*ZP(gI1,1) + 2*((-((2*MuU + 1.4142135623730951*LamSU*vS + LamTU*vT
      )*Conj(LamTU)) + g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gI1,2) - (-2*g2*MDWBT - vT*
      AbsSqr(LamTU) + 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(MuU)
      + vT*Sqr(g2))*ZP(gI1,3))));

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

   std::complex<double> tmp_1312;
   std::complex<double> tmp_1313;
   std::complex<double> tmp_1314;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1314 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1313 += tmp_1314;
   tmp_1312 += (std::complex<double>(0,-0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1313;
   std::complex<double> tmp_1315;
   std::complex<double> tmp_1316;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1316 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1315 += tmp_1316;
   tmp_1312 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1315;
   std::complex<double> tmp_1317;
   std::complex<double> tmp_1318;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1318 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1317 += tmp_1318;
   tmp_1312 += (std::complex<double>(0,-0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1317;
   std::complex<double> tmp_1319;
   std::complex<double> tmp_1320;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1320 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1319 += tmp_1320;
   tmp_1312 += (std::complex<double>(0,0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1319;
   std::complex<double> tmp_1321;
   std::complex<double> tmp_1322;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1322 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1321 += tmp_1322;
   tmp_1312 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1321;
   std::complex<double> tmp_1323;
   std::complex<double> tmp_1324;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1324 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1323 += tmp_1324;
   tmp_1312 += (std::complex<double>(0,0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1323;
   result += (std::complex<double>(0,-1)) * tmp_1312;

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1325;
   std::complex<double> tmp_1326;
   std::complex<double> tmp_1327;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1327 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1326 += tmp_1327;
   tmp_1325 += (std::complex<double>(0,0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1326;
   std::complex<double> tmp_1328;
   std::complex<double> tmp_1329;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1329 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1328 += tmp_1329;
   tmp_1325 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1328;
   std::complex<double> tmp_1330;
   std::complex<double> tmp_1331;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1331 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1330 += tmp_1331;
   tmp_1325 += (std::complex<double>(0,-0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1330;
   std::complex<double> tmp_1332;
   std::complex<double> tmp_1333;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1333 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1332 += tmp_1333;
   tmp_1325 += (std::complex<double>(0,-0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1332;
   std::complex<double> tmp_1334;
   std::complex<double> tmp_1335;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1335 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1334 += tmp_1335;
   tmp_1325 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1334;
   std::complex<double> tmp_1336;
   std::complex<double> tmp_1337;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1337 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1336 += tmp_1337;
   tmp_1325 += (std::complex<double>(0,0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1336;
   result += (std::complex<double>(0,-1)) * tmp_1325;

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1338;
   std::complex<double> tmp_1339;
   std::complex<double> tmp_1340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1340 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1339 += tmp_1340;
   tmp_1338 += (std::complex<double>(0,-0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1339;
   std::complex<double> tmp_1341;
   std::complex<double> tmp_1342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1342 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1341 += tmp_1342;
   tmp_1338 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1341;
   std::complex<double> tmp_1343;
   std::complex<double> tmp_1344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1344 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1343 += tmp_1344;
   tmp_1338 += (std::complex<double>(0,0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1343;
   std::complex<double> tmp_1345;
   std::complex<double> tmp_1346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1346 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1345 += tmp_1346;
   tmp_1338 += (std::complex<double>(0,0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1345;
   std::complex<double> tmp_1347;
   std::complex<double> tmp_1348;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1348 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1347 += tmp_1348;
   tmp_1338 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1347;
   std::complex<double> tmp_1349;
   std::complex<double> tmp_1350;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1350 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1349 += tmp_1350;
   tmp_1338 += (std::complex<double>(0,-0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1349;
   result += (std::complex<double>(0,-1)) * tmp_1338;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1351;
   std::complex<double> tmp_1352;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1353;
      std::complex<double> tmp_1354;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1354 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1353 += tmp_1354;
      tmp_1352 += (Conj(ZD(gI2,j2))) * tmp_1353;
   }
   tmp_1351 += tmp_1352;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(
      MuD))*KroneckerDelta(0,gO2)) * tmp_1351;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1355;
   std::complex<double> tmp_1356;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1357;
      std::complex<double> tmp_1358;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1358 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1357 += tmp_1358;
      tmp_1356 += (Conj(ZE(gI2,j2))) * tmp_1357;
   }
   tmp_1355 += tmp_1356;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(
      MuD))*KroneckerDelta(0,gO2)) * tmp_1355;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1359;
   std::complex<double> tmp_1360;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1361;
      std::complex<double> tmp_1362;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1362 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1361 += tmp_1362;
      tmp_1360 += (Conj(ZU(gI2,j2))) * tmp_1361;
   }
   tmp_1359 += tmp_1360;
   result += (-0.5*(1.4142135623730951*vS*Conj(LamSU) - vT*Conj(LamTU) + 2*Conj
      (MuU))*KroneckerDelta(1,gO2)) * tmp_1359;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhVZRh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(
      KroneckerDelta(0,gO2)*ZHR(gI2,0) - KroneckerDelta(1,gO2)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhSRdpHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(-1.4142135623730951*KroneckerDelta(1,gO2)*(2*LamSD*vu*Conj(
      LamSU)*ZP(gI2,0) + LamTD*Conj(LamTU)*(vu*ZP(gI2,0) + 2*vd*ZP(gI2,1))) -
      KroneckerDelta(0,gO2)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(
      LamTD) + Sqr(g2))*ZP(gI2,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI2,1) + 4*g2
      *MDWBT*ZP(gI2,2) - 2*vT*AbsSqr(LamTD)*ZP(gI2,2) - 2.8284271247461903*LamTD*
      vS*Conj(LamSD)*ZP(gI2,2) - 4*LamTD*Conj(MuD)*ZP(gI2,2) + 2*vT*Sqr(g2)*ZP(gI2
      ,2) + 2*vT*AbsSqr(LamTD)*ZP(gI2,3) - 4*MuD*Conj(LamTD)*ZP(gI2,3) -
      2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gI2,3) + 4*g2*Conj(MDWBT)*ZP(gI2,
      3) - 2*vT*Sqr(g2)*ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSRdpSRdp(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*
      Sqr(g2))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamSD) -
      10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2)) + 10*(-(KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) + KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)
      *(-2*AbsSqr(LamTD) + Sqr(g2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSRumSRum(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*
      Sqr(g2))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-20*AbsSqr(LamSU) -
      10*AbsSqr(LamTU) + 3*Sqr(g1) + 5*Sqr(g2)) - 10*(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2) - KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*
      (-2*AbsSqr(LamTU) + Sqr(g2))));

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

std::complex<double> CLASSNAME::CpconjUHpmconjRhSRum(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = 0.25*(-2.8284271247461903*LamSU*vd*Conj(LamSD)*KroneckerDelta(1,gO2
      )*ZHR(gI1,0) - 1.4142135623730951*LamTU*Conj(LamTD)*(2*vu*KroneckerDelta(0,
      gO2) + vd*KroneckerDelta(1,gO2))*ZHR(gI1,0) - (1.4142135623730951*vd*
      KroneckerDelta(0,gO2)*Sqr(g2) + 1.4142135623730951*vu*KroneckerDelta(1,gO2)*
      (-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2)) + 2*((-((2*MuU +
      1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT + 2*Conj(
      MDWBT)))*KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2)*(-2*g2*MDWBT - vT*
      AbsSqr(LamTU) + 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(MuU)
      + vT*Sqr(g2))))*ZHR(gI1,1));

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

std::complex<double> CLASSNAME::CpconjUHpmconjRhHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(LamTU)*KroneckerDelta(2,gO2)*Mu*ZHR(gI1,1)*ZP(gI2,0)) + Conj
      (LamTD)*KroneckerDelta(1,gO2)*Mu*ZHR(gI1,0)*ZP(gI2,3);

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

   std::complex<double> tmp_1363;
   std::complex<double> tmp_1364;
   std::complex<double> tmp_1365;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1366;
      std::complex<double> tmp_1367;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1368;
         std::complex<double> tmp_1369;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1369 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1368 += tmp_1369;
         tmp_1367 += (Conj(ZV(gI2,j2))) * tmp_1368;
      }
      tmp_1366 += tmp_1367;
      tmp_1365 += (ZV(gI1,j3)) * tmp_1366;
   }
   tmp_1364 += tmp_1365;
   tmp_1363 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1364;
   tmp_1363 += std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1363 += std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1363 += std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1363 += std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1363 += std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*KroneckerDelta
      (2,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1363 += std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   result += (std::complex<double>(0,-1)) * tmp_1363;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1370;
   std::complex<double> tmp_1371;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1372;
      std::complex<double> tmp_1373;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1373 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1372 += tmp_1373;
      tmp_1371 += (ZUL(gI1,j2)) * tmp_1372;
   }
   tmp_1370 += tmp_1371;
   result += (KroneckerDelta(0,gO2)) * tmp_1370;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1374;
   std::complex<double> tmp_1375;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1376;
      std::complex<double> tmp_1377;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1377 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1376 += tmp_1377;
      tmp_1375 += (Conj(ZDL(gI2,j2))) * tmp_1376;
   }
   tmp_1374 += tmp_1375;
   result += (KroneckerDelta(1,gO1)) * tmp_1374;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1378;
   std::complex<double> tmp_1379;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1379 += Conj(Ye(j1,gI1))*ZER(gI2,j1);
   }
   tmp_1378 += tmp_1379;
   result += (KroneckerDelta(0,gO2)) * tmp_1378;

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

   std::complex<double> tmp_1380;
   std::complex<double> tmp_1381;
   std::complex<double> tmp_1382;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1382 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1381 += tmp_1382;
   tmp_1380 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1381;
   std::complex<double> tmp_1383;
   std::complex<double> tmp_1384;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1385;
      std::complex<double> tmp_1386;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1387;
         std::complex<double> tmp_1388;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1388 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1387 += tmp_1388;
         tmp_1386 += (Conj(ZE(gI2,j2))) * tmp_1387;
      }
      tmp_1385 += tmp_1386;
      tmp_1384 += (ZV(gI1,j3)) * tmp_1385;
   }
   tmp_1383 += tmp_1384;
   tmp_1380 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1383;
   std::complex<double> tmp_1389;
   std::complex<double> tmp_1390;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1390 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1389 += tmp_1390;
   tmp_1380 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1389;
   std::complex<double> tmp_1391;
   std::complex<double> tmp_1392;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1393;
      std::complex<double> tmp_1394;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1394 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1393 += tmp_1394;
      tmp_1392 += (ZV(gI1,j2)) * tmp_1393;
   }
   tmp_1391 += tmp_1392;
   tmp_1380 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1391;
   std::complex<double> tmp_1395;
   std::complex<double> tmp_1396;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1396 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1395 += tmp_1396;
   tmp_1380 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1395;
   std::complex<double> tmp_1397;
   std::complex<double> tmp_1398;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1398 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1397 += tmp_1398;
   tmp_1380 += (std::complex<double>(0,-1)*g2*Conj(MDWBT)*KroneckerDelta(2,gO2)
      ) * tmp_1397;
   std::complex<double> tmp_1399;
   std::complex<double> tmp_1400;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1400 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1399 += tmp_1400;
   tmp_1380 += (std::complex<double>(0,-1)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1399;
   std::complex<double> tmp_1401;
   std::complex<double> tmp_1402;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1402 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1401 += tmp_1402;
   tmp_1380 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1401;
   result += (std::complex<double>(0,-1)) * tmp_1380;

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

   std::complex<double> tmp_1403;
   std::complex<double> tmp_1404;
   std::complex<double> tmp_1405;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1405 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1404 += tmp_1405;
   tmp_1403 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1404;
   std::complex<double> tmp_1406;
   std::complex<double> tmp_1407;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1407 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1406 += tmp_1407;
   tmp_1403 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1406;
   std::complex<double> tmp_1408;
   std::complex<double> tmp_1409;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1409 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1408 += tmp_1409;
   tmp_1403 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1408;
   std::complex<double> tmp_1410;
   std::complex<double> tmp_1411;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1412;
      std::complex<double> tmp_1413;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1414;
         std::complex<double> tmp_1415;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1415 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1414 += tmp_1415;
         tmp_1413 += (ZD(gI1,3 + j2)) * tmp_1414;
      }
      tmp_1412 += tmp_1413;
      tmp_1411 += (Conj(ZD(gI2,3 + j3))) * tmp_1412;
   }
   tmp_1410 += tmp_1411;
   tmp_1403 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1410;
   std::complex<double> tmp_1416;
   std::complex<double> tmp_1417;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1417 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1416 += tmp_1417;
   tmp_1403 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1416;
   std::complex<double> tmp_1418;
   std::complex<double> tmp_1419;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1419 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1418 += tmp_1419;
   tmp_1403 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1418;
   std::complex<double> tmp_1420;
   std::complex<double> tmp_1421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1421 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1420 += tmp_1421;
   tmp_1403 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1420;
   std::complex<double> tmp_1422;
   std::complex<double> tmp_1423;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1424;
      std::complex<double> tmp_1425;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1426;
         std::complex<double> tmp_1427;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1427 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1426 += tmp_1427;
         tmp_1425 += (Conj(ZD(gI2,j2))) * tmp_1426;
      }
      tmp_1424 += tmp_1425;
      tmp_1423 += (ZD(gI1,j3)) * tmp_1424;
   }
   tmp_1422 += tmp_1423;
   tmp_1403 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1422;
   std::complex<double> tmp_1428;
   std::complex<double> tmp_1429;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1429 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1428 += tmp_1429;
   tmp_1403 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1428;
   std::complex<double> tmp_1430;
   std::complex<double> tmp_1431;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1431 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1430 += tmp_1431;
   tmp_1403 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1430;
   result += (std::complex<double>(0,-1)) * tmp_1403;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1432;
   std::complex<double> tmp_1433;
   std::complex<double> tmp_1434;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1434 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1433 += tmp_1434;
   tmp_1432 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1433;
   std::complex<double> tmp_1435;
   std::complex<double> tmp_1436;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1436 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1435 += tmp_1436;
   tmp_1432 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1435;
   std::complex<double> tmp_1437;
   std::complex<double> tmp_1438;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1438 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1437 += tmp_1438;
   tmp_1432 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1437;
   std::complex<double> tmp_1439;
   std::complex<double> tmp_1440;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1441;
      std::complex<double> tmp_1442;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1443;
         std::complex<double> tmp_1444;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1444 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1443 += tmp_1444;
         tmp_1442 += (ZE(gI1,3 + j2)) * tmp_1443;
      }
      tmp_1441 += tmp_1442;
      tmp_1440 += (Conj(ZE(gI2,3 + j3))) * tmp_1441;
   }
   tmp_1439 += tmp_1440;
   tmp_1432 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1439;
   std::complex<double> tmp_1445;
   std::complex<double> tmp_1446;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1446 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1445 += tmp_1446;
   tmp_1432 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1445;
   std::complex<double> tmp_1447;
   std::complex<double> tmp_1448;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1448 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1447 += tmp_1448;
   tmp_1432 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1447;
   std::complex<double> tmp_1449;
   std::complex<double> tmp_1450;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1450 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1449 += tmp_1450;
   tmp_1432 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1449;
   std::complex<double> tmp_1451;
   std::complex<double> tmp_1452;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1452 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1451 += tmp_1452;
   tmp_1432 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1451;
   std::complex<double> tmp_1453;
   std::complex<double> tmp_1454;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1454 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1453 += tmp_1454;
   tmp_1432 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1453;
   result += (std::complex<double>(0,-1)) * tmp_1432;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1455;
   std::complex<double> tmp_1456;
   std::complex<double> tmp_1457;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1457 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1456 += tmp_1457;
   tmp_1455 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1456;
   std::complex<double> tmp_1458;
   std::complex<double> tmp_1459;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1459 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1458 += tmp_1459;
   tmp_1455 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1458;
   std::complex<double> tmp_1460;
   std::complex<double> tmp_1461;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1461 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1460 += tmp_1461;
   tmp_1455 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1460;
   std::complex<double> tmp_1462;
   std::complex<double> tmp_1463;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1464;
      std::complex<double> tmp_1465;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1466;
         std::complex<double> tmp_1467;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1467 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1466 += tmp_1467;
         tmp_1465 += (Conj(ZU(gI2,j2))) * tmp_1466;
      }
      tmp_1464 += tmp_1465;
      tmp_1463 += (ZU(gI1,j3)) * tmp_1464;
   }
   tmp_1462 += tmp_1463;
   tmp_1455 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1462;
   std::complex<double> tmp_1468;
   std::complex<double> tmp_1469;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1469 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1468 += tmp_1469;
   tmp_1455 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1468;
   std::complex<double> tmp_1470;
   std::complex<double> tmp_1471;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1471 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1470 += tmp_1471;
   tmp_1455 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1470;
   std::complex<double> tmp_1472;
   std::complex<double> tmp_1473;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1473 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1472 += tmp_1473;
   tmp_1455 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1472;
   std::complex<double> tmp_1474;
   std::complex<double> tmp_1475;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1476;
      std::complex<double> tmp_1477;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1478;
         std::complex<double> tmp_1479;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1479 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1478 += tmp_1479;
         tmp_1477 += (ZU(gI1,3 + j2)) * tmp_1478;
      }
      tmp_1476 += tmp_1477;
      tmp_1475 += (Conj(ZU(gI2,3 + j3))) * tmp_1476;
   }
   tmp_1474 += tmp_1475;
   tmp_1455 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1474;
   std::complex<double> tmp_1480;
   std::complex<double> tmp_1481;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1481 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1480 += tmp_1481;
   tmp_1455 += (std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1480;
   std::complex<double> tmp_1482;
   std::complex<double> tmp_1483;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1483 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1482 += tmp_1483;
   tmp_1455 += (std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1482;
   result += (std::complex<double>(0,-1)) * tmp_1455;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1484;
   std::complex<double> tmp_1485;
   std::complex<double> tmp_1486;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1486 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1485 += tmp_1486;
   tmp_1484 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1485;
   std::complex<double> tmp_1487;
   std::complex<double> tmp_1488;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1489;
      std::complex<double> tmp_1490;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1490 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1489 += tmp_1490;
      tmp_1488 += (Conj(ZD(gI2,j2))) * tmp_1489;
   }
   tmp_1487 += tmp_1488;
   tmp_1484 += (std::complex<double>(0,1)*Conj(Mu)*KroneckerDelta(0,gO2)) *
      tmp_1487;
   std::complex<double> tmp_1491;
   std::complex<double> tmp_1492;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1493;
      std::complex<double> tmp_1494;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1495;
         std::complex<double> tmp_1496;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1496 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1495 += tmp_1496;
         tmp_1494 += (ZU(gI1,3 + j2)) * tmp_1495;
      }
      tmp_1493 += tmp_1494;
      tmp_1492 += (Conj(ZD(gI2,3 + j3))) * tmp_1493;
   }
   tmp_1491 += tmp_1492;
   tmp_1484 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_1491;
   std::complex<double> tmp_1497;
   std::complex<double> tmp_1498;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1499;
      std::complex<double> tmp_1500;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1501;
         std::complex<double> tmp_1502;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1502 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1501 += tmp_1502;
         tmp_1500 += (Conj(ZD(gI2,j2))) * tmp_1501;
      }
      tmp_1499 += tmp_1500;
      tmp_1498 += (ZU(gI1,j3)) * tmp_1499;
   }
   tmp_1497 += tmp_1498;
   tmp_1484 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1497;
   std::complex<double> tmp_1503;
   std::complex<double> tmp_1504;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1504 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1503 += tmp_1504;
   tmp_1484 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1503;
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
   tmp_1484 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1505;
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
   tmp_1484 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_1509;
   std::complex<double> tmp_1515;
   std::complex<double> tmp_1516;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1517;
      std::complex<double> tmp_1518;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1519;
         std::complex<double> tmp_1520;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1520 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1519 += tmp_1520;
         tmp_1518 += (Conj(ZD(gI2,j2))) * tmp_1519;
      }
      tmp_1517 += tmp_1518;
      tmp_1516 += (ZU(gI1,j3)) * tmp_1517;
   }
   tmp_1515 += tmp_1516;
   tmp_1484 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_1515;
   std::complex<double> tmp_1521;
   std::complex<double> tmp_1522;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1522 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1521 += tmp_1522;
   tmp_1484 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1521;
   std::complex<double> tmp_1523;
   std::complex<double> tmp_1524;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1524 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1523 += tmp_1524;
   tmp_1484 += (std::complex<double>(0,-1)*g2*Conj(MDWBT)*KroneckerDelta(2,gO2)
      ) * tmp_1523;
   std::complex<double> tmp_1525;
   std::complex<double> tmp_1526;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1526 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1525 += tmp_1526;
   tmp_1484 += (std::complex<double>(0,-1)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1525;
   std::complex<double> tmp_1527;
   std::complex<double> tmp_1528;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1528 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1527 += tmp_1528;
   tmp_1484 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1527;
   result += (std::complex<double>(0,-1)) * tmp_1484;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSRdpRh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(2*(((2*MuD + 1.4142135623730951*LamSD*vS + LamTD*vT)*Conj(
      LamTD) - g2*(g2*vT + 2*Conj(MDWBT)))*KroneckerDelta(2,gO2) + KroneckerDelta(
      3,gO2)*(-2*g2*MDWBT - vT*AbsSqr(LamTD) + 1.4142135623730951*LamTD*vS*Conj(
      LamSD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZHR(gI2,0) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*(vu*Sqr(g2)*ZHR(gI2,0) + 2*LamTU*vd*Conj(LamTD)*ZHR(
      gI2,1)) - 1.4142135623730951*KroneckerDelta(0,gO2)*(vd*(-2*AbsSqr(LamSD) +
      AbsSqr(LamTD) + Sqr(g2))*ZHR(gI2,0) + vu*(2*LamSU*Conj(LamSD) + LamTU*Conj(
      LamTD))*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmSRumAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*Conj(Mu)*(1.4142135623730951*LamTU*
      KroneckerDelta(3,gO2)*ZA(gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*
      LamSU*ZA(gI2,2) + LamTU*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmSRumhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Conj(Mu)*(-1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*ZH(
      gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*LamSU*ZH(gI2,2) + LamTU*
      ZH(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSRdpAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*Mu*(1.4142135623730951*Conj(LamSD)*
      KroneckerDelta(1,gO2)*ZA(gI2,2) + Conj(LamTD)*(1.4142135623730951*
      KroneckerDelta(2,gO2)*ZA(gI2,1) - KroneckerDelta(1,gO2)*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSRdphh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Mu*(-1.4142135623730951*Conj(LamSD)*KroneckerDelta(1,gO2)*ZH(
      gI2,2) + Conj(LamTD)*(1.4142135623730951*KroneckerDelta(2,gO2)*ZH(gI2,1) +
      KroneckerDelta(1,gO2)*ZH(gI2,3)));

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

std::complex<double> CLASSNAME::CpSRdpconjSRdpVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpconjSRdpSRdp() const
{
   double result = 0.0;

   result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpconjSRumSRum() const
{
   double result = 0.0;

   result = 0.25*(0.6*Sqr(g1) + Sqr(g2));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpconjVWmVWm() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjSRdpVPSRdp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjSRdpVZSRdp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0)) + (3*Sqr(g1)
      - 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.05*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,
      0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSD*ZA(gI2,2) - 2*LamTD*ZA(gI2,3))) + Conj(LamSD)*(1.4142135623730951*
      LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*((-4*AbsSqr(LamSD) - 2*AbsSqr(LamTD) + 0.6*Sqr(g1) + Sqr(g2))*
      ZP(gI1,0)*ZP(gI2,0) - (0.6*Sqr(g1) + Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 4*AbsSqr
      (LamTD)*ZP(gI1,2)*ZP(gI2,2) + 2*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 2*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdphhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,
      0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSD*ZH(gI2,2) - 2*LamTD*ZH(gI2,3))) + Conj(LamSD)*(1.4142135623730951*
      LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSD*ZH(gI2,2) +
      1.4142135623730951*LamTD*ZH(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjHpmRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(-1.4142135623730951*ZHR(gI2,1)*(2*LamSU*vu*Conj(LamSD)*ZP(gI1
      ,0) + LamTU*Conj(LamTD)*(vu*ZP(gI1,0) + 2*vd*ZP(gI1,1))) - ZHR(gI2,0)*(
      1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,0)
      + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1) - 2*vT*AbsSqr(LamTD)*ZP(gI1,2) -
      4*MuD*Conj(LamTD)*ZP(gI1,2) - 2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gI1
      ,2) + 4*g2*Conj(MDWBT)*ZP(gI1,2) + 2*vT*Sqr(g2)*ZP(gI1,2) + 4*g2*MDWBT*ZP(
      gI1,3) + 2*vT*AbsSqr(LamTD)*ZP(gI1,3) - 2.8284271247461903*LamTD*vS*Conj(
      LamSD)*ZP(gI1,3) - 4*LamTD*Conj(MuD)*ZP(gI1,3) - 2*vT*Sqr(g2)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpChiCha1PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(LamSD)*UM1(gI2,1)*ZN2(gI1,0)) + Conj(LamTD)*(
      0.7071067811865475*UM1(gI2,1)*ZN2(gI1,1) - UM1(gI2,0)*ZN2(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpChiCha1PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*Conj(UP1(gI2,1))*(0.7745966692414834*g1*Conj(
      ZN1(gI1,0)) + g2*Conj(ZN1(gI1,1))) - g2*Conj(UP1(gI2,0))*Conj(ZN1(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*Mu*(1.4142135623730951*Conj(LamSD)*ZA(
      gI2,2)*ZP(gI1,1) + Conj(LamTD)*(-(ZA(gI2,3)*ZP(gI1,1)) + 1.4142135623730951*
      ZA(gI2,1)*ZP(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Mu*(-1.4142135623730951*Conj(LamSD)*ZH(gI2,2)*ZP(gI1,1) + Conj(
      LamTD)*(ZH(gI2,3)*ZP(gI1,1) + 1.4142135623730951*ZH(gI2,1)*ZP(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1529;
   std::complex<double> tmp_1530;
   std::complex<double> tmp_1531;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1531 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1530 += tmp_1531;
   tmp_1529 += (-0.6*Sqr(g1) + 3*Sqr(g2)) * tmp_1530;
   std::complex<double> tmp_1532;
   std::complex<double> tmp_1533;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1533 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1532 += tmp_1533;
   tmp_1529 += (-1.2*Sqr(g1)) * tmp_1532;
   result += (0.08333333333333333) * tmp_1529;

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1534;
   std::complex<double> tmp_1535;
   std::complex<double> tmp_1536;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1536 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1535 += tmp_1536;
   tmp_1534 += (0.6*Sqr(g1) + Sqr(g2)) * tmp_1535;
   std::complex<double> tmp_1537;
   std::complex<double> tmp_1538;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1538 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1537 += tmp_1538;
   tmp_1534 += (-1.2*Sqr(g1)) * tmp_1537;
   result += (0.25) * tmp_1534;

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1539;
   std::complex<double> tmp_1540;
   std::complex<double> tmp_1541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1541 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1540 += tmp_1541;
   tmp_1539 += (-0.6*Sqr(g1) - 3*Sqr(g2)) * tmp_1540;
   std::complex<double> tmp_1542;
   std::complex<double> tmp_1543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1543 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1542 += tmp_1543;
   tmp_1539 += (2.4*Sqr(g1)) * tmp_1542;
   result += (0.08333333333333333) * tmp_1539;

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjSeSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1544;
   std::complex<double> tmp_1545;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1546;
      std::complex<double> tmp_1547;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1547 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1546 += tmp_1547;
      tmp_1545 += (Conj(ZV(gI2,j2))) * tmp_1546;
   }
   tmp_1544 += tmp_1545;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(
      MuD))) * tmp_1544;

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjSdSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1548;
   std::complex<double> tmp_1549;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1550;
      std::complex<double> tmp_1551;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1551 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1550 += tmp_1551;
      tmp_1549 += (Conj(ZU(gI2,j2))) * tmp_1550;
   }
   tmp_1548 += tmp_1549;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(
      MuD))) * tmp_1548;

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpconjVWmRh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*ZHR(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpSRdpAh(unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*((1.5491933384829668*g1*MDBS +
      1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(LamSD) - 1.4142135623730951*
      LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903
      *LamSD*Conj(MuD))*ZA(gI2,2) + (2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj
      (LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(
      MDWBT) - 2*LamTD*Conj(MuD))*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRdpSRdphh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gI2,0) + vu*
      (-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2,1) - 1.5491933384829668*g1*MDBS*ZH(gI2,2) -
      4*vS*AbsSqr(LamSD)*ZH(gI2,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gI2,2)
      + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gI2,2) + 1.4142135623730951*
      LamSD*vT*Conj(LamTD)*ZH(gI2,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2)
      - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gI2,2) - 2*g2*MDWBT*ZH(gI2,3) - 2*vT
      *AbsSqr(LamTD)*ZH(gI2,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gI2,3)
      + 2*MuD*Conj(LamTD)*ZH(gI2,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(
      gI2,3) - 2*g2*Conj(MDWBT)*ZH(gI2,3) + 2*LamTD*Conj(MuD)*ZH(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpSRumconjSRumconjSRdpSRdp() const
{
   double result = 0.0;

   result = 0.25*(0.6*Sqr(g1) + Sqr(g2));

   return result;
}

double CLASSNAME::CpSRumconjSRumconjSRumSRum() const
{
   double result = 0.0;

   result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

double CLASSNAME::CpSRumconjSRumconjVWmVWm() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjSRumVPSRum() const
{
   double result = 0.0;

   result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjSRumVZSRum() const
{
   double result = 0.0;

   result = 0.5*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*((3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) +
      5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumRhHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(-2.8284271247461903*LamSD*vd*Conj(LamSU)*ZHR(gI1,0)*ZP(gI2,1)
      - ZHR(gI1,1)*(1.4142135623730951*vd*Sqr(g2)*ZP(gI2,0) + 1.4142135623730951*
      vu*(-2*AbsSqr(LamSU) + Sqr(g2))*ZP(gI2,1) + 4*g2*MDWBT*ZP(gI2,2) -
      2.8284271247461903*LamTU*vS*Conj(LamSU)*ZP(gI2,2) - 4*LamTU*Conj(MuU)*ZP(gI2
      ,2) + 2*vT*Sqr(g2)*ZP(gI2,2) + 4*g2*Conj(MDWBT)*ZP(gI2,3) - 2*vT*Sqr(g2)*ZP(
      gI2,3)) - Conj(LamTU)*(1.4142135623730951*LamTD*ZHR(gI1,0)*(2*vu*ZP(gI2,0) +
      vd*ZP(gI2,1)) + ZHR(gI1,1)*(1.4142135623730951*LamTU*vu*ZP(gI2,1) - 2*(
      LamTU*vT*ZP(gI2,2) + (2*MuU + 1.4142135623730951*LamSU*vS - LamTU*vT)*ZP(gI2
      ,3)))));

   return result;
}

double CLASSNAME::CpSRumconjSRumconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.25*KroneckerDelta(gI1,gI2)*(-0.6*Sqr(g1) + Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (-20*AbsSqr(
      LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) - 5*(Conj(LamSU)*(
      1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(4*LamSU*ZA(gI2,2)
      + 1.4142135623730951*LamTU*ZA(gI2,3))) + Conj(LamTU)*(1.4142135623730951*
      LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSU*ZA(gI2,2) +
      2*LamTU*ZA(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(-((0.6*Sqr(g1) + Sqr(g2))*ZP(gI1,0)*ZP(gI2,0)) + (-4*AbsSqr(
      LamSU) - 2*AbsSqr(LamTU) + 0.6*Sqr(g1) + Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 2*
      Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 4*AbsSqr(LamTU)*ZP(gI1,3)*ZP(gI2,3) + 2*Sqr(g2
      )*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (-20*AbsSqr(
      LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 5*(Conj(LamSU)*(
      1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(4*LamSU*ZH(gI2,2)
      + 1.4142135623730951*LamTU*ZH(gI2,3))) + Conj(LamTU)*(1.4142135623730951*
      LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951*LamSU*ZH(gI2,2) +
      2*LamTU*ZH(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumChiCha2PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(LamSU)*UP2(gI2,1)*ZN2(gI1,0) + Conj(LamTU)*(0.7071067811865475
      *UP2(gI2,1)*ZN2(gI1,1) + UP2(gI2,0)*ZN2(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumChiCha2PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*Conj(UM2(gI2,1))*(0.7745966692414834*g1*Conj(ZN1
      (gI1,0)) + g2*Conj(ZN1(gI1,1))) - g2*Conj(UM2(gI2,0))*Conj(ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*Mu*(1.4142135623730951*Conj(LamSU)*ZA(
      gI2,2)*ZP(gI1,0) + Conj(LamTU)*(ZA(gI2,3)*ZP(gI1,0) + 1.4142135623730951*ZA(
      gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Mu*(1.4142135623730951*Conj(LamSU)*ZH(gI2,2)*ZP(gI1,0) + Conj(
      LamTU)*(ZH(gI2,3)*ZP(gI1,0) - 1.4142135623730951*ZH(gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1552;
   std::complex<double> tmp_1553;
   std::complex<double> tmp_1554;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1554 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1553 += tmp_1554;
   tmp_1552 += (0.6*Sqr(g1) - 3*Sqr(g2)) * tmp_1553;
   std::complex<double> tmp_1555;
   std::complex<double> tmp_1556;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1556 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1555 += tmp_1556;
   tmp_1552 += (1.2*Sqr(g1)) * tmp_1555;
   result += (0.08333333333333333) * tmp_1552;

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1557;
   std::complex<double> tmp_1558;
   std::complex<double> tmp_1559;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1559 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1558 += tmp_1559;
   tmp_1557 += (-0.6*Sqr(g1) - Sqr(g2)) * tmp_1558;
   std::complex<double> tmp_1560;
   std::complex<double> tmp_1561;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1561 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1560 += tmp_1561;
   tmp_1557 += (1.2*Sqr(g1)) * tmp_1560;
   result += (0.25) * tmp_1557;

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1562;
   std::complex<double> tmp_1563;
   std::complex<double> tmp_1564;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1564 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1563 += tmp_1564;
   tmp_1562 += (0.6*Sqr(g1) + 3*Sqr(g2)) * tmp_1563;
   std::complex<double> tmp_1565;
   std::complex<double> tmp_1566;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1566 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1565 += tmp_1566;
   tmp_1562 += (-2.4*Sqr(g1)) * tmp_1565;
   result += (0.08333333333333333) * tmp_1562;

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1567;
   std::complex<double> tmp_1568;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1569;
      std::complex<double> tmp_1570;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1570 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1569 += tmp_1570;
      tmp_1568 += (Conj(ZD(gI2,j2))) * tmp_1569;
   }
   tmp_1567 += tmp_1568;
   result += (0.5*(-1.4142135623730951*vS*Conj(LamSU) - vT*Conj(LamTU) - 2*Conj
      (MuU))) * tmp_1567;

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumVWmRh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumSRumAh(unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU) - 7.0710678118654755*LamSU
      *vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSU
      *Conj(MuU))*ZA(gI2,2) + 5*(2*g2*MDWBT - 1.4142135623730951*LamTU*vS*Conj(
      LamSU) + 2*MuU*Conj(LamTU) + 1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*
      Conj(MDWBT) - 2*LamTU*Conj(MuU))*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjSRumSRumhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*(vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2,0) + vu*(-4*AbsSqr(LamTU)
      + 0.6*Sqr(g1) - Sqr(g2))*ZH(gI2,1) + 1.5491933384829668*g1*MDBS*ZH(gI2,2) -
      4*vS*AbsSqr(LamSU)*ZH(gI2,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gI2,2)
      - 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gI2,2) - 1.4142135623730951*
      LamSU*vT*Conj(LamTU)*ZH(gI2,2) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2)
      - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gI2,2) + 2*g2*MDWBT*ZH(gI2,3) - 2*vT
      *AbsSqr(LamTU)*ZH(gI2,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gI2,3)
      - 2*MuU*Conj(LamTU)*ZH(gI2,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(
      gI2,3) + 2*g2*Conj(MDWBT)*ZH(gI2,3) - 2*LamTU*Conj(MuU)*ZH(gI2,3));

   return result;
}

double CLASSNAME::CpsigmaOsigmaOphiOphiO() const
{
   double result = 0.0;

   result = -6*Sqr(g3);

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOVGsigmaO() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1571;
   std::complex<double> tmp_1572;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1572 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1571 += tmp_1572;
   std::complex<double> tmp_1573;
   std::complex<double> tmp_1574;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1574 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1573 += tmp_1574;
   tmp_1571 += (-1) * tmp_1573;
   result += (std::complex<double>(0,-1)*g3*(MDGoc - Conj(MDGoc))) * tmp_1571;

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1575;
   std::complex<double> tmp_1576;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1576 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1575 += tmp_1576;
   std::complex<double> tmp_1577;
   std::complex<double> tmp_1578;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1578 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1577 += tmp_1578;
   tmp_1575 += (-1) * tmp_1577;
   result += (std::complex<double>(0,-1)*g3*(MDGoc - Conj(MDGoc))) * tmp_1575;

   return result;
}

double CLASSNAME::CpsigmaObarGluGluPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpsigmaObarGluGluPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpphiOphiOsigmaOsigmaO() const
{
   double result = 0.0;

   result = -6*Sqr(g3);

   return result;
}

std::complex<double> CLASSNAME::CpphiOVGphiO() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpphiOconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1579;
   std::complex<double> tmp_1580;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1580 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1579 += tmp_1580;
   std::complex<double> tmp_1581;
   std::complex<double> tmp_1582;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1582 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1581 += tmp_1582;
   tmp_1579 += (-1) * tmp_1581;
   result += (-(g3*(MDGoc + Conj(MDGoc)))) * tmp_1579;

   return result;
}

std::complex<double> CLASSNAME::CpphiOconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1583;
   std::complex<double> tmp_1584;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1584 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1583 += tmp_1584;
   std::complex<double> tmp_1585;
   std::complex<double> tmp_1586;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1586 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1585 += tmp_1586;
   tmp_1583 += (-1) * tmp_1585;
   result += (-(g3*(MDGoc + Conj(MDGoc)))) * tmp_1583;

   return result;
}

double CLASSNAME::CpphiObarGluGluPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpphiObarGluGluPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVGphiOphiO() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpVGsigmaOsigmaO() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

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

double CLASSNAME::CpVGVGphiOphiO() const
{
   double result = 0.0;

   result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGsigmaOsigmaO() const
{
   double result = 0.0;

   result = 16*Sqr(g3);

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

double CLASSNAME::CpVGVGconjSdSd(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGconjSuSu(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGconjSdSd(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVGconjSuSu(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVGbarGluGluPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGbarGluGluPR(unsigned , unsigned ) const
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

double CLASSNAME::CpVPconjSRdpSRdp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPconjSRumSRum() const
{
   double result = 0.0;

   result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSRdpSRdp() const
{
   std::complex<double> result;

   result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSRumSRum() const
{
   std::complex<double> result;

   result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpVPconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarCha1Cha1PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-2*g2*Conj(UP1(gI2,0))*Sin(ThetaW())*UP1(gI1,0) - Conj(UP1(gI2
      ,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UP1(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarCha1Cha1PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-2*g2*Conj(UM1(gI1,0))*Sin(ThetaW())*UM1(gI2,0) - Conj(UM1(gI1
      ,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UM1(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarCha2Cha2PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM2(gI2,0))*Sin(ThetaW())*UM2(gI1,0) + Conj(UM2(gI2,
      1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UM2(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarCha2Cha2PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP2(gI1,0))*Sin(ThetaW())*UP2(gI2,0) + Conj(UP2(gI1,
      1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UP2(gI2,1));

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

std::complex<double> CLASSNAME::CpVPVPconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*((g2*Sin(ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())))*ZP(gI1,0)*ZP(gI2,0) + (g2*Sin
      (ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr
      (g1)*Sqr(Cos(ThetaW())))*ZP(gI1,1)*ZP(gI2,1) + 20*Sqr(g2)*Sqr(Sin(ThetaW()))
      *(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpVPconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-((0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*ZP(
      gI1,0)*ZP(gI2,0)) - (0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))
      *ZP(gI1,1)*ZP(gI2,1) - 2*g2*Sin(ThetaW())*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*
      ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1587;
   std::complex<double> tmp_1588;
   std::complex<double> tmp_1589;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1589 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1588 += tmp_1589;
   tmp_1587 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1588;
   std::complex<double> tmp_1590;
   std::complex<double> tmp_1591;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1591 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1590 += tmp_1591;
   tmp_1587 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1590;
   std::complex<double> tmp_1592;
   std::complex<double> tmp_1593;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1593 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1592 += tmp_1593;
   tmp_1587 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1592;
   std::complex<double> tmp_1594;
   std::complex<double> tmp_1595;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1595 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1594 += tmp_1595;
   tmp_1587 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1594;
   result += (std::complex<double>(0,-1)) * tmp_1587;

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1596;
   std::complex<double> tmp_1597;
   std::complex<double> tmp_1598;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1598 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1597 += tmp_1598;
   tmp_1596 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Cos(ThetaW()))) *
      tmp_1597;
   std::complex<double> tmp_1599;
   std::complex<double> tmp_1600;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1600 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1599 += tmp_1600;
   tmp_1596 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Cos(ThetaW()))) *
      tmp_1599;
   std::complex<double> tmp_1601;
   std::complex<double> tmp_1602;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1602 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1601 += tmp_1602;
   tmp_1596 += (std::complex<double>(0.,0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1601;
   std::complex<double> tmp_1603;
   std::complex<double> tmp_1604;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1604 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1603 += tmp_1604;
   tmp_1596 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1603;
   result += (std::complex<double>(0,-1)) * tmp_1596;

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1605;
   std::complex<double> tmp_1606;
   std::complex<double> tmp_1607;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1607 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1606 += tmp_1607;
   tmp_1605 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1606;
   std::complex<double> tmp_1608;
   std::complex<double> tmp_1609;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1609 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1608 += tmp_1609;
   tmp_1605 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1608;
   std::complex<double> tmp_1610;
   std::complex<double> tmp_1611;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1611 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1610 += tmp_1611;
   tmp_1605 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1610;
   std::complex<double> tmp_1612;
   std::complex<double> tmp_1613;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1613 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1612 += tmp_1613;
   tmp_1605 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1612;
   result += (std::complex<double>(0,-1)) * tmp_1605;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1614;
   std::complex<double> tmp_1615;
   std::complex<double> tmp_1616;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1616 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1615 += tmp_1616;
   tmp_1614 += (-1.5491933384829668*g1*Cos(ThetaW())) * tmp_1615;
   std::complex<double> tmp_1617;
   std::complex<double> tmp_1618;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1618 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1617 += tmp_1618;
   tmp_1614 += (0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW())) *
      tmp_1617;
   result += (0.16666666666666666) * tmp_1614;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1619;
   std::complex<double> tmp_1620;
   std::complex<double> tmp_1621;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1621 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1620 += tmp_1621;
   tmp_1619 += (-1.5491933384829668*g1*Cos(ThetaW())) * tmp_1620;
   std::complex<double> tmp_1622;
   std::complex<double> tmp_1623;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1623 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1622 += tmp_1623;
   tmp_1619 += (-0.7745966692414834*g1*Cos(ThetaW()) - g2*Sin(ThetaW())) *
      tmp_1622;
   result += (0.5) * tmp_1619;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1624;
   std::complex<double> tmp_1625;
   std::complex<double> tmp_1626;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1626 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1625 += tmp_1626;
   tmp_1624 += (3.0983866769659336*g1*Cos(ThetaW())) * tmp_1625;
   std::complex<double> tmp_1627;
   std::complex<double> tmp_1628;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1628 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1627 += tmp_1628;
   tmp_1624 += (0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW())) *
      tmp_1627;
   result += (0.16666666666666666) * tmp_1624;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjVWmHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(0.7745966692414834*g1*vd*Cos(ThetaW())*ZP(gI2,0) -
      0.7745966692414834*g1*vu*Cos(ThetaW())*ZP(gI2,1) + 1.4142135623730951*g2*vT*
      Sin(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

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

double CLASSNAME::CpVZconjSRdpSRdp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjSRumSRum() const
{
   double result = 0.0;

   result = 0.5*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSRdpSRdp() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSRumSRum() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

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

std::complex<double> CLASSNAME::CpVZconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(ZHR(
      gI1,0)*ZHR(gI2,0) - ZHR(gI1,1)*ZHR(gI2,1));

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

   std::complex<double> tmp_1629;
   std::complex<double> tmp_1630;
   std::complex<double> tmp_1631;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1631 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1630 += tmp_1631;
   tmp_1629 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1630;
   std::complex<double> tmp_1632;
   std::complex<double> tmp_1633;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1633 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1632 += tmp_1633;
   tmp_1629 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1632;
   std::complex<double> tmp_1634;
   std::complex<double> tmp_1635;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1635 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1634 += tmp_1635;
   tmp_1629 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1634;
   std::complex<double> tmp_1636;
   std::complex<double> tmp_1637;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1637 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1636 += tmp_1637;
   tmp_1629 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1636;
   result += (std::complex<double>(0,-1)) * tmp_1629;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1638;
   std::complex<double> tmp_1639;
   std::complex<double> tmp_1640;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1640 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1639 += tmp_1640;
   tmp_1638 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1639;
   std::complex<double> tmp_1641;
   std::complex<double> tmp_1642;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1642 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1641 += tmp_1642;
   tmp_1638 += (std::complex<double>(0.,-0.7745966692414834)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1641;
   std::complex<double> tmp_1643;
   std::complex<double> tmp_1644;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1644 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1643 += tmp_1644;
   tmp_1638 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1643;
   std::complex<double> tmp_1645;
   std::complex<double> tmp_1646;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1646 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1645 += tmp_1646;
   tmp_1638 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1645;
   result += (std::complex<double>(0,-1)) * tmp_1638;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1647;
   std::complex<double> tmp_1648;
   std::complex<double> tmp_1649;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1649 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1648 += tmp_1649;
   tmp_1647 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1648;
   std::complex<double> tmp_1650;
   std::complex<double> tmp_1651;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1651 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1650 += tmp_1651;
   tmp_1647 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1650;
   std::complex<double> tmp_1652;
   std::complex<double> tmp_1653;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1653 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1652 += tmp_1653;
   tmp_1647 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1652;
   std::complex<double> tmp_1654;
   std::complex<double> tmp_1655;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1655 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1654 += tmp_1655;
   tmp_1647 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1654;
   result += (std::complex<double>(0,-1)) * tmp_1647;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1656;
   std::complex<double> tmp_1657;
   std::complex<double> tmp_1658;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1658 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1657 += tmp_1658;
   tmp_1656 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1657;
   std::complex<double> tmp_1659;
   std::complex<double> tmp_1660;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1660 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1659 += tmp_1660;
   tmp_1656 += (-3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1659;
   result += (0.16666666666666666) * tmp_1656;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1661;
   std::complex<double> tmp_1662;
   std::complex<double> tmp_1663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1663 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1662 += tmp_1663;
   tmp_1661 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1662;
   std::complex<double> tmp_1664;
   std::complex<double> tmp_1665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1665 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1664 += tmp_1665;
   tmp_1661 += (-(g2*Cos(ThetaW())) + 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1664;
   result += (0.5) * tmp_1661;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1666;
   std::complex<double> tmp_1667;
   std::complex<double> tmp_1668;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1668 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1667 += tmp_1668;
   tmp_1666 += (-3.0983866769659336*g1*Sin(ThetaW())) * tmp_1667;
   std::complex<double> tmp_1669;
   std::complex<double> tmp_1670;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1670 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1669 += tmp_1670;
   tmp_1666 += (3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1669;
   result += (0.16666666666666666) * tmp_1666;

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

double CLASSNAME::CpVWmconjVWmconjSRdpSRdp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjSRumSRum() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

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

std::complex<double> CLASSNAME::CpconjVWmconjRhSRum(unsigned gI1) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*ZHR(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZHR(gI1,0)*ZHR(gI2,0) + ZHR(gI1,1)*ZHR(gI2,1));

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

   std::complex<double> tmp_1671;
   std::complex<double> tmp_1672;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1672 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1671 += tmp_1672;
   result += (-0.7071067811865475*g2) * tmp_1671;

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

   std::complex<double> tmp_1673;
   std::complex<double> tmp_1674;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1674 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1673 += tmp_1674;
   result += (0.7071067811865475*g2) * tmp_1673;

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

   std::complex<double> tmp_1675;
   std::complex<double> tmp_1676;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1676 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1675 += tmp_1676;
   result += (0.5*Sqr(g2)) * tmp_1675;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1677;
   std::complex<double> tmp_1678;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1678 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1677 += tmp_1678;
   result += (0.5*Sqr(g2)) * tmp_1677;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1679;
   std::complex<double> tmp_1680;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1680 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1679 += tmp_1680;
   result += (0.5*Sqr(g2)) * tmp_1679;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1681;
   std::complex<double> tmp_1682;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1682 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1681 += tmp_1682;
   result += (0.7071067811865475*g2) * tmp_1681;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSRdpRh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*ZHR(gI2,0);

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

std::complex<double> CLASSNAME::CpbarUChibarCha1SRdpPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = Conj(UM1(gI1,1))*(-(LamSD*KroneckerDelta(0,gO2)) +
      0.7071067811865475*LamTD*KroneckerDelta(1,gO2)) - LamTD*Conj(UM1(gI1,0))*
      KroneckerDelta(2,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarCha1SRdpPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(2,gO1)*UP1(gI1,0)) - 0.1414213562373095*(
      3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*KroneckerDelta(1,gO1))*UP1
      (gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarCha2SRumPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = Conj(UP2(gI1,1))*(LamSU*KroneckerDelta(0,gO2) + 0.7071067811865475*
      LamTU*KroneckerDelta(1,gO2)) + LamTU*Conj(UP2(gI1,0))*KroneckerDelta(3,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarCha2SRumPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(3,gO1)*UM2(gI1,0)) + 0.5477225575051661*g1*
      KroneckerDelta(0,gO1)*UM2(gI1,1) + 0.7071067811865475*g2*KroneckerDelta(1,
      gO1)*UM2(gI1,1);

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

   std::complex<double> tmp_1683;
   std::complex<double> tmp_1684;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1685;
      std::complex<double> tmp_1686;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1686 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1685 += tmp_1686;
      tmp_1684 += (Conj(ZD(gI2,j2))) * tmp_1685;
   }
   tmp_1683 += tmp_1684;
   result += (-KroneckerDelta(2,gO2)) * tmp_1683;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFdSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1687;
   std::complex<double> tmp_1688;
   std::complex<double> tmp_1689;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1689 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1688 += tmp_1689;
   tmp_1687 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO1)) * tmp_1688;
   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1691 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1690 += tmp_1691;
   tmp_1687 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO1)) * tmp_1690;
   result += (std::complex<double>(0,-1)) * tmp_1687;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1692;
   std::complex<double> tmp_1693;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1694;
      std::complex<double> tmp_1695;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1695 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1694 += tmp_1695;
      tmp_1693 += (Conj(ZE(gI2,j2))) * tmp_1694;
   }
   tmp_1692 += tmp_1693;
   result += (-KroneckerDelta(2,gO2)) * tmp_1692;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1696;
   std::complex<double> tmp_1697;
   std::complex<double> tmp_1698;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1698 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1697 += tmp_1698;
   tmp_1696 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(0
      ,gO1)) * tmp_1697;
   std::complex<double> tmp_1699;
   std::complex<double> tmp_1700;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1700 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1699 += tmp_1700;
   tmp_1696 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO1)) * tmp_1699;
   result += (std::complex<double>(0,-1)) * tmp_1696;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1701;
   std::complex<double> tmp_1702;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1703;
      std::complex<double> tmp_1704;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1704 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1703 += tmp_1704;
      tmp_1702 += (Conj(ZU(gI2,j2))) * tmp_1703;
   }
   tmp_1701 += tmp_1702;
   result += (-KroneckerDelta(3,gO2)) * tmp_1701;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1705;
   std::complex<double> tmp_1706;
   std::complex<double> tmp_1707;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1707 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1706 += tmp_1707;
   tmp_1705 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO1)) * tmp_1706;
   std::complex<double> tmp_1708;
   std::complex<double> tmp_1709;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1709 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1708 += tmp_1709;
   tmp_1705 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      1,gO1)) * tmp_1708;
   result += (std::complex<double>(0,-1)) * tmp_1705;

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

   std::complex<double> tmp_1710;
   std::complex<double> tmp_1711;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1712;
      std::complex<double> tmp_1713;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1713 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1712 += tmp_1713;
      tmp_1711 += (Conj(ZDL(gI2,j2))) * tmp_1712;
   }
   tmp_1710 += tmp_1711;
   result += (-KroneckerDelta(2,gO2)) * tmp_1710;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1714;
   std::complex<double> tmp_1715;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1715 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1714 += tmp_1715;
   result += (-0.3651483716701107*g1*KroneckerDelta(0,gO1)) * tmp_1714;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1716;
   std::complex<double> tmp_1717;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1718;
      std::complex<double> tmp_1719;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1719 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1718 += tmp_1719;
      tmp_1717 += (Conj(ZEL(gI2,j2))) * tmp_1718;
   }
   tmp_1716 += tmp_1717;
   result += (-KroneckerDelta(2,gO2)) * tmp_1716;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1720;
   std::complex<double> tmp_1721;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1721 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_1720 += tmp_1721;
   result += (-1.0954451150103321*g1*KroneckerDelta(0,gO1)) * tmp_1720;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1722;
   std::complex<double> tmp_1723;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1724;
      std::complex<double> tmp_1725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1725 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1724 += tmp_1725;
      tmp_1723 += (Conj(ZUL(gI2,j2))) * tmp_1724;
   }
   tmp_1722 += tmp_1723;
   result += (-KroneckerDelta(3,gO2)) * tmp_1722;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1726;
   std::complex<double> tmp_1727;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1727 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1726 += tmp_1727;
   result += (0.7302967433402214*g1*KroneckerDelta(0,gO1)) * tmp_1726;

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

   std::complex<double> tmp_1728;
   std::complex<double> tmp_1729;
   std::complex<double> tmp_1730;
   for (unsigned gl1429 = 0; gl1429 < 2; ++gl1429) {
      tmp_1730 += Conj(UM1(gl1429,1))*UP1(gl1429,gO2);
   }
   tmp_1729 += tmp_1730;
   tmp_1728 += (std::complex<double>(0,1)*LamTD*Conj(UP2(gI1,0))*ZHR(gI2,0)) *
      tmp_1729;
   std::complex<double> tmp_1731;
   std::complex<double> tmp_1732;
   for (unsigned gl1429 = 0; gl1429 < 2; ++gl1429) {
      tmp_1732 += Conj(UM1(gl1429,0))*UP1(gl1429,gO2);
   }
   tmp_1731 += tmp_1732;
   tmp_1728 += (std::complex<double>(0,-1)*LamTU*Conj(UP2(gI1,1))*ZHR(gI2,1)) *
      tmp_1731;
   result += (std::complex<double>(0,-1)) * tmp_1728;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barCha2RhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1733;
   std::complex<double> tmp_1734;
   std::complex<double> tmp_1735;
   for (unsigned gl1432 = 0; gl1432 < 2; ++gl1432) {
      tmp_1735 += Conj(UM1(gl1432,gO1))*UP1(gl1432,1);
   }
   tmp_1734 += tmp_1735;
   tmp_1733 += (std::complex<double>(0,-1)*g2*UM2(gI1,0)*ZHR(gI2,0)) * tmp_1734
      ;
   std::complex<double> tmp_1736;
   std::complex<double> tmp_1737;
   for (unsigned gl1432 = 0; gl1432 < 2; ++gl1432) {
      tmp_1737 += Conj(UM1(gl1432,gO1))*UP1(gl1432,0);
   }
   tmp_1736 += tmp_1737;
   tmp_1733 += (std::complex<double>(0,-1)*g2*UM2(gI1,1)*ZHR(gI2,1)) * tmp_1736
      ;
   result += (std::complex<double>(0,-1)) * tmp_1733;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1738;
   std::complex<double> tmp_1739;
   std::complex<double> tmp_1740;
   for (unsigned gl1435 = 0; gl1435 < 2; ++gl1435) {
      tmp_1740 += Conj(UM1(gl1435,0))*UP1(gl1435,gO2);
   }
   tmp_1739 += tmp_1740;
   tmp_1738 += (0.7071067811865475*LamTD*Conj(UP1(gI1,1))*ZA(gI2,0)) * tmp_1739
      ;
   std::complex<double> tmp_1741;
   std::complex<double> tmp_1742;
   for (unsigned gl1435 = 0; gl1435 < 2; ++gl1435) {
      tmp_1742 += Conj(UM1(gl1435,1))*UP1(gl1435,gO2);
   }
   tmp_1741 += tmp_1742;
   tmp_1738 += (-0.7071067811865475*g2*Conj(UP1(gI1,0))*ZA(gI2,0)) * tmp_1741;
   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   for (unsigned gl1435 = 0; gl1435 < 2; ++gl1435) {
      tmp_1744 += Conj(UM1(gl1435,1))*UP1(gl1435,gO2);
   }
   tmp_1743 += tmp_1744;
   tmp_1738 += (0.7071067811865475*LamSD*Conj(UP1(gI1,1))*ZA(gI2,2)) * tmp_1743
      ;
   std::complex<double> tmp_1745;
   std::complex<double> tmp_1746;
   for (unsigned gl1435 = 0; gl1435 < 2; ++gl1435) {
      tmp_1746 += Conj(UM1(gl1435,0))*UP1(gl1435,gO2);
   }
   tmp_1745 += tmp_1746;
   tmp_1738 += (-(g2*Conj(UP1(gI1,0))*ZA(gI2,3))) * tmp_1745;
   std::complex<double> tmp_1747;
   std::complex<double> tmp_1748;
   for (unsigned gl1435 = 0; gl1435 < 2; ++gl1435) {
      tmp_1748 += Conj(UM1(gl1435,1))*UP1(gl1435,gO2);
   }
   tmp_1747 += tmp_1748;
   tmp_1738 += (-0.5*LamTD*Conj(UP1(gI1,1))*ZA(gI2,3)) * tmp_1747;
   result += (std::complex<double>(0,-1)) * tmp_1738;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1749;
   std::complex<double> tmp_1750;
   std::complex<double> tmp_1751;
   for (unsigned gl1438 = 0; gl1438 < 2; ++gl1438) {
      tmp_1751 += Conj(UM1(gl1438,gO1))*UP1(gl1438,1);
   }
   tmp_1750 += tmp_1751;
   tmp_1749 += (-0.7071067811865475*Conj(LamTD)*UM1(gI1,0)*ZA(gI2,0)) *
      tmp_1750;
   std::complex<double> tmp_1752;
   std::complex<double> tmp_1753;
   for (unsigned gl1438 = 0; gl1438 < 2; ++gl1438) {
      tmp_1753 += Conj(UM1(gl1438,gO1))*UP1(gl1438,0);
   }
   tmp_1752 += tmp_1753;
   tmp_1749 += (0.7071067811865475*g2*UM1(gI1,1)*ZA(gI2,0)) * tmp_1752;
   std::complex<double> tmp_1754;
   std::complex<double> tmp_1755;
   for (unsigned gl1438 = 0; gl1438 < 2; ++gl1438) {
      tmp_1755 += Conj(UM1(gl1438,gO1))*UP1(gl1438,1);
   }
   tmp_1754 += tmp_1755;
   tmp_1749 += (-0.7071067811865475*Conj(LamSD)*UM1(gI1,1)*ZA(gI2,2)) *
      tmp_1754;
   std::complex<double> tmp_1756;
   std::complex<double> tmp_1757;
   for (unsigned gl1438 = 0; gl1438 < 2; ++gl1438) {
      tmp_1757 += Conj(UM1(gl1438,gO1))*UP1(gl1438,0);
   }
   tmp_1756 += tmp_1757;
   tmp_1749 += (g2*UM1(gI1,0)*ZA(gI2,3)) * tmp_1756;
   std::complex<double> tmp_1758;
   std::complex<double> tmp_1759;
   for (unsigned gl1438 = 0; gl1438 < 2; ++gl1438) {
      tmp_1759 += Conj(UM1(gl1438,gO1))*UP1(gl1438,1);
   }
   tmp_1758 += tmp_1759;
   tmp_1749 += (0.5*Conj(LamTD)*UM1(gI1,1)*ZA(gI2,3)) * tmp_1758;
   result += (std::complex<double>(0,-1)) * tmp_1749;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1760;
   std::complex<double> tmp_1762;
   for (unsigned gl1441 = 0; gl1441 < 2; ++gl1441) {
      tmp_1762 += Conj(UM1(gl1441,1))*UP1(gl1441,gO2);
   }
   tmp_1760 += tmp_1762;
   std::complex<double> tmp_1761;
   std::complex<double> tmp_1763;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1764;
      std::complex<double> tmp_1765;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1765 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1764 += tmp_1765;
      tmp_1763 += (Conj(ZV(gI2,j2))) * tmp_1764;
   }
   tmp_1761 += tmp_1763;
   result += (1) * tmp_1760 * tmp_1761;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1766;
   std::complex<double> tmp_1768;
   for (unsigned gl1444 = 0; gl1444 < 2; ++gl1444) {
      tmp_1768 += Conj(UM1(gl1444,gO1))*UP1(gl1444,0);
   }
   tmp_1766 += tmp_1768;
   std::complex<double> tmp_1767;
   std::complex<double> tmp_1769;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1769 += Conj(ZV(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1767 += tmp_1769;
   result += (-g2) * tmp_1766 * tmp_1767;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1770;
   std::complex<double> tmp_1772;
   for (unsigned gl1447 = 0; gl1447 < 2; ++gl1447) {
      tmp_1772 += Conj(UM1(gl1447,1))*UP1(gl1447,gO2);
   }
   tmp_1770 += tmp_1772;
   std::complex<double> tmp_1771;
   std::complex<double> tmp_1773;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1774;
      std::complex<double> tmp_1775;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1775 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1774 += tmp_1775;
      tmp_1773 += (Conj(ZU(gI2,j2))) * tmp_1774;
   }
   tmp_1771 += tmp_1773;
   result += (1) * tmp_1770 * tmp_1771;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1776;
   std::complex<double> tmp_1778;
   for (unsigned gl1450 = 0; gl1450 < 2; ++gl1450) {
      tmp_1778 += Conj(UM1(gl1450,gO1))*UP1(gl1450,0);
   }
   tmp_1776 += tmp_1778;
   std::complex<double> tmp_1777;
   std::complex<double> tmp_1779;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1779 += Conj(ZU(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1777 += tmp_1779;
   result += (-g2) * tmp_1776 * tmp_1777;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1hhCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1780;
   std::complex<double> tmp_1781;
   std::complex<double> tmp_1782;
   for (unsigned gl1453 = 0; gl1453 < 2; ++gl1453) {
      tmp_1782 += Conj(UM1(gl1453,0))*UP1(gl1453,gO2);
   }
   tmp_1781 += tmp_1782;
   tmp_1780 += (std::complex<double>(0.,-0.7071067811865475)*LamTD*Conj(UP1(gI2
      ,1))*ZH(gI1,0)) * tmp_1781;
   std::complex<double> tmp_1783;
   std::complex<double> tmp_1784;
   for (unsigned gl1453 = 0; gl1453 < 2; ++gl1453) {
      tmp_1784 += Conj(UM1(gl1453,1))*UP1(gl1453,gO2);
   }
   tmp_1783 += tmp_1784;
   tmp_1780 += (std::complex<double>(0.,-0.7071067811865475)*g2*Conj(UP1(gI2,0)
      )*ZH(gI1,0)) * tmp_1783;
   std::complex<double> tmp_1785;
   std::complex<double> tmp_1786;
   for (unsigned gl1453 = 0; gl1453 < 2; ++gl1453) {
      tmp_1786 += Conj(UM1(gl1453,1))*UP1(gl1453,gO2);
   }
   tmp_1785 += tmp_1786;
   tmp_1780 += (std::complex<double>(0.,-0.7071067811865475)*LamSD*Conj(UP1(gI2
      ,1))*ZH(gI1,2)) * tmp_1785;
   std::complex<double> tmp_1787;
   std::complex<double> tmp_1788;
   for (unsigned gl1453 = 0; gl1453 < 2; ++gl1453) {
      tmp_1788 += Conj(UM1(gl1453,0))*UP1(gl1453,gO2);
   }
   tmp_1787 += tmp_1788;
   tmp_1780 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*ZH(gI1,3)) *
      tmp_1787;
   std::complex<double> tmp_1789;
   std::complex<double> tmp_1790;
   for (unsigned gl1453 = 0; gl1453 < 2; ++gl1453) {
      tmp_1790 += Conj(UM1(gl1453,1))*UP1(gl1453,gO2);
   }
   tmp_1789 += tmp_1790;
   tmp_1780 += (std::complex<double>(0,0.5)*LamTD*Conj(UP1(gI2,1))*ZH(gI1,3)) *
      tmp_1789;
   result += (std::complex<double>(0,-1)) * tmp_1780;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1hhCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1791;
   std::complex<double> tmp_1792;
   std::complex<double> tmp_1793;
   for (unsigned gl1456 = 0; gl1456 < 2; ++gl1456) {
      tmp_1793 += Conj(UM1(gl1456,gO1))*UP1(gl1456,1);
   }
   tmp_1792 += tmp_1793;
   tmp_1791 += (std::complex<double>(0.,-0.7071067811865475)*Conj(LamTD)*UM1(
      gI2,0)*ZH(gI1,0)) * tmp_1792;
   std::complex<double> tmp_1794;
   std::complex<double> tmp_1795;
   for (unsigned gl1456 = 0; gl1456 < 2; ++gl1456) {
      tmp_1795 += Conj(UM1(gl1456,gO1))*UP1(gl1456,0);
   }
   tmp_1794 += tmp_1795;
   tmp_1791 += (std::complex<double>(0.,-0.7071067811865475)*g2*UM1(gI2,1)*ZH(
      gI1,0)) * tmp_1794;
   std::complex<double> tmp_1796;
   std::complex<double> tmp_1797;
   for (unsigned gl1456 = 0; gl1456 < 2; ++gl1456) {
      tmp_1797 += Conj(UM1(gl1456,gO1))*UP1(gl1456,1);
   }
   tmp_1796 += tmp_1797;
   tmp_1791 += (std::complex<double>(0.,-0.7071067811865475)*Conj(LamSD)*UM1(
      gI2,1)*ZH(gI1,2)) * tmp_1796;
   std::complex<double> tmp_1798;
   std::complex<double> tmp_1799;
   for (unsigned gl1456 = 0; gl1456 < 2; ++gl1456) {
      tmp_1799 += Conj(UM1(gl1456,gO1))*UP1(gl1456,0);
   }
   tmp_1798 += tmp_1799;
   tmp_1791 += (std::complex<double>(0,-1)*g2*UM1(gI2,0)*ZH(gI1,3)) * tmp_1798;
   std::complex<double> tmp_1800;
   std::complex<double> tmp_1801;
   for (unsigned gl1456 = 0; gl1456 < 2; ++gl1456) {
      tmp_1801 += Conj(UM1(gl1456,gO1))*UP1(gl1456,1);
   }
   tmp_1800 += tmp_1801;
   tmp_1791 += (std::complex<double>(0,0.5)*Conj(LamTD)*UM1(gI2,1)*ZH(gI1,3)) *
      tmp_1800;
   result += (std::complex<double>(0,-1)) * tmp_1791;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1802;
   std::complex<double> tmp_1803;
   std::complex<double> tmp_1804;
   for (unsigned gl1459 = 0; gl1459 < 2; ++gl1459) {
      tmp_1804 += Conj(UM1(gl1459,1))*UP1(gl1459,gO2);
   }
   tmp_1803 += tmp_1804;
   tmp_1802 += (std::complex<double>(0.,0.5477225575051661)*g1*Conj(ZN1(gI2,0))
      *ZP(gI1,0)) * tmp_1803;
   std::complex<double> tmp_1805;
   std::complex<double> tmp_1806;
   for (unsigned gl1459 = 0; gl1459 < 2; ++gl1459) {
      tmp_1806 += Conj(UM1(gl1459,1))*UP1(gl1459,gO2);
   }
   tmp_1805 += tmp_1806;
   tmp_1802 += (std::complex<double>(0.,0.7071067811865475)*g2*Conj(ZN1(gI2,1))
      *ZP(gI1,0)) * tmp_1805;
   std::complex<double> tmp_1807;
   std::complex<double> tmp_1808;
   for (unsigned gl1459 = 0; gl1459 < 2; ++gl1459) {
      tmp_1808 += Conj(UM1(gl1459,0))*UP1(gl1459,gO2);
   }
   tmp_1807 += tmp_1808;
   tmp_1802 += (std::complex<double>(0,-1)*LamTU*Conj(ZN1(gI2,3))*ZP(gI1,1)) *
      tmp_1807;
   std::complex<double> tmp_1809;
   std::complex<double> tmp_1810;
   for (unsigned gl1459 = 0; gl1459 < 2; ++gl1459) {
      tmp_1810 += Conj(UM1(gl1459,0))*UP1(gl1459,gO2);
   }
   tmp_1809 += tmp_1810;
   tmp_1802 += (std::complex<double>(0.,1.4142135623730951)*g2*Conj(ZN1(gI2,1))
      *ZP(gI1,2)) * tmp_1809;
   std::complex<double> tmp_1811;
   std::complex<double> tmp_1812;
   for (unsigned gl1459 = 0; gl1459 < 2; ++gl1459) {
      tmp_1812 += Conj(UM1(gl1459,1))*UP1(gl1459,gO2);
   }
   tmp_1811 += tmp_1812;
   tmp_1802 += (std::complex<double>(0,1)*LamTD*Conj(ZN1(gI2,2))*ZP(gI1,3)) *
      tmp_1811;
   result += (std::complex<double>(0,-1)) * tmp_1802;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1813;
   std::complex<double> tmp_1814;
   std::complex<double> tmp_1815;
   for (unsigned gl1462 = 0; gl1462 < 2; ++gl1462) {
      tmp_1815 += Conj(UM1(gl1462,gO1))*UP1(gl1462,1);
   }
   tmp_1814 += tmp_1815;
   tmp_1813 += (std::complex<double>(0,-1)*Conj(LamSD)*ZN2(gI2,0)*ZP(gI1,0)) *
      tmp_1814;
   std::complex<double> tmp_1816;
   std::complex<double> tmp_1817;
   for (unsigned gl1462 = 0; gl1462 < 2; ++gl1462) {
      tmp_1817 += Conj(UM1(gl1462,gO1))*UP1(gl1462,1);
   }
   tmp_1816 += tmp_1817;
   tmp_1813 += (std::complex<double>(0.,0.7071067811865475)*Conj(LamTD)*ZN2(gI2
      ,1)*ZP(gI1,0)) * tmp_1816;
   std::complex<double> tmp_1818;
   std::complex<double> tmp_1819;
   for (unsigned gl1462 = 0; gl1462 < 2; ++gl1462) {
      tmp_1819 += Conj(UM1(gl1462,gO1))*UP1(gl1462,0);
   }
   tmp_1818 += tmp_1819;
   tmp_1813 += (std::complex<double>(0,-1)*g2*ZN2(gI2,3)*ZP(gI1,1)) * tmp_1818;
   std::complex<double> tmp_1820;
   std::complex<double> tmp_1821;
   for (unsigned gl1462 = 0; gl1462 < 2; ++gl1462) {
      tmp_1821 += Conj(UM1(gl1462,gO1))*UP1(gl1462,1);
   }
   tmp_1820 += tmp_1821;
   tmp_1813 += (std::complex<double>(0,-1)*Conj(LamTD)*ZN2(gI2,2)*ZP(gI1,2)) *
      tmp_1820;
   std::complex<double> tmp_1822;
   std::complex<double> tmp_1823;
   for (unsigned gl1462 = 0; gl1462 < 2; ++gl1462) {
      tmp_1823 += Conj(UM1(gl1462,gO1))*UP1(gl1462,0);
   }
   tmp_1822 += tmp_1823;
   tmp_1813 += (std::complex<double>(0.,1.4142135623730951)*g2*ZN2(gI2,1)*ZP(
      gI1,3)) * tmp_1822;
   result += (std::complex<double>(0,-1)) * tmp_1813;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barChiSRdpPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_1824;
   std::complex<double> tmp_1825;
   std::complex<double> tmp_1826;
   for (unsigned gl1465 = 0; gl1465 < 2; ++gl1465) {
      tmp_1826 += Conj(UM1(gl1465,0))*UP1(gl1465,gO2);
   }
   tmp_1825 += tmp_1826;
   tmp_1824 += (std::complex<double>(0,-1)*LamTD*Conj(ZN2(gI1,2))) * tmp_1825;
   std::complex<double> tmp_1827;
   std::complex<double> tmp_1828;
   for (unsigned gl1465 = 0; gl1465 < 2; ++gl1465) {
      tmp_1828 += Conj(UM1(gl1465,1))*UP1(gl1465,gO2);
   }
   tmp_1827 += tmp_1828;
   tmp_1824 += (std::complex<double>(0,-1)*LamSD*Conj(ZN2(gI1,0))) * tmp_1827;
   std::complex<double> tmp_1829;
   std::complex<double> tmp_1830;
   for (unsigned gl1465 = 0; gl1465 < 2; ++gl1465) {
      tmp_1830 += Conj(UM1(gl1465,1))*UP1(gl1465,gO2);
   }
   tmp_1829 += tmp_1830;
   tmp_1824 += (std::complex<double>(0.,0.7071067811865475)*LamTD*Conj(ZN2(gI1,
      1))) * tmp_1829;
   result += (std::complex<double>(0,-1)) * tmp_1824;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barChiSRdpPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_1831;
   std::complex<double> tmp_1832;
   std::complex<double> tmp_1833;
   for (unsigned gl1468 = 0; gl1468 < 2; ++gl1468) {
      tmp_1833 += Conj(UM1(gl1468,gO1))*UP1(gl1468,1);
   }
   tmp_1832 += tmp_1833;
   tmp_1831 += (std::complex<double>(0.,-0.5477225575051661)*g1*ZN1(gI1,0)) *
      tmp_1832;
   std::complex<double> tmp_1834;
   std::complex<double> tmp_1835;
   for (unsigned gl1468 = 0; gl1468 < 2; ++gl1468) {
      tmp_1835 += Conj(UM1(gl1468,gO1))*UP1(gl1468,1);
   }
   tmp_1834 += tmp_1835;
   tmp_1831 += (std::complex<double>(0.,-0.7071067811865475)*g2*ZN1(gI1,1)) *
      tmp_1834;
   std::complex<double> tmp_1836;
   std::complex<double> tmp_1837;
   for (unsigned gl1468 = 0; gl1468 < 2; ++gl1468) {
      tmp_1837 += Conj(UM1(gl1468,gO1))*UP1(gl1468,0);
   }
   tmp_1836 += tmp_1837;
   tmp_1831 += (std::complex<double>(0,-1)*g2*ZN1(gI1,2)) * tmp_1836;
   result += (std::complex<double>(0,-1)) * tmp_1831;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjSdFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1838;
   std::complex<double> tmp_1840;
   for (unsigned gl1471 = 0; gl1471 < 2; ++gl1471) {
      tmp_1840 += Conj(UM1(gl1471,1))*UP1(gl1471,gO2);
   }
   tmp_1838 += tmp_1840;
   std::complex<double> tmp_1839;
   std::complex<double> tmp_1841;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1842;
      std::complex<double> tmp_1843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1843 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1842 += tmp_1843;
      tmp_1841 += (Conj(ZUL(gI2,j2))) * tmp_1842;
   }
   tmp_1839 += tmp_1841;
   result += (1) * tmp_1838 * tmp_1839;

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

   std::complex<double> tmp_1844;
   std::complex<double> tmp_1846;
   for (unsigned gl1477 = 0; gl1477 < 2; ++gl1477) {
      tmp_1846 += Conj(UM1(gl1477,1))*UP1(gl1477,gO2);
   }
   tmp_1844 += tmp_1846;
   std::complex<double> tmp_1845;
   std::complex<double> tmp_1847;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1847 += Ye(j1,gI2)*ZE(gI1,3 + j1);
   }
   tmp_1845 += tmp_1847;
   result += (1) * tmp_1844 * tmp_1845;

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

   std::complex<double> tmp_1848;
   std::complex<double> tmp_1849;
   std::complex<double> tmp_1850;
   for (unsigned gl1483 = 0; gl1483 < 2; ++gl1483) {
      tmp_1850 += Conj(UM1(gl1483,0))*UP1(gl1483,gO2);
   }
   tmp_1849 += tmp_1850;
   tmp_1848 += (std::complex<double>(0,-1)*g2*Sin(ThetaW())*UM1(gI2,0)) *
      tmp_1849;
   std::complex<double> tmp_1851;
   std::complex<double> tmp_1852;
   for (unsigned gl1483 = 0; gl1483 < 2; ++gl1483) {
      tmp_1852 += Conj(UM1(gl1483,1))*UP1(gl1483,gO2);
   }
   tmp_1851 += tmp_1852;
   tmp_1848 += (std::complex<double>(0.,-0.3872983346207417)*g1*Cos(ThetaW())*
      UM1(gI2,1)) * tmp_1851;
   std::complex<double> tmp_1853;
   std::complex<double> tmp_1854;
   for (unsigned gl1483 = 0; gl1483 < 2; ++gl1483) {
      tmp_1854 += Conj(UM1(gl1483,1))*UP1(gl1483,gO2);
   }
   tmp_1853 += tmp_1854;
   tmp_1848 += (std::complex<double>(0,-0.5)*g2*Sin(ThetaW())*UM1(gI2,1)) *
      tmp_1853;
   result += (std::complex<double>(0,-1)) * tmp_1848;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VPCha1PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1855;
   std::complex<double> tmp_1856;
   std::complex<double> tmp_1857;
   for (unsigned gl1486 = 0; gl1486 < 2; ++gl1486) {
      tmp_1857 += Conj(UM1(gl1486,gO1))*UP1(gl1486,1);
   }
   tmp_1856 += tmp_1857;
   tmp_1855 += (std::complex<double>(0.,-0.3872983346207417)*g1*Conj(UP1(gI2,1)
      )*Cos(ThetaW())) * tmp_1856;
   std::complex<double> tmp_1858;
   std::complex<double> tmp_1859;
   for (unsigned gl1486 = 0; gl1486 < 2; ++gl1486) {
      tmp_1859 += Conj(UM1(gl1486,gO1))*UP1(gl1486,0);
   }
   tmp_1858 += tmp_1859;
   tmp_1855 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*Sin(ThetaW())) *
      tmp_1858;
   std::complex<double> tmp_1860;
   std::complex<double> tmp_1861;
   for (unsigned gl1486 = 0; gl1486 < 2; ++gl1486) {
      tmp_1861 += Conj(UM1(gl1486,gO1))*UP1(gl1486,1);
   }
   tmp_1860 += tmp_1861;
   tmp_1855 += (std::complex<double>(0,-0.5)*g2*Conj(UP1(gI2,1))*Sin(ThetaW()))
      * tmp_1860;
   result += (std::complex<double>(0,-1)) * tmp_1855;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VZCha1PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1862;
   std::complex<double> tmp_1863;
   std::complex<double> tmp_1864;
   for (unsigned gl1489 = 0; gl1489 < 2; ++gl1489) {
      tmp_1864 += Conj(UM1(gl1489,0))*UP1(gl1489,gO2);
   }
   tmp_1863 += tmp_1864;
   tmp_1862 += (std::complex<double>(0,-1)*g2*Cos(ThetaW())*UM1(gI2,0)) *
      tmp_1863;
   std::complex<double> tmp_1865;
   std::complex<double> tmp_1866;
   for (unsigned gl1489 = 0; gl1489 < 2; ++gl1489) {
      tmp_1866 += Conj(UM1(gl1489,1))*UP1(gl1489,gO2);
   }
   tmp_1865 += tmp_1866;
   tmp_1862 += (std::complex<double>(0,-0.5)*g2*Cos(ThetaW())*UM1(gI2,1)) *
      tmp_1865;
   std::complex<double> tmp_1867;
   std::complex<double> tmp_1868;
   for (unsigned gl1489 = 0; gl1489 < 2; ++gl1489) {
      tmp_1868 += Conj(UM1(gl1489,1))*UP1(gl1489,gO2);
   }
   tmp_1867 += tmp_1868;
   tmp_1862 += (std::complex<double>(0.,0.3872983346207417)*g1*Sin(ThetaW())*
      UM1(gI2,1)) * tmp_1867;
   result += (std::complex<double>(0,-1)) * tmp_1862;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VZCha1PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1869;
   std::complex<double> tmp_1870;
   std::complex<double> tmp_1871;
   for (unsigned gl1492 = 0; gl1492 < 2; ++gl1492) {
      tmp_1871 += Conj(UM1(gl1492,gO1))*UP1(gl1492,0);
   }
   tmp_1870 += tmp_1871;
   tmp_1869 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*Cos(ThetaW())) *
      tmp_1870;
   std::complex<double> tmp_1872;
   std::complex<double> tmp_1873;
   for (unsigned gl1492 = 0; gl1492 < 2; ++gl1492) {
      tmp_1873 += Conj(UM1(gl1492,gO1))*UP1(gl1492,1);
   }
   tmp_1872 += tmp_1873;
   tmp_1869 += (std::complex<double>(0,-0.5)*g2*Conj(UP1(gI2,1))*Cos(ThetaW()))
      * tmp_1872;
   std::complex<double> tmp_1874;
   std::complex<double> tmp_1875;
   for (unsigned gl1492 = 0; gl1492 < 2; ++gl1492) {
      tmp_1875 += Conj(UM1(gl1492,gO1))*UP1(gl1492,1);
   }
   tmp_1874 += tmp_1875;
   tmp_1869 += (std::complex<double>(0.,0.3872983346207417)*g1*Conj(UP1(gI2,1))
      *Sin(ThetaW())) * tmp_1874;
   result += (std::complex<double>(0,-1)) * tmp_1869;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1876;
   std::complex<double> tmp_1877;
   std::complex<double> tmp_1878;
   for (unsigned gl1495 = 0; gl1495 < 2; ++gl1495) {
      tmp_1878 += Conj(UM1(gl1495,0))*UP1(gl1495,gO2);
   }
   tmp_1877 += tmp_1878;
   tmp_1876 += (std::complex<double>(0,1)*g2*ZN2(gI2,1)) * tmp_1877;
   std::complex<double> tmp_1879;
   std::complex<double> tmp_1880;
   for (unsigned gl1495 = 0; gl1495 < 2; ++gl1495) {
      tmp_1880 += Conj(UM1(gl1495,1))*UP1(gl1495,gO2);
   }
   tmp_1879 += tmp_1880;
   tmp_1876 += (std::complex<double>(0.,0.7071067811865475)*g2*ZN2(gI2,2)) *
      tmp_1879;
   result += (std::complex<double>(0,-1)) * tmp_1876;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1881;
   std::complex<double> tmp_1882;
   std::complex<double> tmp_1883;
   for (unsigned gl1498 = 0; gl1498 < 2; ++gl1498) {
      tmp_1883 += Conj(UM1(gl1498,gO1))*UP1(gl1498,0);
   }
   tmp_1882 += tmp_1883;
   tmp_1881 += (std::complex<double>(0,1)*g2*Conj(ZN1(gI2,1))) * tmp_1882;
   std::complex<double> tmp_1884;
   std::complex<double> tmp_1885;
   for (unsigned gl1498 = 0; gl1498 < 2; ++gl1498) {
      tmp_1885 += Conj(UM1(gl1498,gO1))*UP1(gl1498,1);
   }
   tmp_1884 += tmp_1885;
   tmp_1881 += (std::complex<double>(0.,-0.7071067811865475)*g2*Conj(ZN1(gI2,2)
      )) * tmp_1884;
   result += (std::complex<double>(0,-1)) * tmp_1881;

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

   std::complex<double> tmp_1886;
   std::complex<double> tmp_1887;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1888;
      std::complex<double> tmp_1889;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1889 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1888 += tmp_1889;
      tmp_1887 += (Conj(ZD(gI2,j2))) * tmp_1888;
   }
   tmp_1886 += tmp_1887;
   result += (KroneckerDelta(1,gO2)) * tmp_1886;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1890;
   std::complex<double> tmp_1891;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1891 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1890 += tmp_1891;
   result += (-(g2*KroneckerDelta(0,gO1))) * tmp_1890;

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

std::complex<double> CLASSNAME::CpbarUCha2barChiSRumPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = LamTU*Conj(ZN2(gI1,3))*KroneckerDelta(0,gO2) + LamSU*Conj(ZN2(gI1,0
      ))*KroneckerDelta(1,gO2) + 0.7071067811865475*LamTU*Conj(ZN2(gI1,1))*
      KroneckerDelta(1,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barChiSRumPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = KroneckerDelta(1,gO1)*(0.5477225575051661*g1*ZN1(gI1,0) +
      0.7071067811865475*g2*ZN1(gI1,1)) - g2*KroneckerDelta(0,gO1)*ZN1(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2conjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1892;
   std::complex<double> tmp_1893;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1894;
      std::complex<double> tmp_1895;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1895 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1894 += tmp_1895;
      tmp_1893 += (Conj(ZDL(gI2,j2))) * tmp_1894;
   }
   tmp_1892 += tmp_1893;
   result += (KroneckerDelta(1,gO2)) * tmp_1892;

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
      std::complex<double> tmp_1896;
      std::complex<double> tmp_1897;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1897 += Conj(ZV(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1896 += tmp_1897;
      result += (Conj(UM1(gI1,1))) * tmp_1896;
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
      std::complex<double> tmp_1898;
      std::complex<double> tmp_1899;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1899 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1898 += tmp_1899;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1898;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1900;
      std::complex<double> tmp_1901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1901 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_1900 += tmp_1901;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1900;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1902;
      std::complex<double> tmp_1903;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1903 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1902 += tmp_1903;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1902;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1904;
      std::complex<double> tmp_1905;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1905 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_1904 += tmp_1905;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1904;
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
      std::complex<double> tmp_1906;
      std::complex<double> tmp_1907;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1907 += Conj(ZE(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1906 += tmp_1907;
      result += (-Conj(ZN2(gI1,2))) * tmp_1906;
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
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1908 += tmp_1909;
      result += (-ZN2(gI2,2)) * tmp_1908;
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
      std::complex<double> tmp_1910;
      std::complex<double> tmp_1911;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1911 += Conj(ZU(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1910 += tmp_1911;
      result += (Conj(UM1(gI1,1))) * tmp_1910;
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
      std::complex<double> tmp_1912;
      std::complex<double> tmp_1913;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1913 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1912 += tmp_1913;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1912;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1914;
      std::complex<double> tmp_1915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1915 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_1914 += tmp_1915;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1914;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1916;
      std::complex<double> tmp_1917;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1917 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1916 += tmp_1917;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1916;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1918;
      std::complex<double> tmp_1919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1919 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1918 += tmp_1919;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1918;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1920;
      std::complex<double> tmp_1921;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1921 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1920 += tmp_1921;
      result += (ZP(gI1,0)) * tmp_1920;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1922;
      std::complex<double> tmp_1923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1923 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1922 += tmp_1923;
      result += (ZP(gI1,1)) * tmp_1922;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1924;
      std::complex<double> tmp_1925;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1925 += Conj(ZD(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1924 += tmp_1925;
      result += (-Conj(ZN2(gI1,2))) * tmp_1924;
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
      std::complex<double> tmp_1926;
      std::complex<double> tmp_1927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1927 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1926 += tmp_1927;
      result += (UP2(gI2,1)) * tmp_1926;
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
      std::complex<double> tmp_1928;
      std::complex<double> tmp_1929;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1929 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1928 += tmp_1929;
      result += (-ZN2(gI2,2)) * tmp_1928;
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
      std::complex<double> tmp_1930;
      std::complex<double> tmp_1931;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1931 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1930 += tmp_1931;
      result += (Conj(UP2(gI1,1))) * tmp_1930;
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
      std::complex<double> tmp_1932;
      std::complex<double> tmp_1933;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1933 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1932 += tmp_1933;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_1932;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1934;
      std::complex<double> tmp_1935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1935 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_1934 += tmp_1935;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_1934;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1936;
      std::complex<double> tmp_1937;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1937 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1936 += tmp_1937;
      result += (ZP(gI1,1)) * tmp_1936;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1938;
      std::complex<double> tmp_1939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1939 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1938 += tmp_1939;
      result += (ZP(gI1,0)) * tmp_1938;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1940;
      std::complex<double> tmp_1941;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1941 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1940 += tmp_1941;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1940;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1942;
      std::complex<double> tmp_1943;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1943 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1942 += tmp_1943;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1942;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1944;
      std::complex<double> tmp_1945;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1945 += Conj(ZU(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1944 += tmp_1945;
      result += (-Conj(ZN2(gI1,3))) * tmp_1944;
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
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1946 += tmp_1947;
      result += (UM1(gI2,1)) * tmp_1946;
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
      std::complex<double> tmp_1948;
      std::complex<double> tmp_1949;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1949 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1948 += tmp_1949;
      result += (-ZN2(gI2,3)) * tmp_1948;
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

   std::complex<double> tmp_1950;
   std::complex<double> tmp_1951;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1951 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1950 += tmp_1951;
   result += (-1.4142135623730951*g3) * tmp_1950;

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

   std::complex<double> tmp_1952;
   std::complex<double> tmp_1953;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1953 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1952 += tmp_1953;
   result += (-1.4142135623730951*g3) * tmp_1952;

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

   std::complex<double> tmp_1954;
   std::complex<double> tmp_1955;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1955 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1954 += tmp_1955;
   result += (1.4142135623730951*g3) * tmp_1954;

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

   std::complex<double> tmp_1956;
   std::complex<double> tmp_1957;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1957 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1956 += tmp_1957;
   result += (1.4142135623730951*g3) * tmp_1956;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluphiOGluPL() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluphiOGluPR() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpbarGlusigmaOGluPL() const
{
   double result = 0.0;

   result = -g3;

   return result;
}

double CLASSNAME::CpbarGlusigmaOGluPR() const
{
   double result = 0.0;

   result = g3;

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

double CLASSNAME::CpbarFvbarCha2SePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvbarCha2SePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZE(gI2,gO1))*UM2(gI1,0));
   }

   return result;
}

double CLASSNAME::CpbarFvbarChiSvPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvbarChiSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZV(gI2,gO1))*ZN1(gI1,0);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZV(gI2,gO1))*ZN1(gI1,1);
   }

   return result;
}

double CLASSNAME::CpbarFvconjHpmFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvconjHpmFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1958;
   std::complex<double> tmp_1959;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1959 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
   }
   tmp_1958 += tmp_1959;
   result += (ZP(gI1,0)) * tmp_1958;

   return result;
}

double CLASSNAME::CpbarFvSeCha1PL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvSeCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1960;
   std::complex<double> tmp_1961;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1961 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
   }
   tmp_1960 += tmp_1961;
   result += (UM1(gI2,1)) * tmp_1960;

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
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarCha1SvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
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
      tmp_1963 += (Conj(ZV(gI2,j2))) * tmp_1964;
   }
   tmp_1962 += tmp_1963;
   result += (Conj(UM1(gI1,1))) * tmp_1962;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarCha1SvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1966;
   std::complex<double> tmp_1967;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1967 += Conj(ZV(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_1966 += tmp_1967;
   result += (-(g2*UP1(gI1,0))) * tmp_1966;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1968;
   std::complex<double> tmp_1969;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1970;
      std::complex<double> tmp_1971;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1971 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1970 += tmp_1971;
      tmp_1969 += (Conj(ZEL(gI1,j2))) * tmp_1970;
   }
   tmp_1968 += tmp_1969;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_1968;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1972;
   std::complex<double> tmp_1973;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1974;
      std::complex<double> tmp_1975;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1975 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_1974 += tmp_1975;
      tmp_1973 += (ZEL(gO1,j2)) * tmp_1974;
   }
   tmp_1972 += tmp_1973;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_1972
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1976;
   std::complex<double> tmp_1977;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1978;
      std::complex<double> tmp_1979;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1979 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1978 += tmp_1979;
      tmp_1977 += (Conj(ZEL(gI2,j2))) * tmp_1978;
   }
   tmp_1976 += tmp_1977;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1976;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1980;
   std::complex<double> tmp_1981;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1982;
      std::complex<double> tmp_1983;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1983 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1982 += tmp_1983;
      tmp_1981 += (ZEL(gO1,j2)) * tmp_1982;
   }
   tmp_1980 += tmp_1981;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1980;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1984;
   std::complex<double> tmp_1985;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1985 += Conj(ZER(gO2,j1))*Ye(j1,gI2);
   }
   tmp_1984 += tmp_1985;
   result += (ZP(gI1,0)) * tmp_1984;

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

   std::complex<double> tmp_1986;
   std::complex<double> tmp_1987;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1988;
      std::complex<double> tmp_1989;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1989 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1988 += tmp_1989;
      tmp_1987 += (Conj(ZE(gI2,j2))) * tmp_1988;
   }
   tmp_1986 += tmp_1987;
   result += (-Conj(ZN2(gI1,2))) * tmp_1986;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarChiSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1990;
   std::complex<double> tmp_1991;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1991 += Conj(ZE(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_1990 += tmp_1991;
   result += (0.7071067811865475*(0.7745966692414834*g1*ZN1(gI1,0) + g2*ZN1(gI1
      ,1))) * tmp_1990;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1992;
   std::complex<double> tmp_1993;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1993 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_1992 += tmp_1993;
   result += (-1.0954451150103321*g1*Conj(ZN1(gI2,0))) * tmp_1992;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1994;
   std::complex<double> tmp_1995;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1996;
      std::complex<double> tmp_1997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1997 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1996 += tmp_1997;
      tmp_1995 += (ZEL(gO1,j2)) * tmp_1996;
   }
   tmp_1994 += tmp_1995;
   result += (-ZN2(gI2,2)) * tmp_1994;

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

   std::complex<double> tmp_1998;
   std::complex<double> tmp_1999;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2000;
      std::complex<double> tmp_2001;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2001 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2000 += tmp_2001;
      tmp_1999 += (Conj(ZU(gI2,j2))) * tmp_2000;
   }
   tmp_1998 += tmp_1999;
   result += (Conj(UM1(gI1,1))) * tmp_1998;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarCha1SuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2002;
   std::complex<double> tmp_2003;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2003 += Conj(ZU(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2002 += tmp_2003;
   result += (-(g2*UP1(gI1,0))) * tmp_2002;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
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
      tmp_2005 += (Conj(ZDL(gI1,j2))) * tmp_2006;
   }
   tmp_2004 += tmp_2005;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_2004;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2008;
   std::complex<double> tmp_2009;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2010;
      std::complex<double> tmp_2011;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2011 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_2010 += tmp_2011;
      tmp_2009 += (ZDL(gO1,j2)) * tmp_2010;
   }
   tmp_2008 += tmp_2009;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_2008
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2012;
   std::complex<double> tmp_2013;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2014;
      std::complex<double> tmp_2015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2015 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2014 += tmp_2015;
      tmp_2013 += (Conj(ZDL(gI2,j2))) * tmp_2014;
   }
   tmp_2012 += tmp_2013;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2012;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2016;
   std::complex<double> tmp_2017;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2018;
      std::complex<double> tmp_2019;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2019 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2018 += tmp_2019;
      tmp_2017 += (ZDL(gO1,j2)) * tmp_2018;
   }
   tmp_2016 += tmp_2017;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2016;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2020;
   std::complex<double> tmp_2021;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2022;
      std::complex<double> tmp_2023;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2023 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2022 += tmp_2023;
      tmp_2021 += (Conj(ZUL(gI2,j2))) * tmp_2022;
   }
   tmp_2020 += tmp_2021;
   result += (ZP(gI1,0)) * tmp_2020;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2024;
   std::complex<double> tmp_2025;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2026;
      std::complex<double> tmp_2027;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2027 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2026 += tmp_2027;
      tmp_2025 += (ZDL(gO1,j2)) * tmp_2026;
   }
   tmp_2024 += tmp_2025;
   result += (ZP(gI1,1)) * tmp_2024;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2028;
   std::complex<double> tmp_2029;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2030;
      std::complex<double> tmp_2031;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2031 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2030 += tmp_2031;
      tmp_2029 += (Conj(ZD(gI2,j2))) * tmp_2030;
   }
   tmp_2028 += tmp_2029;
   result += (-Conj(ZN2(gI1,2))) * tmp_2028;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarChiSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2032;
   std::complex<double> tmp_2033;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2033 += Conj(ZD(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2032 += tmp_2033;
   result += (-0.2357022603955158*(0.7745966692414834*g1*ZN1(gI1,0) - 3*g2*ZN1(
      gI1,1))) * tmp_2032;

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

   std::complex<double> tmp_2034;
   std::complex<double> tmp_2035;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2036;
      std::complex<double> tmp_2037;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2037 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2036 += tmp_2037;
      tmp_2035 += (ZDL(gO1,j2)) * tmp_2036;
   }
   tmp_2034 += tmp_2035;
   result += (UP2(gI2,1)) * tmp_2034;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2038;
   std::complex<double> tmp_2039;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2039 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2038 += tmp_2039;
   result += (-0.3651483716701107*g1*Conj(ZN1(gI2,0))) * tmp_2038;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2040;
   std::complex<double> tmp_2041;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2042;
      std::complex<double> tmp_2043;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2043 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2042 += tmp_2043;
      tmp_2041 += (ZDL(gO1,j2)) * tmp_2042;
   }
   tmp_2040 += tmp_2041;
   result += (-ZN2(gI2,2)) * tmp_2040;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2044;
   std::complex<double> tmp_2045;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2045 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2044 += tmp_2045;
   result += (1.4142135623730951*g3) * tmp_2044;

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

   std::complex<double> tmp_2046;
   std::complex<double> tmp_2047;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2047 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2046 += tmp_2047;
   result += (-0.7071067811865475*g2) * tmp_2046;

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

   std::complex<double> tmp_2048;
   std::complex<double> tmp_2049;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2049 += Conj(ZD(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2048 += tmp_2049;
   result += (-1.4142135623730951*g3) * tmp_2048;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarCha2SdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2050;
   std::complex<double> tmp_2051;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2052;
      std::complex<double> tmp_2053;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2053 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2052 += tmp_2053;
      tmp_2051 += (Conj(ZD(gI2,j2))) * tmp_2052;
   }
   tmp_2050 += tmp_2051;
   result += (Conj(UP2(gI1,1))) * tmp_2050;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarCha2SdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2055 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2054 += tmp_2055;
   result += (-(g2*UM2(gI1,0))) * tmp_2054;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
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
      tmp_2057 += (Conj(ZUL(gI1,j2))) * tmp_2058;
   }
   tmp_2056 += tmp_2057;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
      tmp_2056;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2060;
   std::complex<double> tmp_2061;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2062;
      std::complex<double> tmp_2063;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2063 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2062 += tmp_2063;
      tmp_2061 += (ZUL(gO1,j2)) * tmp_2062;
   }
   tmp_2060 += tmp_2061;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_2060
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2064;
   std::complex<double> tmp_2065;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2066;
      std::complex<double> tmp_2067;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2067 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2066 += tmp_2067;
      tmp_2065 += (Conj(ZDL(gI2,j2))) * tmp_2066;
   }
   tmp_2064 += tmp_2065;
   result += (ZP(gI1,1)) * tmp_2064;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2068;
   std::complex<double> tmp_2069;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2070;
      std::complex<double> tmp_2071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2071 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2070 += tmp_2071;
      tmp_2069 += (ZUL(gO1,j2)) * tmp_2070;
   }
   tmp_2068 += tmp_2069;
   result += (ZP(gI1,0)) * tmp_2068;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2074;
      std::complex<double> tmp_2075;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2075 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2074 += tmp_2075;
      tmp_2073 += (Conj(ZUL(gI2,j2))) * tmp_2074;
   }
   tmp_2072 += tmp_2073;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2072;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2076;
   std::complex<double> tmp_2077;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2078;
      std::complex<double> tmp_2079;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2079 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2078 += tmp_2079;
      tmp_2077 += (ZUL(gO1,j2)) * tmp_2078;
   }
   tmp_2076 += tmp_2077;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2076;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2080;
   std::complex<double> tmp_2081;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2082;
      std::complex<double> tmp_2083;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2083 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2082 += tmp_2083;
      tmp_2081 += (Conj(ZU(gI2,j2))) * tmp_2082;
   }
   tmp_2080 += tmp_2081;
   result += (-Conj(ZN2(gI1,3))) * tmp_2080;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChiSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2084;
   std::complex<double> tmp_2085;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2085 += Conj(ZU(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2084 += tmp_2085;
   result += (-0.2357022603955158*(0.7745966692414834*g1*ZN1(gI1,0) + 3*g2*ZN1(
      gI1,1))) * tmp_2084;

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

   std::complex<double> tmp_2086;
   std::complex<double> tmp_2087;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2088;
      std::complex<double> tmp_2089;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2089 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2088 += tmp_2089;
      tmp_2087 += (ZUL(gO1,j2)) * tmp_2088;
   }
   tmp_2086 += tmp_2087;
   result += (UM1(gI2,1)) * tmp_2086;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2090;
   std::complex<double> tmp_2091;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2091 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2090 += tmp_2091;
   result += (0.7302967433402214*g1*Conj(ZN1(gI2,0))) * tmp_2090;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2092;
   std::complex<double> tmp_2093;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2094;
      std::complex<double> tmp_2095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2095 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2094 += tmp_2095;
      tmp_2093 += (ZUL(gO1,j2)) * tmp_2094;
   }
   tmp_2092 += tmp_2093;
   result += (-ZN2(gI2,3)) * tmp_2092;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2096;
   std::complex<double> tmp_2097;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2097 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2096 += tmp_2097;
   result += (1.4142135623730951*g3) * tmp_2096;

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

   std::complex<double> tmp_2098;
   std::complex<double> tmp_2099;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2099 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2098 += tmp_2099;
   result += (-0.7071067811865475*g2) * tmp_2098;

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

   std::complex<double> tmp_2100;
   std::complex<double> tmp_2101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2101 += Conj(ZU(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2100 += tmp_2101;
   result += (-1.4142135623730951*g3) * tmp_2100;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpUSdconjUSdconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUSdconjUSdconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2102;
   std::complex<double> tmp_2103;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2103 += A0(MRh(gI1))*CpUSdconjUSdconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2102 += tmp_2103;
   result += (-1) * tmp_2102;
   std::complex<double> tmp_2104;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2105;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2105 += (Conj(CpconjUSdbarCha1FuPL(gO2,gI1,gI2))*
            CpconjUSdbarCha1FuPL(gO1,gI1,gI2) + Conj(CpconjUSdbarCha1FuPR(gO2,gI1,
            gI2))*CpconjUSdbarCha1FuPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MFu(gI2));
      }
      tmp_2104 += tmp_2105;
   }
   result += tmp_2104;
   std::complex<double> tmp_2106;
   std::complex<double> tmp_2107;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2108;
      std::complex<double> tmp_2109;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2109 += B0(p,MCha1(gI1),MFu(gI2))*(Conj(CpconjUSdbarCha1FuPR
            (gO2,gI1,gI2))*CpconjUSdbarCha1FuPL(gO1,gI1,gI2) + Conj(
            CpconjUSdbarCha1FuPL(gO2,gI1,gI2))*CpconjUSdbarCha1FuPR(gO1,gI1,gI2))*
            MFu(gI2);
      }
      tmp_2108 += tmp_2109;
      tmp_2107 += (MCha1(gI1)) * tmp_2108;
   }
   tmp_2106 += tmp_2107;
   result += (-2) * tmp_2106;
   std::complex<double> tmp_2110;
   std::complex<double> tmp_2111;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2111 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2110 += tmp_2111;
   result += (-1) * tmp_2110;
   std::complex<double> tmp_2112;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2113;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2113 += (Conj(CpconjUSdFuCha2PL(gO2,gI1,gI2))*
            CpconjUSdFuCha2PL(gO1,gI1,gI2) + Conj(CpconjUSdFuCha2PR(gO2,gI1,gI2))*
            CpconjUSdFuCha2PR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha2(gI2));
      }
      tmp_2112 += tmp_2113;
   }
   result += tmp_2112;
   std::complex<double> tmp_2114;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2115;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2115 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2114 += tmp_2115;
   }
   result += tmp_2114;
   std::complex<double> tmp_2116;
   std::complex<double> tmp_2117;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2118;
      std::complex<double> tmp_2119;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2119 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2118 += tmp_2119;
      tmp_2117 += (MFd(gI1)) * tmp_2118;
   }
   tmp_2116 += tmp_2117;
   result += (-2) * tmp_2116;
   std::complex<double> tmp_2120;
   std::complex<double> tmp_2121;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2122;
      std::complex<double> tmp_2123;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2123 += B0(p,MFu(gI1),MCha2(gI2))*(Conj(CpconjUSdFuCha2PR(
            gO2,gI1,gI2))*CpconjUSdFuCha2PL(gO1,gI1,gI2) + Conj(CpconjUSdFuCha2PL(
            gO2,gI1,gI2))*CpconjUSdFuCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2122 += tmp_2123;
      tmp_2121 += (MFu(gI1)) * tmp_2122;
   }
   tmp_2120 += tmp_2121;
   result += (-2) * tmp_2120;
   std::complex<double> tmp_2124;
   std::complex<double> tmp_2125;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2125 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2124 += tmp_2125;
   result += (-0.5) * tmp_2124;
   std::complex<double> tmp_2126;
   std::complex<double> tmp_2127;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2127 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2126 += tmp_2127;
   result += (-1) * tmp_2126;
   std::complex<double> tmp_2128;
   std::complex<double> tmp_2129;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2129 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2128 += tmp_2129;
   result += (-0.5) * tmp_2128;
   std::complex<double> tmp_2130;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2131;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2131 += (Conj(CpconjUSdbarChiFdPL(gO2,gI1,gI2))*
            CpconjUSdbarChiFdPL(gO1,gI1,gI2) + Conj(CpconjUSdbarChiFdPR(gO2,gI1,
            gI2))*CpconjUSdbarChiFdPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFd(gI2));
      }
      tmp_2130 += tmp_2131;
   }
   result += tmp_2130;
   std::complex<double> tmp_2132;
   std::complex<double> tmp_2133;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2134;
      std::complex<double> tmp_2135;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2135 += B0(p,MChi(gI1),MFd(gI2))*(Conj(CpconjUSdbarChiFdPR(
            gO2,gI1,gI2))*CpconjUSdbarChiFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSdbarChiFdPL(gO2,gI1,gI2))*CpconjUSdbarChiFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2134 += tmp_2135;
      tmp_2133 += (MChi(gI1)) * tmp_2134;
   }
   tmp_2132 += tmp_2133;
   result += (-2) * tmp_2132;
   std::complex<double> tmp_2136;
   std::complex<double> tmp_2137;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2137 += B0(p,MSd(gI1),MphiO)*Conj(CpconjUSdSdphiO(gO2,gI1))*
         CpconjUSdSdphiO(gO1,gI1);
   }
   tmp_2136 += tmp_2137;
   result += (1.3333333333333333) * tmp_2136;
   std::complex<double> tmp_2138;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2138 += B0(p,MSu(gI1),MSRum)*Conj(CpconjUSdSuSRum(gO2,gI1))*
         CpconjUSdSuSRum(gO1,gI1);
   }
   result += tmp_2138;
   std::complex<double> tmp_2139;
   std::complex<double> tmp_2140;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2140 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2139 += tmp_2140;
   result += (-1) * tmp_2139;
   std::complex<double> tmp_2141;
   std::complex<double> tmp_2142;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2142 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2141 += tmp_2142;
   result += (-1) * tmp_2141;
   std::complex<double> tmp_2143;
   std::complex<double> tmp_2144;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2144 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2143 += tmp_2144;
   result += (-1) * tmp_2143;
   std::complex<double> tmp_2145;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2146;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2146 += B0(p,MSd(gI1),MRh(gI2))*Conj(CpconjUSdSdRh(gO2,gI1,
            gI2))*CpconjUSdSdRh(gO1,gI1,gI2);
      }
      tmp_2145 += tmp_2146;
   }
   result += tmp_2145;
   std::complex<double> tmp_2147;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2148;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2148 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2147 += tmp_2148;
   }
   result += tmp_2147;
   std::complex<double> tmp_2149;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2150;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2150 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2149 += tmp_2150;
   }
   result += tmp_2149;
   std::complex<double> tmp_2151;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2152;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2152 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2151 += tmp_2152;
   }
   result += tmp_2151;
   std::complex<double> tmp_2153;
   std::complex<double> tmp_2154;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2154 += (Conj(CpconjUSdbarGluFdPL(gO2,0,gI2))*CpconjUSdbarGluFdPL(
         gO1,0,gI2) + Conj(CpconjUSdbarGluFdPR(gO2,0,gI2))*CpconjUSdbarGluFdPR(gO1
         ,0,gI2))*G0(p,MGlu,MFd(gI2));
   }
   tmp_2153 += tmp_2154;
   result += (1.3333333333333333) * tmp_2153;
   std::complex<double> tmp_2155;
   std::complex<double> tmp_2156;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2156 += (Conj(CpconjUSdGluFdPL(gO2,0,gI2))*CpconjUSdGluFdPL(gO1,0,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,0,gI2))*CpconjUSdGluFdPR(gO1,0,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2155 += tmp_2156;
   result += (1.3333333333333333) * tmp_2155;
   std::complex<double> tmp_2157;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2157 += B0(p,MSRdp,MSu(gI2))*Conj(CpconjUSdconjSRdpSu(gO2,gI2))*
         CpconjUSdconjSRdpSu(gO1,gI2);
   }
   result += tmp_2157;
   std::complex<double> tmp_2158;
   std::complex<double> tmp_2159;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2159 += B0(p,MsigmaO,MSd(gI2))*Conj(CpconjUSdsigmaOSd(gO2,gI2))*
         CpconjUSdsigmaOSd(gO1,gI2);
   }
   tmp_2158 += tmp_2159;
   result += (1.3333333333333333) * tmp_2158;
   std::complex<double> tmp_2160;
   std::complex<double> tmp_2161;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2161 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2160 += tmp_2161;
   result += (1.3333333333333333) * tmp_2160;
   std::complex<double> tmp_2162;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2162 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2162;
   std::complex<double> tmp_2163;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2163 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2163;
   std::complex<double> tmp_2164;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2164 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2164;
   std::complex<double> tmp_2165;
   std::complex<double> tmp_2166;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2166 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdbarGluFdPR(gO2,0,gI2))*
         CpconjUSdbarGluFdPL(gO1,0,gI2) + Conj(CpconjUSdbarGluFdPL(gO2,0,gI2))*
         CpconjUSdbarGluFdPR(gO1,0,gI2))*MFd(gI2);
   }
   tmp_2165 += tmp_2166;
   result += (-2.6666666666666665*MGlu) * tmp_2165;
   std::complex<double> tmp_2167;
   std::complex<double> tmp_2168;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2168 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,0,gI2))*
         CpconjUSdGluFdPL(gO1,0,gI2) + Conj(CpconjUSdGluFdPL(gO2,0,gI2))*
         CpconjUSdGluFdPR(gO1,0,gI2))*MFd(gI2);
   }
   tmp_2167 += tmp_2168;
   result += (-2.6666666666666665*MGlu) * tmp_2167;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpUSvconjUSvconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUSvconjUSvconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2169;
   std::complex<double> tmp_2170;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2170 += A0(MRh(gI1))*CpUSvconjUSvconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2169 += tmp_2170;
   result += (-1) * tmp_2169;
   std::complex<double> tmp_2171;
   std::complex<double> tmp_2172;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2172 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2171 += tmp_2172;
   result += (-1) * tmp_2171;
   std::complex<double> tmp_2173;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2174;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2174 += (Conj(CpconjUSvFeCha1PL(gO2,gI1,gI2))*
            CpconjUSvFeCha1PL(gO1,gI1,gI2) + Conj(CpconjUSvFeCha1PR(gO2,gI1,gI2))*
            CpconjUSvFeCha1PR(gO1,gI1,gI2))*G0(p,MFe(gI1),MCha1(gI2));
      }
      tmp_2173 += tmp_2174;
   }
   result += tmp_2173;
   std::complex<double> tmp_2175;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2176;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2176 += B0(p,MSv(gI1),MAh(gI2))*Conj(CpconjUSvSvAh(gO2,gI1,
            gI2))*CpconjUSvSvAh(gO1,gI1,gI2);
      }
      tmp_2175 += tmp_2176;
   }
   result += tmp_2175;
   std::complex<double> tmp_2177;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2178;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2178 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2177 += tmp_2178;
   }
   result += tmp_2177;
   std::complex<double> tmp_2179;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2180;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2180 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2179 += tmp_2180;
   }
   result += tmp_2179;
   std::complex<double> tmp_2181;
   std::complex<double> tmp_2182;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2183;
      std::complex<double> tmp_2184;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2184 += B0(p,MFe(gI1),MCha1(gI2))*(Conj(CpconjUSvFeCha1PR(
            gO2,gI1,gI2))*CpconjUSvFeCha1PL(gO1,gI1,gI2) + Conj(CpconjUSvFeCha1PL(
            gO2,gI1,gI2))*CpconjUSvFeCha1PR(gO1,gI1,gI2))*MCha1(gI2);
      }
      tmp_2183 += tmp_2184;
      tmp_2182 += (MFe(gI1)) * tmp_2183;
   }
   tmp_2181 += tmp_2182;
   result += (-2) * tmp_2181;
   std::complex<double> tmp_2185;
   std::complex<double> tmp_2186;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2187;
      std::complex<double> tmp_2188;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2188 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2187 += tmp_2188;
      tmp_2186 += (MFv(gI1)) * tmp_2187;
   }
   tmp_2185 += tmp_2186;
   result += (-2) * tmp_2185;
   std::complex<double> tmp_2189;
   std::complex<double> tmp_2190;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2190 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2189 += tmp_2190;
   result += (-0.5) * tmp_2189;
   std::complex<double> tmp_2191;
   std::complex<double> tmp_2192;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2192 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2191 += tmp_2192;
   result += (-1) * tmp_2191;
   std::complex<double> tmp_2193;
   std::complex<double> tmp_2194;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2194 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2193 += tmp_2194;
   result += (-0.5) * tmp_2193;
   std::complex<double> tmp_2195;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2196;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2196 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2195 += tmp_2196;
   }
   result += tmp_2195;
   std::complex<double> tmp_2197;
   std::complex<double> tmp_2198;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2198 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2197 += tmp_2198;
   result += (-3) * tmp_2197;
   std::complex<double> tmp_2199;
   std::complex<double> tmp_2200;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2200 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2199 += tmp_2200;
   result += (-1) * tmp_2199;
   std::complex<double> tmp_2201;
   std::complex<double> tmp_2202;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2202 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2201 += tmp_2202;
   result += (-3) * tmp_2201;
   std::complex<double> tmp_2203;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2203 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2203;
   std::complex<double> tmp_2204;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2204 += B0(p,MSRdp,MSe(gI2))*Conj(CpconjUSvSRdpSe(gO2,gI2))*
         CpconjUSvSRdpSe(gO1,gI2);
   }
   result += tmp_2204;
   std::complex<double> tmp_2205;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2205 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2205;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpUSuconjUSuconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUSuconjUSuconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2206;
   std::complex<double> tmp_2207;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2207 += A0(MRh(gI1))*CpUSuconjUSuconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2206 += tmp_2207;
   result += (-1) * tmp_2206;
   std::complex<double> tmp_2208;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2209;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2209 += (Conj(CpconjUSubarCha2FdPL(gO2,gI1,gI2))*
            CpconjUSubarCha2FdPL(gO1,gI1,gI2) + Conj(CpconjUSubarCha2FdPR(gO2,gI1,
            gI2))*CpconjUSubarCha2FdPR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MFd(gI2));
      }
      tmp_2208 += tmp_2209;
   }
   result += tmp_2208;
   std::complex<double> tmp_2210;
   std::complex<double> tmp_2211;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2212;
      std::complex<double> tmp_2213;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2213 += B0(p,MCha2(gI1),MFd(gI2))*(Conj(CpconjUSubarCha2FdPR
            (gO2,gI1,gI2))*CpconjUSubarCha2FdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarCha2FdPL(gO2,gI1,gI2))*CpconjUSubarCha2FdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2212 += tmp_2213;
      tmp_2211 += (MCha2(gI1)) * tmp_2212;
   }
   tmp_2210 += tmp_2211;
   result += (-2) * tmp_2210;
   std::complex<double> tmp_2214;
   std::complex<double> tmp_2215;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2215 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2214 += tmp_2215;
   result += (-1) * tmp_2214;
   std::complex<double> tmp_2216;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2217;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2217 += (Conj(CpconjUSuFdCha1PL(gO2,gI1,gI2))*
            CpconjUSuFdCha1PL(gO1,gI1,gI2) + Conj(CpconjUSuFdCha1PR(gO2,gI1,gI2))*
            CpconjUSuFdCha1PR(gO1,gI1,gI2))*G0(p,MFd(gI1),MCha1(gI2));
      }
      tmp_2216 += tmp_2217;
   }
   result += tmp_2216;
   std::complex<double> tmp_2218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2219;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2219 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2218 += tmp_2219;
   }
   result += tmp_2218;
   std::complex<double> tmp_2220;
   std::complex<double> tmp_2221;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2222;
      std::complex<double> tmp_2223;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2223 += B0(p,MFd(gI1),MCha1(gI2))*(Conj(CpconjUSuFdCha1PR(
            gO2,gI1,gI2))*CpconjUSuFdCha1PL(gO1,gI1,gI2) + Conj(CpconjUSuFdCha1PL(
            gO2,gI1,gI2))*CpconjUSuFdCha1PR(gO1,gI1,gI2))*MCha1(gI2);
      }
      tmp_2222 += tmp_2223;
      tmp_2221 += (MFd(gI1)) * tmp_2222;
   }
   tmp_2220 += tmp_2221;
   result += (-2) * tmp_2220;
   std::complex<double> tmp_2224;
   std::complex<double> tmp_2225;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2226;
      std::complex<double> tmp_2227;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2227 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2226 += tmp_2227;
      tmp_2225 += (MFu(gI1)) * tmp_2226;
   }
   tmp_2224 += tmp_2225;
   result += (-2) * tmp_2224;
   std::complex<double> tmp_2228;
   std::complex<double> tmp_2229;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2229 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2228 += tmp_2229;
   result += (-0.5) * tmp_2228;
   std::complex<double> tmp_2230;
   std::complex<double> tmp_2231;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2231 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2230 += tmp_2231;
   result += (-1) * tmp_2230;
   std::complex<double> tmp_2232;
   std::complex<double> tmp_2233;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2233 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2232 += tmp_2233;
   result += (-0.5) * tmp_2232;
   std::complex<double> tmp_2234;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2235;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2235 += (Conj(CpconjUSubarChiFuPL(gO2,gI1,gI2))*
            CpconjUSubarChiFuPL(gO1,gI1,gI2) + Conj(CpconjUSubarChiFuPR(gO2,gI1,
            gI2))*CpconjUSubarChiFuPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFu(gI2));
      }
      tmp_2234 += tmp_2235;
   }
   result += tmp_2234;
   std::complex<double> tmp_2236;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2237;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2237 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2236 += tmp_2237;
   }
   result += tmp_2236;
   std::complex<double> tmp_2238;
   std::complex<double> tmp_2239;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2240;
      std::complex<double> tmp_2241;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2241 += B0(p,MChi(gI1),MFu(gI2))*(Conj(CpconjUSubarChiFuPR(
            gO2,gI1,gI2))*CpconjUSubarChiFuPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChiFuPL(gO2,gI1,gI2))*CpconjUSubarChiFuPR(gO1,gI1,gI2))*
            MFu(gI2);
      }
      tmp_2240 += tmp_2241;
      tmp_2239 += (MChi(gI1)) * tmp_2240;
   }
   tmp_2238 += tmp_2239;
   result += (-2) * tmp_2238;
   std::complex<double> tmp_2242;
   std::complex<double> tmp_2243;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2243 += B0(p,MSu(gI1),MphiO)*Conj(CpconjUSuSuphiO(gO2,gI1))*
         CpconjUSuSuphiO(gO1,gI1);
   }
   tmp_2242 += tmp_2243;
   result += (1.3333333333333333) * tmp_2242;
   std::complex<double> tmp_2244;
   std::complex<double> tmp_2245;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2245 += B0(p,MSu(gI1),MsigmaO)*Conj(CpconjUSuSusigmaO(gO2,gI1))*
         CpconjUSuSusigmaO(gO1,gI1);
   }
   tmp_2244 += tmp_2245;
   result += (1.3333333333333333) * tmp_2244;
   std::complex<double> tmp_2246;
   std::complex<double> tmp_2247;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2247 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2246 += tmp_2247;
   result += (-1) * tmp_2246;
   std::complex<double> tmp_2248;
   std::complex<double> tmp_2249;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2249 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2248 += tmp_2249;
   result += (-1) * tmp_2248;
   std::complex<double> tmp_2250;
   std::complex<double> tmp_2251;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2251 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2250 += tmp_2251;
   result += (-1) * tmp_2250;
   std::complex<double> tmp_2252;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2253;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2253 += B0(p,MSu(gI1),MRh(gI2))*Conj(CpconjUSuSuRh(gO2,gI1,
            gI2))*CpconjUSuSuRh(gO1,gI1,gI2);
      }
      tmp_2252 += tmp_2253;
   }
   result += tmp_2252;
   std::complex<double> tmp_2254;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2255;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2255 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2254 += tmp_2255;
   }
   result += tmp_2254;
   std::complex<double> tmp_2256;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2257;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2257 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2256 += tmp_2257;
   }
   result += tmp_2256;
   std::complex<double> tmp_2258;
   std::complex<double> tmp_2259;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2259 += (Conj(CpconjUSubarGluFuPL(gO2,0,gI2))*CpconjUSubarGluFuPL(
         gO1,0,gI2) + Conj(CpconjUSubarGluFuPR(gO2,0,gI2))*CpconjUSubarGluFuPR(gO1
         ,0,gI2))*G0(p,MGlu,MFu(gI2));
   }
   tmp_2258 += tmp_2259;
   result += (1.3333333333333333) * tmp_2258;
   std::complex<double> tmp_2260;
   std::complex<double> tmp_2261;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2261 += (Conj(CpconjUSuGluFuPL(gO2,0,gI2))*CpconjUSuGluFuPL(gO1,0,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,0,gI2))*CpconjUSuGluFuPR(gO1,0,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2260 += tmp_2261;
   result += (1.3333333333333333) * tmp_2260;
   std::complex<double> tmp_2262;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2262 += B0(p,MSRum,MSd(gI2))*Conj(CpconjUSuconjSRumSd(gO2,gI2))*
         CpconjUSuconjSRumSd(gO1,gI2);
   }
   result += tmp_2262;
   std::complex<double> tmp_2263;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2263 += B0(p,MSRdp,MSd(gI2))*Conj(CpconjUSuSRdpSd(gO2,gI2))*
         CpconjUSuSRdpSd(gO1,gI2);
   }
   result += tmp_2263;
   std::complex<double> tmp_2264;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2264 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2264;
   std::complex<double> tmp_2265;
   std::complex<double> tmp_2266;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2266 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2265 += tmp_2266;
   result += (1.3333333333333333) * tmp_2265;
   std::complex<double> tmp_2267;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2267 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2267;
   std::complex<double> tmp_2268;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2268 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2268;
   std::complex<double> tmp_2269;
   std::complex<double> tmp_2270;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2270 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSubarGluFuPR(gO2,0,gI2))*
         CpconjUSubarGluFuPL(gO1,0,gI2) + Conj(CpconjUSubarGluFuPL(gO2,0,gI2))*
         CpconjUSubarGluFuPR(gO1,0,gI2))*MFu(gI2);
   }
   tmp_2269 += tmp_2270;
   result += (-2.6666666666666665*MGlu) * tmp_2269;
   std::complex<double> tmp_2271;
   std::complex<double> tmp_2272;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2272 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,0,gI2))*
         CpconjUSuGluFuPL(gO1,0,gI2) + Conj(CpconjUSuGluFuPL(gO2,0,gI2))*
         CpconjUSuGluFuPR(gO1,0,gI2))*MFu(gI2);
   }
   tmp_2271 += tmp_2272;
   result += (-2.6666666666666665*MGlu) * tmp_2271;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpUSeconjUSeconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUSeconjUSeconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2273;
   std::complex<double> tmp_2274;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2274 += A0(MRh(gI1))*CpUSeconjUSeconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2273 += tmp_2274;
   result += (-1) * tmp_2273;
   std::complex<double> tmp_2275;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2276;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2276 += (Conj(CpconjUSebarCha1FvPL(gO2,gI1,gI2))*
            CpconjUSebarCha1FvPL(gO1,gI1,gI2) + Conj(CpconjUSebarCha1FvPR(gO2,gI1,
            gI2))*CpconjUSebarCha1FvPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MFv(gI2));
      }
      tmp_2275 += tmp_2276;
   }
   result += tmp_2275;
   std::complex<double> tmp_2277;
   std::complex<double> tmp_2278;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2279;
      std::complex<double> tmp_2280;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2280 += B0(p,MCha1(gI1),MFv(gI2))*(Conj(CpconjUSebarCha1FvPR
            (gO2,gI1,gI2))*CpconjUSebarCha1FvPL(gO1,gI1,gI2) + Conj(
            CpconjUSebarCha1FvPL(gO2,gI1,gI2))*CpconjUSebarCha1FvPR(gO1,gI1,gI2))*
            MFv(gI2);
      }
      tmp_2279 += tmp_2280;
      tmp_2278 += (MCha1(gI1)) * tmp_2279;
   }
   tmp_2277 += tmp_2278;
   result += (-2) * tmp_2277;
   std::complex<double> tmp_2281;
   std::complex<double> tmp_2282;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2282 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2281 += tmp_2282;
   result += (-1) * tmp_2281;
   std::complex<double> tmp_2283;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2284;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2284 += (Conj(CpconjUSeFvCha2PL(gO2,gI1,gI2))*
            CpconjUSeFvCha2PL(gO1,gI1,gI2) + Conj(CpconjUSeFvCha2PR(gO2,gI1,gI2))*
            CpconjUSeFvCha2PR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha2(gI2));
      }
      tmp_2283 += tmp_2284;
   }
   result += tmp_2283;
   std::complex<double> tmp_2285;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2286;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2286 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_2285 += tmp_2286;
   }
   result += tmp_2285;
   std::complex<double> tmp_2287;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2288;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2288 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_2287 += tmp_2288;
   }
   result += tmp_2287;
   std::complex<double> tmp_2289;
   std::complex<double> tmp_2290;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2291;
      std::complex<double> tmp_2292;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2292 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2291 += tmp_2292;
      tmp_2290 += (MFe(gI1)) * tmp_2291;
   }
   tmp_2289 += tmp_2290;
   result += (-2) * tmp_2289;
   std::complex<double> tmp_2293;
   std::complex<double> tmp_2294;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2295;
      std::complex<double> tmp_2296;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2296 += B0(p,MFv(gI1),MCha2(gI2))*(Conj(CpconjUSeFvCha2PR(
            gO2,gI1,gI2))*CpconjUSeFvCha2PL(gO1,gI1,gI2) + Conj(CpconjUSeFvCha2PL(
            gO2,gI1,gI2))*CpconjUSeFvCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2295 += tmp_2296;
      tmp_2294 += (MFv(gI1)) * tmp_2295;
   }
   tmp_2293 += tmp_2294;
   result += (-2) * tmp_2293;
   std::complex<double> tmp_2297;
   std::complex<double> tmp_2298;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2298 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2297 += tmp_2298;
   result += (-0.5) * tmp_2297;
   std::complex<double> tmp_2299;
   std::complex<double> tmp_2300;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2300 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2299 += tmp_2300;
   result += (-1) * tmp_2299;
   std::complex<double> tmp_2301;
   std::complex<double> tmp_2302;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2302 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2301 += tmp_2302;
   result += (-0.5) * tmp_2301;
   std::complex<double> tmp_2303;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2304;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2304 += (Conj(CpconjUSebarChiFePL(gO2,gI1,gI2))*
            CpconjUSebarChiFePL(gO1,gI1,gI2) + Conj(CpconjUSebarChiFePR(gO2,gI1,
            gI2))*CpconjUSebarChiFePR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFe(gI2));
      }
      tmp_2303 += tmp_2304;
   }
   result += tmp_2303;
   std::complex<double> tmp_2305;
   std::complex<double> tmp_2306;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2307;
      std::complex<double> tmp_2308;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2308 += B0(p,MChi(gI1),MFe(gI2))*(Conj(CpconjUSebarChiFePR(
            gO2,gI1,gI2))*CpconjUSebarChiFePL(gO1,gI1,gI2) + Conj(
            CpconjUSebarChiFePL(gO2,gI1,gI2))*CpconjUSebarChiFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2307 += tmp_2308;
      tmp_2306 += (MChi(gI1)) * tmp_2307;
   }
   tmp_2305 += tmp_2306;
   result += (-2) * tmp_2305;
   std::complex<double> tmp_2309;
   std::complex<double> tmp_2310;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2310 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2309 += tmp_2310;
   result += (-3) * tmp_2309;
   std::complex<double> tmp_2311;
   std::complex<double> tmp_2312;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2312 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2311 += tmp_2312;
   result += (-1) * tmp_2311;
   std::complex<double> tmp_2313;
   std::complex<double> tmp_2314;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2314 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2313 += tmp_2314;
   result += (-3) * tmp_2313;
   std::complex<double> tmp_2315;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2316;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2316 += B0(p,MSe(gI1),MRh(gI2))*Conj(CpconjUSeSeRh(gO2,gI1,
            gI2))*CpconjUSeSeRh(gO1,gI1,gI2);
      }
      tmp_2315 += tmp_2316;
   }
   result += tmp_2315;
   std::complex<double> tmp_2317;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2318;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2318 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2317 += tmp_2318;
   }
   result += tmp_2317;
   std::complex<double> tmp_2319;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2320;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2320 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2319 += tmp_2320;
   }
   result += tmp_2319;
   std::complex<double> tmp_2321;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2321 += B0(p,MSRdp,MSv(gI2))*Conj(CpconjUSeconjSRdpSv(gO2,gI2))*
         CpconjUSeconjSRdpSv(gO1,gI2);
   }
   result += tmp_2321;
   std::complex<double> tmp_2322;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2322 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2322;
   std::complex<double> tmp_2323;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2323 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2323;
   std::complex<double> tmp_2324;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2324 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2324;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmCgWmC(gO1)*CpUhhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmgWm(gO1)*CpUhhbargWmgWm(gO2));
   result += -(B0(p,MVZ,MVZ)*CpUhhbargZgZ(gO1)*CpUhhbargZgZ(gO2));
   result += B0(p,MSRdp,MSRdp)*Conj(CpUhhconjSRdpSRdp(gO2))*CpUhhconjSRdpSRdp(
      gO1);
   result += B0(p,MSRum,MSRum)*Conj(CpUhhconjSRumSRum(gO2))*CpUhhconjSRumSRum(
      gO1);
   result += 4*B0(p,MVWm,MVWm)*Conj(CpUhhconjVWmVWm(gO2))*CpUhhconjVWmVWm(gO1);
   result += -(A0(MSRdp)*CpUhhUhhconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUhhUhhconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(p,MVZ,MVZ)*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   std::complex<double> tmp_2325;
   std::complex<double> tmp_2326;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2326 += A0(MRh(gI1))*CpUhhUhhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2325 += tmp_2326;
   result += (-1) * tmp_2325;
   std::complex<double> tmp_2327;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2328;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2328 += B0(p,MRh(gI1),MRh(gI2))*Conj(CpUhhconjRhRh(gO2,gI1,
            gI2))*CpUhhconjRhRh(gO1,gI1,gI2);
      }
      tmp_2327 += tmp_2328;
   }
   result += tmp_2327;
   std::complex<double> tmp_2329;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2330;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2330 += (Conj(CpUhhbarCha1Cha1PL(gO2,gI1,gI2))*
            CpUhhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(CpUhhbarCha1Cha1PR(gO2,gI1,gI2)
            )*CpUhhbarCha1Cha1PR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MCha1(gI2));
      }
      tmp_2329 += tmp_2330;
   }
   result += tmp_2329;
   std::complex<double> tmp_2331;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2332;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2332 += (Conj(CpUhhbarCha2Cha2PL(gO2,gI1,gI2))*
            CpUhhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(CpUhhbarCha2Cha2PR(gO2,gI1,gI2)
            )*CpUhhbarCha2Cha2PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha2(gI2));
      }
      tmp_2331 += tmp_2332;
   }
   result += tmp_2331;
   std::complex<double> tmp_2333;
   std::complex<double> tmp_2334;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2335;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2335 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpUhhconjRhAh(gO2,gI1,
            gI2))*CpUhhconjRhAh(gO1,gI1,gI2);
      }
      tmp_2334 += tmp_2335;
   }
   tmp_2333 += tmp_2334;
   result += (2) * tmp_2333;
   std::complex<double> tmp_2336;
   std::complex<double> tmp_2337;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2338;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2338 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpUhhconjRhhh(gO2,gI1,
            gI2))*CpUhhconjRhhh(gO1,gI1,gI2);
      }
      tmp_2337 += tmp_2338;
   }
   tmp_2336 += tmp_2337;
   result += (2) * tmp_2336;
   std::complex<double> tmp_2339;
   std::complex<double> tmp_2340;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2341;
      std::complex<double> tmp_2342;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2342 += B0(p,MCha1(gI1),MCha1(gI2))*(Conj(CpUhhbarCha1Cha1PR
            (gO2,gI1,gI2))*CpUhhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(
            CpUhhbarCha1Cha1PL(gO2,gI1,gI2))*CpUhhbarCha1Cha1PR(gO1,gI1,gI2))*
            MCha1(gI2);
      }
      tmp_2341 += tmp_2342;
      tmp_2340 += (MCha1(gI1)) * tmp_2341;
   }
   tmp_2339 += tmp_2340;
   result += (-2) * tmp_2339;
   std::complex<double> tmp_2343;
   std::complex<double> tmp_2344;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2345;
      std::complex<double> tmp_2346;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2346 += B0(p,MCha2(gI1),MCha2(gI2))*(Conj(CpUhhbarCha2Cha2PR
            (gO2,gI1,gI2))*CpUhhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(
            CpUhhbarCha2Cha2PL(gO2,gI1,gI2))*CpUhhbarCha2Cha2PR(gO1,gI1,gI2))*
            MCha2(gI2);
      }
      tmp_2345 += tmp_2346;
      tmp_2344 += (MCha2(gI1)) * tmp_2345;
   }
   tmp_2343 += tmp_2344;
   result += (-2) * tmp_2343;
   std::complex<double> tmp_2347;
   std::complex<double> tmp_2348;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2348 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2347 += tmp_2348;
   result += (-1) * tmp_2347;
   std::complex<double> tmp_2349;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2350;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2350 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2349 += tmp_2350;
   }
   result += tmp_2349;
   std::complex<double> tmp_2351;
   std::complex<double> tmp_2352;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2353;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2353 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2352 += tmp_2353;
   }
   tmp_2351 += tmp_2352;
   result += (3) * tmp_2351;
   std::complex<double> tmp_2354;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2355;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2355 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2354 += tmp_2355;
   }
   result += tmp_2354;
   std::complex<double> tmp_2356;
   std::complex<double> tmp_2357;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2358;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2358 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2357 += tmp_2358;
   }
   tmp_2356 += tmp_2357;
   result += (3) * tmp_2356;
   std::complex<double> tmp_2359;
   std::complex<double> tmp_2360;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2361;
      std::complex<double> tmp_2362;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2362 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2361 += tmp_2362;
      tmp_2360 += (MFd(gI1)) * tmp_2361;
   }
   tmp_2359 += tmp_2360;
   result += (-6) * tmp_2359;
   std::complex<double> tmp_2363;
   std::complex<double> tmp_2364;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2365;
      std::complex<double> tmp_2366;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2366 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2365 += tmp_2366;
      tmp_2364 += (MFe(gI1)) * tmp_2365;
   }
   tmp_2363 += tmp_2364;
   result += (-2) * tmp_2363;
   std::complex<double> tmp_2367;
   std::complex<double> tmp_2368;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2369;
      std::complex<double> tmp_2370;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2370 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2369 += tmp_2370;
      tmp_2368 += (MFu(gI1)) * tmp_2369;
   }
   tmp_2367 += tmp_2368;
   result += (-6) * tmp_2367;
   std::complex<double> tmp_2371;
   std::complex<double> tmp_2372;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2372 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2371 += tmp_2372;
   result += (-0.5) * tmp_2371;
   std::complex<double> tmp_2373;
   std::complex<double> tmp_2374;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2374 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2373 += tmp_2374;
   result += (-1) * tmp_2373;
   std::complex<double> tmp_2375;
   std::complex<double> tmp_2376;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2376 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2375 += tmp_2376;
   result += (-0.5) * tmp_2375;
   std::complex<double> tmp_2377;
   std::complex<double> tmp_2378;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2379;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2379 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2378 += tmp_2379;
   }
   tmp_2377 += tmp_2378;
   result += (0.5) * tmp_2377;
   std::complex<double> tmp_2380;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2381;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2381 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2380 += tmp_2381;
   }
   result += tmp_2380;
   std::complex<double> tmp_2382;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2383;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2383 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_2382 += tmp_2383;
   }
   result += tmp_2382;
   std::complex<double> tmp_2384;
   std::complex<double> tmp_2385;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2386;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2386 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2385 += tmp_2386;
   }
   tmp_2384 += tmp_2385;
   result += (0.5) * tmp_2384;
   std::complex<double> tmp_2387;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2388;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2388 += (Conj(CpUhhbarChiChiPL(gO2,gI1,gI2))*
            CpUhhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUhhbarChiChiPR(gO2,gI1,gI2))*
            CpUhhbarChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2387 += tmp_2388;
   }
   result += tmp_2387;
   std::complex<double> tmp_2389;
   std::complex<double> tmp_2390;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2391;
      std::complex<double> tmp_2392;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2392 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhbarChiChiPR(gO2
            ,gI1,gI2))*CpUhhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUhhbarChiChiPL(gO2,
            gI1,gI2))*CpUhhbarChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2391 += tmp_2392;
      tmp_2390 += (MChi(gI1)) * tmp_2391;
   }
   tmp_2389 += tmp_2390;
   result += (-2) * tmp_2389;
   std::complex<double> tmp_2393;
   std::complex<double> tmp_2394;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2394 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2393 += tmp_2394;
   result += (-3) * tmp_2393;
   std::complex<double> tmp_2395;
   std::complex<double> tmp_2396;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2396 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2395 += tmp_2396;
   result += (-1) * tmp_2395;
   std::complex<double> tmp_2397;
   std::complex<double> tmp_2398;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2398 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2397 += tmp_2398;
   result += (-3) * tmp_2397;
   std::complex<double> tmp_2399;
   std::complex<double> tmp_2400;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2401;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2401 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2400 += tmp_2401;
   }
   tmp_2399 += tmp_2400;
   result += (3) * tmp_2399;
   std::complex<double> tmp_2402;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2403;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2403 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2402 += tmp_2403;
   }
   result += tmp_2402;
   std::complex<double> tmp_2404;
   std::complex<double> tmp_2405;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2406;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2406 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2405 += tmp_2406;
   }
   tmp_2404 += tmp_2405;
   result += (3) * tmp_2404;
   std::complex<double> tmp_2407;
   std::complex<double> tmp_2408;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2408 += B0(p,MSRum,MHpm(gI2))*Conj(CpUhhconjSRumHpm(gO2,gI2))*
         CpUhhconjSRumHpm(gO1,gI2);
   }
   tmp_2407 += tmp_2408;
   result += (2) * tmp_2407;
   std::complex<double> tmp_2409;
   std::complex<double> tmp_2410;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2410 += B0(p,MSRdp,MHpm(gI2))*Conj(CpUhhSRdpHpm(gO2,gI2))*
         CpUhhSRdpHpm(gO1,gI2);
   }
   tmp_2409 += tmp_2410;
   result += (2) * tmp_2409;
   std::complex<double> tmp_2411;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2411 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2411;
   std::complex<double> tmp_2412;
   std::complex<double> tmp_2413;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2413 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2412 += tmp_2413;
   result += (2) * tmp_2412;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += B0(p,MSRdp,MSRdp)*Conj(CpUAhconjSRdpSRdp(gO2))*CpUAhconjSRdpSRdp(
      gO1);
   result += B0(p,MSRum,MSRum)*Conj(CpUAhconjSRumSRum(gO2))*CpUAhconjSRumSRum(
      gO1);
   result += -(A0(MSRdp)*CpUAhUAhconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUAhUAhconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2414;
   std::complex<double> tmp_2415;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2415 += A0(MRh(gI1))*CpUAhUAhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2414 += tmp_2415;
   result += (-1) * tmp_2414;
   std::complex<double> tmp_2416;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2417;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2417 += B0(p,MRh(gI1),MRh(gI2))*Conj(CpUAhconjRhRh(gO2,gI1,
            gI2))*CpUAhconjRhRh(gO1,gI1,gI2);
      }
      tmp_2416 += tmp_2417;
   }
   result += tmp_2416;
   std::complex<double> tmp_2418;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2419;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2419 += (Conj(CpUAhbarCha1Cha1PL(gO2,gI1,gI2))*
            CpUAhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(CpUAhbarCha1Cha1PR(gO2,gI1,gI2)
            )*CpUAhbarCha1Cha1PR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MCha1(gI2));
      }
      tmp_2418 += tmp_2419;
   }
   result += tmp_2418;
   std::complex<double> tmp_2420;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2421;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2421 += (Conj(CpUAhbarCha2Cha2PL(gO2,gI1,gI2))*
            CpUAhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(CpUAhbarCha2Cha2PR(gO2,gI1,gI2)
            )*CpUAhbarCha2Cha2PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha2(gI2));
      }
      tmp_2420 += tmp_2421;
   }
   result += tmp_2420;
   std::complex<double> tmp_2422;
   std::complex<double> tmp_2423;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2424;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2424 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpUAhconjRhAh(gO2,gI1,
            gI2))*CpUAhconjRhAh(gO1,gI1,gI2);
      }
      tmp_2423 += tmp_2424;
   }
   tmp_2422 += tmp_2423;
   result += (2) * tmp_2422;
   std::complex<double> tmp_2425;
   std::complex<double> tmp_2426;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2427;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2427 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpUAhconjRhhh(gO2,gI1,
            gI2))*CpUAhconjRhhh(gO1,gI1,gI2);
      }
      tmp_2426 += tmp_2427;
   }
   tmp_2425 += tmp_2426;
   result += (2) * tmp_2425;
   std::complex<double> tmp_2428;
   std::complex<double> tmp_2429;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2430;
      std::complex<double> tmp_2431;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2431 += B0(p,MCha1(gI1),MCha1(gI2))*(Conj(CpUAhbarCha1Cha1PR
            (gO2,gI1,gI2))*CpUAhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(
            CpUAhbarCha1Cha1PL(gO2,gI1,gI2))*CpUAhbarCha1Cha1PR(gO1,gI1,gI2))*
            MCha1(gI2);
      }
      tmp_2430 += tmp_2431;
      tmp_2429 += (MCha1(gI1)) * tmp_2430;
   }
   tmp_2428 += tmp_2429;
   result += (-2) * tmp_2428;
   std::complex<double> tmp_2432;
   std::complex<double> tmp_2433;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2434;
      std::complex<double> tmp_2435;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2435 += B0(p,MCha2(gI1),MCha2(gI2))*(Conj(CpUAhbarCha2Cha2PR
            (gO2,gI1,gI2))*CpUAhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(
            CpUAhbarCha2Cha2PL(gO2,gI1,gI2))*CpUAhbarCha2Cha2PR(gO1,gI1,gI2))*
            MCha2(gI2);
      }
      tmp_2434 += tmp_2435;
      tmp_2433 += (MCha2(gI1)) * tmp_2434;
   }
   tmp_2432 += tmp_2433;
   result += (-2) * tmp_2432;
   std::complex<double> tmp_2436;
   std::complex<double> tmp_2437;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2437 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2436 += tmp_2437;
   result += (-1) * tmp_2436;
   std::complex<double> tmp_2438;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2439;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2439 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUAhconjSvSv(gO2,gI1,
            gI2))*CpUAhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2438 += tmp_2439;
   }
   result += tmp_2438;
   std::complex<double> tmp_2440;
   std::complex<double> tmp_2441;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2442;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2442 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2441 += tmp_2442;
   }
   tmp_2440 += tmp_2441;
   result += (3) * tmp_2440;
   std::complex<double> tmp_2443;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2444;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2444 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2443 += tmp_2444;
   }
   result += tmp_2443;
   std::complex<double> tmp_2445;
   std::complex<double> tmp_2446;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2447;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2447 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2446 += tmp_2447;
   }
   tmp_2445 += tmp_2446;
   result += (3) * tmp_2445;
   std::complex<double> tmp_2448;
   std::complex<double> tmp_2449;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2450;
      std::complex<double> tmp_2451;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2451 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2450 += tmp_2451;
      tmp_2449 += (MFd(gI1)) * tmp_2450;
   }
   tmp_2448 += tmp_2449;
   result += (-6) * tmp_2448;
   std::complex<double> tmp_2452;
   std::complex<double> tmp_2453;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2454;
      std::complex<double> tmp_2455;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2455 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2454 += tmp_2455;
      tmp_2453 += (MFe(gI1)) * tmp_2454;
   }
   tmp_2452 += tmp_2453;
   result += (-2) * tmp_2452;
   std::complex<double> tmp_2456;
   std::complex<double> tmp_2457;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2458;
      std::complex<double> tmp_2459;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2459 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2458 += tmp_2459;
      tmp_2457 += (MFu(gI1)) * tmp_2458;
   }
   tmp_2456 += tmp_2457;
   result += (-6) * tmp_2456;
   std::complex<double> tmp_2460;
   std::complex<double> tmp_2461;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2461 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2460 += tmp_2461;
   result += (-0.5) * tmp_2460;
   std::complex<double> tmp_2462;
   std::complex<double> tmp_2463;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2463 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2462 += tmp_2463;
   result += (-1) * tmp_2462;
   std::complex<double> tmp_2464;
   std::complex<double> tmp_2465;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2465 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2464 += tmp_2465;
   result += (-0.5) * tmp_2464;
   std::complex<double> tmp_2466;
   std::complex<double> tmp_2467;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2468;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2468 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_2467 += tmp_2468;
   }
   tmp_2466 += tmp_2467;
   result += (0.5) * tmp_2466;
   std::complex<double> tmp_2469;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2470;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2470 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2469 += tmp_2470;
   }
   result += tmp_2469;
   std::complex<double> tmp_2471;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2472;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2472 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2471 += tmp_2472;
   }
   result += tmp_2471;
   std::complex<double> tmp_2473;
   std::complex<double> tmp_2474;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2475;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2475 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_2474 += tmp_2475;
   }
   tmp_2473 += tmp_2474;
   result += (0.5) * tmp_2473;
   std::complex<double> tmp_2476;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2477;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2477 += (Conj(CpUAhbarChiChiPL(gO2,gI1,gI2))*
            CpUAhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUAhbarChiChiPR(gO2,gI1,gI2))*
            CpUAhbarChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2476 += tmp_2477;
   }
   result += tmp_2476;
   std::complex<double> tmp_2478;
   std::complex<double> tmp_2479;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2480;
      std::complex<double> tmp_2481;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2481 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhbarChiChiPR(gO2
            ,gI1,gI2))*CpUAhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUAhbarChiChiPL(gO2,
            gI1,gI2))*CpUAhbarChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2480 += tmp_2481;
      tmp_2479 += (MChi(gI1)) * tmp_2480;
   }
   tmp_2478 += tmp_2479;
   result += (-2) * tmp_2478;
   std::complex<double> tmp_2482;
   std::complex<double> tmp_2483;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2483 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2482 += tmp_2483;
   result += (-3) * tmp_2482;
   std::complex<double> tmp_2484;
   std::complex<double> tmp_2485;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2485 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2484 += tmp_2485;
   result += (-1) * tmp_2484;
   std::complex<double> tmp_2486;
   std::complex<double> tmp_2487;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2487 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2486 += tmp_2487;
   result += (-3) * tmp_2486;
   std::complex<double> tmp_2488;
   std::complex<double> tmp_2489;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2490;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2490 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2489 += tmp_2490;
   }
   tmp_2488 += tmp_2489;
   result += (3) * tmp_2488;
   std::complex<double> tmp_2491;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2492;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2492 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2491 += tmp_2492;
   }
   result += tmp_2491;
   std::complex<double> tmp_2493;
   std::complex<double> tmp_2494;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2495;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2495 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2494 += tmp_2495;
   }
   tmp_2493 += tmp_2494;
   result += (3) * tmp_2493;
   std::complex<double> tmp_2496;
   std::complex<double> tmp_2497;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2497 += B0(p,MSRum,MHpm(gI2))*Conj(CpUAhconjSRumHpm(gO2,gI2))*
         CpUAhconjSRumHpm(gO1,gI2);
   }
   tmp_2496 += tmp_2497;
   result += (2) * tmp_2496;
   std::complex<double> tmp_2498;
   std::complex<double> tmp_2499;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2499 += B0(p,MSRdp,MHpm(gI2))*Conj(CpUAhSRdpHpm(gO2,gI2))*
         CpUAhSRdpHpm(gO1,gI2);
   }
   tmp_2498 += tmp_2499;
   result += (2) * tmp_2498;
   std::complex<double> tmp_2500;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2500 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2500;
   std::complex<double> tmp_2501;
   std::complex<double> tmp_2502;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2502 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2501 += tmp_2502;
   result += (2) * tmp_2501;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Rh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpURhconjURhconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpURhconjURhconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpURhconjURhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpURhconjURhVZVZ(gO1,gO2);
   result += Conj(CpconjURhVWmSRdp(gO2))*CpconjURhVWmSRdp(gO1)*F0(p,MSRdp,MVWm)
      ;
   result += Conj(CpconjURhconjVWmSRum(gO2))*CpconjURhconjVWmSRum(gO1)*F0(p,
      MSRum,MVWm);
   std::complex<double> tmp_2503;
   std::complex<double> tmp_2504;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2504 += A0(MRh(gI1))*CpURhconjURhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2503 += tmp_2504;
   result += (-1) * tmp_2503;
   std::complex<double> tmp_2505;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2506;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2506 += (Conj(CpconjURhCha2Cha1PL(gO2,gI1,gI2))*
            CpconjURhCha2Cha1PL(gO1,gI1,gI2) + Conj(CpconjURhCha2Cha1PR(gO2,gI1,
            gI2))*CpconjURhCha2Cha1PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha1(gI2));
      }
      tmp_2505 += tmp_2506;
   }
   result += tmp_2505;
   std::complex<double> tmp_2507;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2508;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2508 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpconjURhRhAh(gO2,gI1,
            gI2))*CpconjURhRhAh(gO1,gI1,gI2);
      }
      tmp_2507 += tmp_2508;
   }
   result += tmp_2507;
   std::complex<double> tmp_2509;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2510;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2510 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpconjURhRhhh(gO2,gI1,
            gI2))*CpconjURhRhhh(gO1,gI1,gI2);
      }
      tmp_2509 += tmp_2510;
   }
   result += tmp_2509;
   std::complex<double> tmp_2511;
   std::complex<double> tmp_2512;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2513;
      std::complex<double> tmp_2514;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2514 += B0(p,MCha2(gI1),MCha1(gI2))*(Conj(
            CpconjURhCha2Cha1PR(gO2,gI1,gI2))*CpconjURhCha2Cha1PL(gO1,gI1,gI2) +
            Conj(CpconjURhCha2Cha1PL(gO2,gI1,gI2))*CpconjURhCha2Cha1PR(gO1,gI1,gI2
            ))*MCha1(gI2);
      }
      tmp_2513 += tmp_2514;
      tmp_2512 += (MCha2(gI1)) * tmp_2513;
   }
   tmp_2511 += tmp_2512;
   result += (-2) * tmp_2511;
   std::complex<double> tmp_2515;
   std::complex<double> tmp_2516;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2516 += A0(MSv(gI1))*CpURhconjURhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2515 += tmp_2516;
   result += (-1) * tmp_2515;
   std::complex<double> tmp_2517;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2517 += B0(p,MHpm(gI1),MSRum)*Conj(CpconjURhconjHpmSRum(gO2,gI1))*
         CpconjURhconjHpmSRum(gO1,gI1);
   }
   result += tmp_2517;
   std::complex<double> tmp_2518;
   std::complex<double> tmp_2519;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2519 += A0(MAh(gI1))*CpURhconjURhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2518 += tmp_2519;
   result += (-0.5) * tmp_2518;
   std::complex<double> tmp_2520;
   std::complex<double> tmp_2521;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2521 += A0(MHpm(gI1))*CpURhconjURhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2520 += tmp_2521;
   result += (-1) * tmp_2520;
   std::complex<double> tmp_2522;
   std::complex<double> tmp_2523;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2523 += A0(Mhh(gI1))*CpURhconjURhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2522 += tmp_2523;
   result += (-0.5) * tmp_2522;
   std::complex<double> tmp_2524;
   std::complex<double> tmp_2525;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2526;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2526 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpconjURhAhAh(gO2,gI1,
            gI2))*CpconjURhAhAh(gO1,gI1,gI2);
      }
      tmp_2525 += tmp_2526;
   }
   tmp_2524 += tmp_2525;
   result += (0.25) * tmp_2524;
   std::complex<double> tmp_2527;
   std::complex<double> tmp_2528;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2529;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2529 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpconjURhconjHpmHpm(
            gO2,gI1,gI2))*CpconjURhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2528 += tmp_2529;
   }
   tmp_2527 += tmp_2528;
   result += (0.5) * tmp_2527;
   std::complex<double> tmp_2530;
   std::complex<double> tmp_2531;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2532;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2532 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpconjURhhhAh(gO2,gI1,
            gI2))*CpconjURhhhAh(gO1,gI1,gI2);
      }
      tmp_2531 += tmp_2532;
   }
   tmp_2530 += tmp_2531;
   result += (0.5) * tmp_2530;
   std::complex<double> tmp_2533;
   std::complex<double> tmp_2534;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2535;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2535 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpconjURhhhhh(gO2,gI1,
            gI2))*CpconjURhhhhh(gO1,gI1,gI2);
      }
      tmp_2534 += tmp_2535;
   }
   tmp_2533 += tmp_2534;
   result += (0.25) * tmp_2533;
   std::complex<double> tmp_2536;
   std::complex<double> tmp_2537;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2538;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2538 += (Conj(CpconjURhChiChiPL(gO2,gI1,gI2))*
            CpconjURhChiChiPL(gO1,gI1,gI2) + Conj(CpconjURhChiChiPR(gO2,gI1,gI2))*
            CpconjURhChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2537 += tmp_2538;
   }
   tmp_2536 += tmp_2537;
   result += (0.5) * tmp_2536;
   std::complex<double> tmp_2539;
   std::complex<double> tmp_2540;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2541;
      std::complex<double> tmp_2542;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2542 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpconjURhChiChiPR(
            gO2,gI1,gI2))*CpconjURhChiChiPL(gO1,gI1,gI2) + Conj(CpconjURhChiChiPL(
            gO2,gI1,gI2))*CpconjURhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2541 += tmp_2542;
      tmp_2540 += (MChi(gI1)) * tmp_2541;
   }
   tmp_2539 += tmp_2540;
   result += (-1) * tmp_2539;
   std::complex<double> tmp_2543;
   std::complex<double> tmp_2544;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2544 += A0(MSd(gI1))*CpURhconjURhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2543 += tmp_2544;
   result += (-3) * tmp_2543;
   std::complex<double> tmp_2545;
   std::complex<double> tmp_2546;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2546 += A0(MSe(gI1))*CpURhconjURhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2545 += tmp_2546;
   result += (-1) * tmp_2545;
   std::complex<double> tmp_2547;
   std::complex<double> tmp_2548;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2548 += A0(MSu(gI1))*CpURhconjURhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2547 += tmp_2548;
   result += (-3) * tmp_2547;
   std::complex<double> tmp_2549;
   std::complex<double> tmp_2550;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2551;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2551 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpconjURhconjSdSd(gO2,
            gI1,gI2))*CpconjURhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2550 += tmp_2551;
   }
   tmp_2549 += tmp_2550;
   result += (1.5) * tmp_2549;
   std::complex<double> tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2554;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2554 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpconjURhconjSeSe(gO2,
            gI1,gI2))*CpconjURhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2553 += tmp_2554;
   }
   tmp_2552 += tmp_2553;
   result += (0.5) * tmp_2552;
   std::complex<double> tmp_2555;
   std::complex<double> tmp_2556;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2557;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2557 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpconjURhconjSuSu(gO2,
            gI1,gI2))*CpconjURhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2556 += tmp_2557;
   }
   tmp_2555 += tmp_2556;
   result += (1.5) * tmp_2555;
   std::complex<double> tmp_2558;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2558 += Conj(CpconjURhVZRh(gO2,gI2))*CpconjURhVZRh(gO1,gI2)*F0(p,
         MRh(gI2),MVZ);
   }
   result += tmp_2558;
   std::complex<double> tmp_2559;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2559 += B0(p,MSRdp,MHpm(gI2))*Conj(CpconjURhSRdpHpm(gO2,gI2))*
         CpconjURhSRdpHpm(gO1,gI2);
   }
   result += tmp_2559;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*B0(p,0,MVWm)*Conj(CpconjUHpmVWmVP(gO2))*CpconjUHpmVWmVP(gO1);
   result += 4*B0(p,MVWm,MVZ)*Conj(CpconjUHpmVZVWm(gO2))*CpconjUHpmVZVWm(gO1);
   result += -(A0(MSRdp)*CpUHpmconjUHpmconjSRdpSRdp(gO1,gO2));
   result += -(A0(MSRum)*CpUHpmconjUHpmconjSRumSRum(gO1,gO2));
   result += 4*A0(MVWm)*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += -(B0(p,MVZ,MVWm)*CpconjUHpmbargWmCgZ(gO1)*CpUHpmgWmCbargZ(gO2));
   result += -(B0(p,MVWm,MVZ)*CpconjUHpmbargZgWm(gO1)*CpUHpmgZbargWm(gO2));
   std::complex<double> tmp_2560;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2560 += B0(p,MRh(gI1),MSRum)*Conj(CpconjUHpmconjRhSRum(gO2,gI1))*
         CpconjUHpmconjRhSRum(gO1,gI1);
   }
   result += tmp_2560;
   std::complex<double> tmp_2561;
   std::complex<double> tmp_2562;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2562 += A0(MRh(gI1))*CpUHpmconjUHpmconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2561 += tmp_2562;
   result += (-1) * tmp_2561;
   std::complex<double> tmp_2563;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2564;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2564 += B0(p,MRh(gI1),MHpm(gI2))*Conj(CpconjUHpmconjRhHpm(
            gO2,gI1,gI2))*CpconjUHpmconjRhHpm(gO1,gI1,gI2);
      }
      tmp_2563 += tmp_2564;
   }
   result += tmp_2563;
   std::complex<double> tmp_2565;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2566;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2566 += (Conj(CpconjUHpmbarCha1ChiPL(gO2,gI1,gI2))*
            CpconjUHpmbarCha1ChiPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarCha1ChiPR(gO2,
            gI1,gI2))*CpconjUHpmbarCha1ChiPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MChi(
            gI2));
      }
      tmp_2565 += tmp_2566;
   }
   result += tmp_2565;
   std::complex<double> tmp_2567;
   std::complex<double> tmp_2568;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2569;
      std::complex<double> tmp_2570;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2570 += B0(p,MCha1(gI1),MChi(gI2))*(Conj(
            CpconjUHpmbarCha1ChiPR(gO2,gI1,gI2))*CpconjUHpmbarCha1ChiPL(gO1,gI1,
            gI2) + Conj(CpconjUHpmbarCha1ChiPL(gO2,gI1,gI2))*
            CpconjUHpmbarCha1ChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2569 += tmp_2570;
      tmp_2568 += (MCha1(gI1)) * tmp_2569;
   }
   tmp_2567 += tmp_2568;
   result += (-2) * tmp_2567;
   std::complex<double> tmp_2571;
   std::complex<double> tmp_2572;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2572 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2571 += tmp_2572;
   result += (-1) * tmp_2571;
   std::complex<double> tmp_2573;
   std::complex<double> tmp_2574;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2575;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2575 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2574 += tmp_2575;
   }
   tmp_2573 += tmp_2574;
   result += (3) * tmp_2573;
   std::complex<double> tmp_2576;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2577;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2577 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2576 += tmp_2577;
   }
   result += tmp_2576;
   std::complex<double> tmp_2578;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2579;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2579 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2578 += tmp_2579;
   }
   result += tmp_2578;
   std::complex<double> tmp_2580;
   std::complex<double> tmp_2581;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2582;
      std::complex<double> tmp_2583;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2583 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2582 += tmp_2583;
      tmp_2581 += (MFu(gI1)) * tmp_2582;
   }
   tmp_2580 += tmp_2581;
   result += (-6) * tmp_2580;
   std::complex<double> tmp_2584;
   std::complex<double> tmp_2585;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2586;
      std::complex<double> tmp_2587;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2587 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2586 += tmp_2587;
      tmp_2585 += (MFv(gI1)) * tmp_2586;
   }
   tmp_2584 += tmp_2585;
   result += (-2) * tmp_2584;
   std::complex<double> tmp_2588;
   std::complex<double> tmp_2589;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2589 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2588 += tmp_2589;
   result += (-0.5) * tmp_2588;
   std::complex<double> tmp_2590;
   std::complex<double> tmp_2591;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2591 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2590 += tmp_2591;
   result += (-1) * tmp_2590;
   std::complex<double> tmp_2592;
   std::complex<double> tmp_2593;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2593 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2592 += tmp_2593;
   result += (-0.5) * tmp_2592;
   std::complex<double> tmp_2594;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2595;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2595 += (Conj(CpconjUHpmbarChiCha2PL(gO2,gI1,gI2))*
            CpconjUHpmbarChiCha2PL(gO1,gI1,gI2) + Conj(CpconjUHpmbarChiCha2PR(gO2,
            gI1,gI2))*CpconjUHpmbarChiCha2PR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha2(
            gI2));
      }
      tmp_2594 += tmp_2595;
   }
   result += tmp_2594;
   std::complex<double> tmp_2596;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2597;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2597 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2596 += tmp_2597;
   }
   result += tmp_2596;
   std::complex<double> tmp_2598;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2599;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2599 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2598 += tmp_2599;
   }
   result += tmp_2598;
   std::complex<double> tmp_2600;
   std::complex<double> tmp_2601;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2602;
      std::complex<double> tmp_2603;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2603 += B0(p,MChi(gI1),MCha2(gI2))*(Conj(
            CpconjUHpmbarChiCha2PR(gO2,gI1,gI2))*CpconjUHpmbarChiCha2PL(gO1,gI1,
            gI2) + Conj(CpconjUHpmbarChiCha2PL(gO2,gI1,gI2))*
            CpconjUHpmbarChiCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2602 += tmp_2603;
      tmp_2601 += (MChi(gI1)) * tmp_2602;
   }
   tmp_2600 += tmp_2601;
   result += (-2) * tmp_2600;
   std::complex<double> tmp_2604;
   std::complex<double> tmp_2605;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2605 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2604 += tmp_2605;
   result += (-3) * tmp_2604;
   std::complex<double> tmp_2606;
   std::complex<double> tmp_2607;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2607 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2606 += tmp_2607;
   result += (-1) * tmp_2606;
   std::complex<double> tmp_2608;
   std::complex<double> tmp_2609;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2609 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2608 += tmp_2609;
   result += (-3) * tmp_2608;
   std::complex<double> tmp_2610;
   std::complex<double> tmp_2611;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2612;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2612 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2611 += tmp_2612;
   }
   tmp_2610 += tmp_2611;
   result += (3) * tmp_2610;
   std::complex<double> tmp_2613;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2613 += B0(p,MSRdp,MRh(gI2))*Conj(CpconjUHpmconjSRdpRh(gO2,gI2))*
         CpconjUHpmconjSRdpRh(gO1,gI2);
   }
   result += tmp_2613;
   std::complex<double> tmp_2614;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2614 += B0(p,MSRdp,MAh(gI2))*Conj(CpconjUHpmconjSRdpAh(gO2,gI2))*
         CpconjUHpmconjSRdpAh(gO1,gI2);
   }
   result += tmp_2614;
   std::complex<double> tmp_2615;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2615 += B0(p,MSRdp,Mhh(gI2))*Conj(CpconjUHpmconjSRdphh(gO2,gI2))*
         CpconjUHpmconjSRdphh(gO1,gI2);
   }
   result += tmp_2615;
   std::complex<double> tmp_2616;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2616 += B0(p,MSRum,MAh(gI2))*Conj(CpconjUHpmSRumAh(gO2,gI2))*
         CpconjUHpmSRumAh(gO1,gI2);
   }
   result += tmp_2616;
   std::complex<double> tmp_2617;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2617 += B0(p,MSRum,Mhh(gI2))*Conj(CpconjUHpmSRumhh(gO2,gI2))*
         CpconjUHpmSRumhh(gO1,gI2);
   }
   result += tmp_2617;
   std::complex<double> tmp_2618;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2618 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2618;
   std::complex<double> tmp_2619;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2619 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2619;
   std::complex<double> tmp_2620;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2620 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2620;
   std::complex<double> tmp_2621;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2621 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2621;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_SRdp(double p ) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpSRdpconjSRdpconjSRdpSRdp());
   result += -(A0(MSRum)*CpSRdpconjSRdpconjSRumSRum());
   result += 4*A0(MVWm)*CpSRdpconjSRdpconjVWmVWm();
   result += 2*A0(MVZ)*CpSRdpconjSRdpVZVZ();
   result += AbsSqr(CpconjSRdpVPSRdp())*F0(p,MSRdp,0);
   result += AbsSqr(CpconjSRdpVZSRdp())*F0(p,MSRdp,MVZ);
   std::complex<double> tmp_2622;
   std::complex<double> tmp_2623;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2623 += A0(MRh(gI1))*CpSRdpconjSRdpconjRhRh(gI1,gI1);
   }
   tmp_2622 += tmp_2623;
   result += (-1) * tmp_2622;
   std::complex<double> tmp_2624;
   std::complex<double> tmp_2625;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2625 += A0(MSv(gI1))*CpSRdpconjSRdpconjSvSv(gI1,gI1);
   }
   tmp_2624 += tmp_2625;
   result += (-1) * tmp_2624;
   std::complex<double> tmp_2626;
   std::complex<double> tmp_2627;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2627 += A0(MAh(gI1))*CpSRdpconjSRdpAhAh(gI1,gI1);
   }
   tmp_2626 += tmp_2627;
   result += (-0.5) * tmp_2626;
   std::complex<double> tmp_2628;
   std::complex<double> tmp_2629;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2629 += A0(MHpm(gI1))*CpSRdpconjSRdpconjHpmHpm(gI1,gI1);
   }
   tmp_2628 += tmp_2629;
   result += (-1) * tmp_2628;
   std::complex<double> tmp_2630;
   std::complex<double> tmp_2631;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2631 += A0(Mhh(gI1))*CpSRdpconjSRdphhhh(gI1,gI1);
   }
   tmp_2630 += tmp_2631;
   result += (-0.5) * tmp_2630;
   std::complex<double> tmp_2632;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2633;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2633 += AbsSqr(CpconjSRdpconjHpmRh(gI1,gI2))*B0(p,MHpm(gI1),
            MRh(gI2));
      }
      tmp_2632 += tmp_2633;
   }
   result += tmp_2632;
   std::complex<double> tmp_2634;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2635;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2635 += (AbsSqr(CpconjSRdpChiCha1PL(gI1,gI2)) + AbsSqr(
            CpconjSRdpChiCha1PR(gI1,gI2)))*G0(p,MChi(gI1),MCha1(gI2));
      }
      tmp_2634 += tmp_2635;
   }
   result += tmp_2634;
   std::complex<double> tmp_2636;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2637;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2637 += AbsSqr(CpconjSRdpconjHpmAh(gI1,gI2))*B0(p,MHpm(gI1),
            MAh(gI2));
      }
      tmp_2636 += tmp_2637;
   }
   result += tmp_2636;
   std::complex<double> tmp_2638;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2639;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2639 += AbsSqr(CpconjSRdpconjHpmhh(gI1,gI2))*B0(p,MHpm(gI1),
            Mhh(gI2));
      }
      tmp_2638 += tmp_2639;
   }
   result += tmp_2638;
   std::complex<double> tmp_2640;
   std::complex<double> tmp_2641;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2642;
      std::complex<double> tmp_2643;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2643 += B0(p,MChi(gI1),MCha1(gI2))*(Conj(CpconjSRdpChiCha1PR
            (gI1,gI2))*CpconjSRdpChiCha1PL(gI1,gI2) + Conj(CpconjSRdpChiCha1PL(gI1
            ,gI2))*CpconjSRdpChiCha1PR(gI1,gI2))*MCha1(gI2);
      }
      tmp_2642 += tmp_2643;
      tmp_2641 += (MChi(gI1)) * tmp_2642;
   }
   tmp_2640 += tmp_2641;
   result += (-2) * tmp_2640;
   std::complex<double> tmp_2644;
   std::complex<double> tmp_2645;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2645 += A0(MSd(gI1))*CpSRdpconjSRdpconjSdSd(gI1,gI1);
   }
   tmp_2644 += tmp_2645;
   result += (-3) * tmp_2644;
   std::complex<double> tmp_2646;
   std::complex<double> tmp_2647;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2647 += A0(MSe(gI1))*CpSRdpconjSRdpconjSeSe(gI1,gI1);
   }
   tmp_2646 += tmp_2647;
   result += (-1) * tmp_2646;
   std::complex<double> tmp_2648;
   std::complex<double> tmp_2649;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2649 += A0(MSu(gI1))*CpSRdpconjSRdpconjSuSu(gI1,gI1);
   }
   tmp_2648 += tmp_2649;
   result += (-3) * tmp_2648;
   std::complex<double> tmp_2650;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2651;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2651 += AbsSqr(CpconjSRdpconjSeSv(gI1,gI2))*B0(p,MSe(gI1),
            MSv(gI2));
      }
      tmp_2650 += tmp_2651;
   }
   result += tmp_2650;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2654;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2654 += AbsSqr(CpconjSRdpconjSdSu(gI1,gI2))*B0(p,MSd(gI1),
            MSu(gI2));
      }
      tmp_2653 += tmp_2654;
   }
   tmp_2652 += tmp_2653;
   result += (3) * tmp_2652;
   std::complex<double> tmp_2655;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2655 += AbsSqr(CpconjSRdpconjVWmRh(gI2))*F0(p,MRh(gI2),MVWm);
   }
   result += tmp_2655;
   std::complex<double> tmp_2656;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2656 += AbsSqr(CpconjSRdpSRdpAh(gI2))*B0(p,MSRdp,MAh(gI2));
   }
   result += tmp_2656;
   std::complex<double> tmp_2657;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2657 += AbsSqr(CpconjSRdpSRdphh(gI2))*B0(p,MSRdp,Mhh(gI2));
   }
   result += tmp_2657;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_SRum(double p ) const
{
   std::complex<double> result;

   result += -(A0(MSRdp)*CpSRumconjSRumconjSRdpSRdp());
   result += -(A0(MSRum)*CpSRumconjSRumconjSRumSRum());
   result += 4*A0(MVWm)*CpSRumconjSRumconjVWmVWm();
   result += 2*A0(MVZ)*CpSRumconjSRumVZVZ();
   result += AbsSqr(CpconjSRumVPSRum())*F0(p,MSRum,0);
   result += AbsSqr(CpconjSRumVZSRum())*F0(p,MSRum,MVZ);
   std::complex<double> tmp_2658;
   std::complex<double> tmp_2659;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2659 += A0(MRh(gI1))*CpSRumconjSRumconjRhRh(gI1,gI1);
   }
   tmp_2658 += tmp_2659;
   result += (-1) * tmp_2658;
   std::complex<double> tmp_2660;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2661;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2661 += AbsSqr(CpconjSRumRhHpm(gI1,gI2))*B0(p,MRh(gI1),MHpm(
            gI2));
      }
      tmp_2660 += tmp_2661;
   }
   result += tmp_2660;
   std::complex<double> tmp_2662;
   std::complex<double> tmp_2663;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2663 += A0(MSv(gI1))*CpSRumconjSRumconjSvSv(gI1,gI1);
   }
   tmp_2662 += tmp_2663;
   result += (-1) * tmp_2662;
   std::complex<double> tmp_2664;
   std::complex<double> tmp_2665;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2665 += A0(MAh(gI1))*CpSRumconjSRumAhAh(gI1,gI1);
   }
   tmp_2664 += tmp_2665;
   result += (-0.5) * tmp_2664;
   std::complex<double> tmp_2666;
   std::complex<double> tmp_2667;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2667 += A0(MHpm(gI1))*CpSRumconjSRumconjHpmHpm(gI1,gI1);
   }
   tmp_2666 += tmp_2667;
   result += (-1) * tmp_2666;
   std::complex<double> tmp_2668;
   std::complex<double> tmp_2669;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2669 += A0(Mhh(gI1))*CpSRumconjSRumhhhh(gI1,gI1);
   }
   tmp_2668 += tmp_2669;
   result += (-0.5) * tmp_2668;
   std::complex<double> tmp_2670;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2671;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2671 += (AbsSqr(CpconjSRumChiCha2PL(gI1,gI2)) + AbsSqr(
            CpconjSRumChiCha2PR(gI1,gI2)))*G0(p,MChi(gI1),MCha2(gI2));
      }
      tmp_2670 += tmp_2671;
   }
   result += tmp_2670;
   std::complex<double> tmp_2672;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2673;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2673 += AbsSqr(CpconjSRumHpmAh(gI1,gI2))*B0(p,MHpm(gI1),MAh(
            gI2));
      }
      tmp_2672 += tmp_2673;
   }
   result += tmp_2672;
   std::complex<double> tmp_2674;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2675;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2675 += AbsSqr(CpconjSRumHpmhh(gI1,gI2))*B0(p,MHpm(gI1),Mhh(
            gI2));
      }
      tmp_2674 += tmp_2675;
   }
   result += tmp_2674;
   std::complex<double> tmp_2676;
   std::complex<double> tmp_2677;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2678;
      std::complex<double> tmp_2679;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2679 += B0(p,MChi(gI1),MCha2(gI2))*(Conj(CpconjSRumChiCha2PR
            (gI1,gI2))*CpconjSRumChiCha2PL(gI1,gI2) + Conj(CpconjSRumChiCha2PL(gI1
            ,gI2))*CpconjSRumChiCha2PR(gI1,gI2))*MCha2(gI2);
      }
      tmp_2678 += tmp_2679;
      tmp_2677 += (MChi(gI1)) * tmp_2678;
   }
   tmp_2676 += tmp_2677;
   result += (-2) * tmp_2676;
   std::complex<double> tmp_2680;
   std::complex<double> tmp_2681;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2681 += A0(MSd(gI1))*CpSRumconjSRumconjSdSd(gI1,gI1);
   }
   tmp_2680 += tmp_2681;
   result += (-3) * tmp_2680;
   std::complex<double> tmp_2682;
   std::complex<double> tmp_2683;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2683 += A0(MSe(gI1))*CpSRumconjSRumconjSeSe(gI1,gI1);
   }
   tmp_2682 += tmp_2683;
   result += (-1) * tmp_2682;
   std::complex<double> tmp_2684;
   std::complex<double> tmp_2685;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2685 += A0(MSu(gI1))*CpSRumconjSRumconjSuSu(gI1,gI1);
   }
   tmp_2684 += tmp_2685;
   result += (-3) * tmp_2684;
   std::complex<double> tmp_2686;
   std::complex<double> tmp_2687;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2688;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2688 += AbsSqr(CpconjSRumconjSuSd(gI1,gI2))*B0(p,MSu(gI1),
            MSd(gI2));
      }
      tmp_2687 += tmp_2688;
   }
   tmp_2686 += tmp_2687;
   result += (3) * tmp_2686;
   std::complex<double> tmp_2689;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2689 += AbsSqr(CpconjSRumVWmRh(gI2))*F0(p,MRh(gI2),MVWm);
   }
   result += tmp_2689;
   std::complex<double> tmp_2690;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2690 += AbsSqr(CpconjSRumSRumAh(gI2))*B0(p,MSRum,MAh(gI2));
   }
   result += tmp_2690;
   std::complex<double> tmp_2691;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2691 += AbsSqr(CpconjSRumSRumhh(gI2))*B0(p,MSRum,Mhh(gI2));
   }
   result += tmp_2691;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_sigmaO(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MphiO)*CpsigmaOsigmaOphiOphiO();
   result += 3*AbsSqr(CpsigmaOVGsigmaO())*F0(p,MsigmaO,0);
   std::complex<double> tmp_2692;
   std::complex<double> tmp_2693;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2694;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2694 += AbsSqr(CpsigmaOconjSdSd(gI1,gI2))*B0(p,MSd(gI1),MSd(
            gI2));
      }
      tmp_2693 += tmp_2694;
   }
   tmp_2692 += tmp_2693;
   result += (0.5) * tmp_2692;
   std::complex<double> tmp_2695;
   std::complex<double> tmp_2696;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2697;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2697 += AbsSqr(CpsigmaOconjSuSu(gI1,gI2))*B0(p,MSu(gI1),MSu(
            gI2));
      }
      tmp_2696 += tmp_2697;
   }
   tmp_2695 += tmp_2696;
   result += (0.5) * tmp_2695;
   result += 3*(AbsSqr(CpsigmaObarGluGluPL(0,0)) + AbsSqr(CpsigmaObarGluGluPR(0
      ,0)))*G0(p,MGlu,MGlu);
   result += -6*B0(p,MGlu,MGlu)*(Conj(CpsigmaObarGluGluPR(0,0))*
      CpsigmaObarGluGluPL(0,0) + Conj(CpsigmaObarGluGluPL(0,0))*
      CpsigmaObarGluGluPR(0,0))*Sqr(MGlu);

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_phiO(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MsigmaO)*CpphiOphiOsigmaOsigmaO();
   result += 3*AbsSqr(CpphiOVGphiO())*F0(p,MphiO,0);
   std::complex<double> tmp_2698;
   std::complex<double> tmp_2699;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2700;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2700 += AbsSqr(CpphiOconjSdSd(gI1,gI2))*B0(p,MSd(gI1),MSd(
            gI2));
      }
      tmp_2699 += tmp_2700;
   }
   tmp_2698 += tmp_2699;
   result += (0.5) * tmp_2698;
   std::complex<double> tmp_2701;
   std::complex<double> tmp_2702;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2703;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2703 += AbsSqr(CpphiOconjSuSu(gI1,gI2))*B0(p,MSu(gI1),MSu(
            gI2));
      }
      tmp_2702 += tmp_2703;
   }
   tmp_2701 += tmp_2702;
   result += (0.5) * tmp_2701;
   result += 3*(AbsSqr(CpphiObarGluGluPL(0,0)) + AbsSqr(CpphiObarGluGluPR(0,0))
      )*G0(p,MGlu,MGlu);
   result += -6*B0(p,MGlu,MGlu)*(Conj(CpphiObarGluGluPR(0,0))*CpphiObarGluGluPL
      (0,0) + Conj(CpphiObarGluGluPL(0,0))*CpphiObarGluGluPR(0,0))*Sqr(MGlu);

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VG(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpVGbargGgG())*B00(p,MVG,MVG);
   result += -6*AbsSqr(CpVGphiOphiO())*B00(p,MphiO,MphiO);
   result += -6*AbsSqr(CpVGsigmaOsigmaO())*B00(p,MsigmaO,MsigmaO);
   result += 499.5*A0(MphiO)*CpVGVGphiOphiO();
   result += 499.5*A0(MsigmaO)*CpVGVGsigmaOsigmaO();
   result += -1.5*AbsSqr(CpVGVGVG())*(10*B00(p,0,0) + 4*B0(p,0,0)*Sqr(p));
   result += 0;
   std::complex<double> tmp_2704;
   std::complex<double> tmp_2705;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2706;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2706 += (AbsSqr(CpVGbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVGbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2706 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVGbarFdFdPL(gI1,gI2))*CpVGbarFdFdPR(gI1,gI2));
      }
      tmp_2705 += tmp_2706;
   }
   tmp_2704 += tmp_2705;
   result += (0.5) * tmp_2704;
   std::complex<double> tmp_2707;
   std::complex<double> tmp_2708;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2709;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2709 += (AbsSqr(CpVGbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVGbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2709 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVGbarFuFuPL(gI1,gI2))*CpVGbarFuFuPR(gI1,gI2));
      }
      tmp_2708 += tmp_2709;
   }
   tmp_2707 += tmp_2708;
   result += (0.5) * tmp_2707;
   std::complex<double> tmp_2710;
   std::complex<double> tmp_2711;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2711 += A0(MSd(gI1))*CpVGVGconjSdSd(gI1,gI1);
   }
   tmp_2710 += tmp_2711;
   result += (999) * tmp_2710;
   std::complex<double> tmp_2712;
   std::complex<double> tmp_2713;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2713 += A0(MSu(gI1))*CpVGVGconjSuSu(gI1,gI1);
   }
   tmp_2712 += tmp_2713;
   result += (999) * tmp_2712;
   std::complex<double> tmp_2714;
   std::complex<double> tmp_2715;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2716;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2716 += AbsSqr(CpVGconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2715 += tmp_2716;
   }
   tmp_2714 += tmp_2715;
   result += (-2) * tmp_2714;
   std::complex<double> tmp_2717;
   std::complex<double> tmp_2718;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2719;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2719 += AbsSqr(CpVGconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2718 += tmp_2719;
   }
   tmp_2717 += tmp_2718;
   result += (-2) * tmp_2717;
   result += 3*((AbsSqr(CpVGbarGluGluPL(0,0)) + AbsSqr(CpVGbarGluGluPR(0,0)))*
      H0(p,MGlu,MGlu) + 4*B0(p,MGlu,MGlu)*Re(Conj(CpVGbarGluGluPL(0,0))*
      CpVGbarGluGluPR(0,0))*Sqr(MGlu));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VP(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVPbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVPbargWmgWm())*B00(p,MVWm,MVWm);
   result += -4*AbsSqr(CpVPconjSRdpSRdp())*B00(p,MSRdp,MSRdp);
   result += -4*AbsSqr(CpVPconjSRumSRum())*B00(p,MSRum,MSRum);
   result += A0(MSRdp)*CpVPVPconjSRdpSRdp();
   result += A0(MSRum)*CpVPVPconjSRumSRum();
   result += -(A0(MVWm)*(4*CpVPVPconjVWmVWm1() + CpVPVPconjVWmVWm2() +
      CpVPVPconjVWmVWm3()));
   std::complex<double> tmp_2720;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2721;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2721 += (AbsSqr(CpVPbarCha1Cha1PL(gI1,gI2)) + AbsSqr(
            CpVPbarCha1Cha1PR(gI1,gI2)))*H0(p,MCha1(gI1),MCha1(gI2));
         tmp_2721 += 4*B0(p,MCha1(gI1),MCha1(gI2))*MCha1(gI1)*MCha1(gI2)*
            Re(Conj(CpVPbarCha1Cha1PL(gI1,gI2))*CpVPbarCha1Cha1PR(gI1,gI2));
      }
      tmp_2720 += tmp_2721;
   }
   result += tmp_2720;
   std::complex<double> tmp_2722;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2723;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2723 += (AbsSqr(CpVPbarCha2Cha2PL(gI1,gI2)) + AbsSqr(
            CpVPbarCha2Cha2PR(gI1,gI2)))*H0(p,MCha2(gI1),MCha2(gI2));
         tmp_2723 += 4*B0(p,MCha2(gI1),MCha2(gI2))*MCha2(gI1)*MCha2(gI2)*
            Re(Conj(CpVPbarCha2Cha2PL(gI1,gI2))*CpVPbarCha2Cha2PR(gI1,gI2));
      }
      tmp_2722 += tmp_2723;
   }
   result += tmp_2722;
   std::complex<double> tmp_2724;
   std::complex<double> tmp_2725;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2726;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2726 += (AbsSqr(CpVPbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVPbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2726 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVPbarFdFdPL(gI1,gI2))*CpVPbarFdFdPR(gI1,gI2));
      }
      tmp_2725 += tmp_2726;
   }
   tmp_2724 += tmp_2725;
   result += (3) * tmp_2724;
   std::complex<double> tmp_2727;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2728;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2728 += (AbsSqr(CpVPbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVPbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2728 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVPbarFeFePL(gI1,gI2))*CpVPbarFeFePR(gI1,gI2));
      }
      tmp_2727 += tmp_2728;
   }
   result += tmp_2727;
   std::complex<double> tmp_2729;
   std::complex<double> tmp_2730;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2731;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2731 += (AbsSqr(CpVPbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVPbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2731 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVPbarFuFuPL(gI1,gI2))*CpVPbarFuFuPR(gI1,gI2));
      }
      tmp_2730 += tmp_2731;
   }
   tmp_2729 += tmp_2730;
   result += (3) * tmp_2729;
   std::complex<double> tmp_2732;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2732 += A0(MHpm(gI1))*CpVPVPconjHpmHpm(gI1,gI1);
   }
   result += tmp_2732;
   std::complex<double> tmp_2733;
   std::complex<double> tmp_2734;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2735;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2735 += AbsSqr(CpVPconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2734 += tmp_2735;
   }
   tmp_2733 += tmp_2734;
   result += (-4) * tmp_2733;
   std::complex<double> tmp_2736;
   std::complex<double> tmp_2737;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2737 += A0(MSd(gI1))*CpVPVPconjSdSd(gI1,gI1);
   }
   tmp_2736 += tmp_2737;
   result += (3) * tmp_2736;
   std::complex<double> tmp_2738;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2738 += A0(MSe(gI1))*CpVPVPconjSeSe(gI1,gI1);
   }
   result += tmp_2738;
   std::complex<double> tmp_2739;
   std::complex<double> tmp_2740;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2740 += A0(MSu(gI1))*CpVPVPconjSuSu(gI1,gI1);
   }
   tmp_2739 += tmp_2740;
   result += (3) * tmp_2739;
   std::complex<double> tmp_2741;
   std::complex<double> tmp_2742;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2743;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2743 += AbsSqr(CpVPconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2742 += tmp_2743;
   }
   tmp_2741 += tmp_2742;
   result += (-12) * tmp_2741;
   std::complex<double> tmp_2744;
   std::complex<double> tmp_2745;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2746;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2746 += AbsSqr(CpVPconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2745 += tmp_2746;
   }
   tmp_2744 += tmp_2745;
   result += (-4) * tmp_2744;
   std::complex<double> tmp_2747;
   std::complex<double> tmp_2748;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2749;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2749 += AbsSqr(CpVPconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2748 += tmp_2749;
   }
   tmp_2747 += tmp_2748;
   result += (-12) * tmp_2747;
   std::complex<double> tmp_2750;
   std::complex<double> tmp_2751;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2751 += AbsSqr(CpVPconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2750 += tmp_2751;
   result += (2) * tmp_2750;
   result += -(AbsSqr(CpVPconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0(
      p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -4*AbsSqr(CpVZconjSRdpSRdp())*B00(p,MSRdp,MSRdp);
   result += -4*AbsSqr(CpVZconjSRumSRum())*B00(p,MSRum,MSRum);
   result += A0(MSRdp)*CpVZVZconjSRdpSRdp();
   result += A0(MSRum)*CpVZVZconjSRumSRum();
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2752;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2752 += A0(MRh(gI1))*CpVZVZconjRhRh(gI1,gI1);
   }
   result += tmp_2752;
   std::complex<double> tmp_2753;
   std::complex<double> tmp_2754;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2755;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2755 += AbsSqr(CpVZconjRhRh(gI1,gI2))*B00(p,MRh(gI1),MRh(gI2
            ));
      }
      tmp_2754 += tmp_2755;
   }
   tmp_2753 += tmp_2754;
   result += (-4) * tmp_2753;
   std::complex<double> tmp_2756;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2757;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2757 += (AbsSqr(CpVZbarCha1Cha1PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha1Cha1PR(gI1,gI2)))*H0(p,MCha1(gI1),MCha1(gI2));
         tmp_2757 += 4*B0(p,MCha1(gI1),MCha1(gI2))*MCha1(gI1)*MCha1(gI2)*
            Re(Conj(CpVZbarCha1Cha1PL(gI1,gI2))*CpVZbarCha1Cha1PR(gI1,gI2));
      }
      tmp_2756 += tmp_2757;
   }
   result += tmp_2756;
   std::complex<double> tmp_2758;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2759;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2759 += (AbsSqr(CpVZbarCha2Cha2PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha2Cha2PR(gI1,gI2)))*H0(p,MCha2(gI1),MCha2(gI2));
         tmp_2759 += 4*B0(p,MCha2(gI1),MCha2(gI2))*MCha2(gI1)*MCha2(gI2)*
            Re(Conj(CpVZbarCha2Cha2PL(gI1,gI2))*CpVZbarCha2Cha2PR(gI1,gI2));
      }
      tmp_2758 += tmp_2759;
   }
   result += tmp_2758;
   std::complex<double> tmp_2760;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2760 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2760;
   std::complex<double> tmp_2761;
   std::complex<double> tmp_2762;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2763;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2763 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2762 += tmp_2763;
   }
   tmp_2761 += tmp_2762;
   result += (-4) * tmp_2761;
   std::complex<double> tmp_2764;
   std::complex<double> tmp_2765;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2766;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2766 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2766 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2765 += tmp_2766;
   }
   tmp_2764 += tmp_2765;
   result += (3) * tmp_2764;
   std::complex<double> tmp_2767;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2768;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2768 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2768 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2767 += tmp_2768;
   }
   result += tmp_2767;
   std::complex<double> tmp_2769;
   std::complex<double> tmp_2770;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2771;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2771 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2771 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2770 += tmp_2771;
   }
   tmp_2769 += tmp_2770;
   result += (3) * tmp_2769;
   std::complex<double> tmp_2772;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2773;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2773 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2773 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_2772 += tmp_2773;
   }
   result += tmp_2772;
   std::complex<double> tmp_2774;
   std::complex<double> tmp_2775;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2775 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2774 += tmp_2775;
   result += (0.5) * tmp_2774;
   std::complex<double> tmp_2776;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2776 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2776;
   std::complex<double> tmp_2777;
   std::complex<double> tmp_2778;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2778 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2777 += tmp_2778;
   result += (0.5) * tmp_2777;
   std::complex<double> tmp_2779;
   std::complex<double> tmp_2780;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2781;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2781 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2780 += tmp_2781;
   }
   tmp_2779 += tmp_2780;
   result += (-4) * tmp_2779;
   std::complex<double> tmp_2782;
   std::complex<double> tmp_2783;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2784;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2784 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2783 += tmp_2784;
   }
   tmp_2782 += tmp_2783;
   result += (-4) * tmp_2782;
   std::complex<double> tmp_2785;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2786;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2786 += (AbsSqr(CpVZbarChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZbarChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2786 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZbarChiChiPL(gI1,gI2))*CpVZbarChiChiPR(gI1,gI2));
      }
      tmp_2785 += tmp_2786;
   }
   result += tmp_2785;
   std::complex<double> tmp_2787;
   std::complex<double> tmp_2788;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2788 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2787 += tmp_2788;
   result += (3) * tmp_2787;
   std::complex<double> tmp_2789;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2789 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2789;
   std::complex<double> tmp_2790;
   std::complex<double> tmp_2791;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2791 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2790 += tmp_2791;
   result += (3) * tmp_2790;
   std::complex<double> tmp_2792;
   std::complex<double> tmp_2793;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2794;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2794 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2793 += tmp_2794;
   }
   tmp_2792 += tmp_2793;
   result += (-12) * tmp_2792;
   std::complex<double> tmp_2795;
   std::complex<double> tmp_2796;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2797;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2797 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2796 += tmp_2797;
   }
   tmp_2795 += tmp_2796;
   result += (-4) * tmp_2795;
   std::complex<double> tmp_2798;
   std::complex<double> tmp_2799;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2800;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2800 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2799 += tmp_2800;
   }
   tmp_2798 += tmp_2799;
   result += (-12) * tmp_2798;
   std::complex<double> tmp_2801;
   std::complex<double> tmp_2802;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2802 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2801 += tmp_2802;
   result += (2) * tmp_2801;
   std::complex<double> tmp_2803;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2803 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2803;
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
   result += A0(MSRdp)*CpVWmconjVWmconjSRdpSRdp();
   result += A0(MSRum)*CpVWmconjVWmconjSRumSRum();
   result += -(A0(MVWm)*(4*CpVWmconjVWmconjVWmVWm1() + CpVWmconjVWmconjVWmVWm2(
      ) + CpVWmconjVWmconjVWmVWm3()));
   result += 0;
   result += -0.5*A0(MVZ)*(4*CpVWmconjVWmVZVZ1() + CpVWmconjVWmVZVZ2() +
      CpVWmconjVWmVZVZ3());
   std::complex<double> tmp_2804;
   std::complex<double> tmp_2805;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2805 += AbsSqr(CpconjVWmconjRhSRum(gI1))*B00(p,MSRum,MRh(gI1));
   }
   tmp_2804 += tmp_2805;
   result += (-4) * tmp_2804;
   std::complex<double> tmp_2806;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2806 += A0(MRh(gI1))*CpVWmconjVWmconjRhRh(gI1,gI1);
   }
   result += tmp_2806;
   std::complex<double> tmp_2807;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2808;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2808 += (AbsSqr(CpconjVWmbarCha1ChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarCha1ChiPR(gI1,gI2)))*H0(p,MCha1(gI1),MChi(gI2));
         tmp_2808 += 4*B0(p,MCha1(gI1),MChi(gI2))*MCha1(gI1)*MChi(gI2)*Re
            (Conj(CpconjVWmbarCha1ChiPL(gI1,gI2))*CpconjVWmbarCha1ChiPR(gI1,gI2));
      }
      tmp_2807 += tmp_2808;
   }
   result += tmp_2807;
   std::complex<double> tmp_2809;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2809 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2809;
   std::complex<double> tmp_2810;
   std::complex<double> tmp_2811;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2812;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2812 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_2812 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_2811 += tmp_2812;
   }
   tmp_2810 += tmp_2811;
   result += (3) * tmp_2810;
   std::complex<double> tmp_2813;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2814;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2814 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_2814 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_2813 += tmp_2814;
   }
   result += tmp_2813;
   std::complex<double> tmp_2815;
   std::complex<double> tmp_2816;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2817;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2817 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2816 += tmp_2817;
   }
   tmp_2815 += tmp_2816;
   result += (-4) * tmp_2815;
   std::complex<double> tmp_2818;
   std::complex<double> tmp_2819;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2819 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_2818 += tmp_2819;
   result += (0.5) * tmp_2818;
   std::complex<double> tmp_2820;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2820 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_2820;
   std::complex<double> tmp_2821;
   std::complex<double> tmp_2822;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2822 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2821 += tmp_2822;
   result += (0.5) * tmp_2821;
   std::complex<double> tmp_2823;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2824;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2824 += (AbsSqr(CpconjVWmbarChiCha2PL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarChiCha2PR(gI1,gI2)))*H0(p,MChi(gI1),MCha2(gI2));
         tmp_2824 += 4*B0(p,MChi(gI1),MCha2(gI2))*MCha2(gI2)*MChi(gI1)*Re
            (Conj(CpconjVWmbarChiCha2PL(gI1,gI2))*CpconjVWmbarChiCha2PR(gI1,gI2));
      }
      tmp_2823 += tmp_2824;
   }
   result += tmp_2823;
   std::complex<double> tmp_2825;
   std::complex<double> tmp_2826;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2827;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2827 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_2826 += tmp_2827;
   }
   tmp_2825 += tmp_2826;
   result += (-4) * tmp_2825;
   std::complex<double> tmp_2828;
   std::complex<double> tmp_2829;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2830;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2830 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_2829 += tmp_2830;
   }
   tmp_2828 += tmp_2829;
   result += (-4) * tmp_2828;
   std::complex<double> tmp_2831;
   std::complex<double> tmp_2832;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2832 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2831 += tmp_2832;
   result += (3) * tmp_2831;
   std::complex<double> tmp_2833;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2833 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2833;
   std::complex<double> tmp_2834;
   std::complex<double> tmp_2835;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2835 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2834 += tmp_2835;
   result += (3) * tmp_2834;
   std::complex<double> tmp_2836;
   std::complex<double> tmp_2837;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2838;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2838 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2837 += tmp_2838;
   }
   tmp_2836 += tmp_2837;
   result += (-12) * tmp_2836;
   std::complex<double> tmp_2839;
   std::complex<double> tmp_2840;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2840 += AbsSqr(CpconjVWmconjSRdpRh(gI2))*B00(p,MRh(gI2),MSRdp);
   }
   tmp_2839 += tmp_2840;
   result += (-4) * tmp_2839;
   std::complex<double> tmp_2841;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2841 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_2841;
   std::complex<double> tmp_2842;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2842 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2842;
   std::complex<double> tmp_2843;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2843 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_2843;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2844;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2844 += B0(p,MCha1(gI1),MSRdp)*Conj(CpbarUChibarCha1SRdpPL(gO2,gI1
         ))*CpbarUChibarCha1SRdpPR(gO1,gI1)*MCha1(gI1);
   }
   result += tmp_2844;
   std::complex<double> tmp_2845;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2845 += B0(p,MCha2(gI1),MSRum)*Conj(CpbarUChibarCha2SRumPL(gO2,gI1
         ))*CpbarUChibarCha2SRumPR(gO1,gI1)*MCha2(gI1);
   }
   result += tmp_2845;
   std::complex<double> tmp_2846;
   std::complex<double> tmp_2847;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2848;
      std::complex<double> tmp_2849;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2849 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPL(gO2,
            gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2848 += tmp_2849;
      tmp_2847 += (MFd(gI1)) * tmp_2848;
   }
   tmp_2846 += tmp_2847;
   result += (3) * tmp_2846;
   std::complex<double> tmp_2850;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2851;
      std::complex<double> tmp_2852;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2852 += B0(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePL(gO2,
            gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2851 += tmp_2852;
      tmp_2850 += (MFe(gI1)) * tmp_2851;
   }
   result += tmp_2850;
   std::complex<double> tmp_2853;
   std::complex<double> tmp_2854;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2855;
      std::complex<double> tmp_2856;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2856 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPL(gO2,
            gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2855 += tmp_2856;
      tmp_2854 += (MFu(gI1)) * tmp_2855;
   }
   tmp_2853 += tmp_2854;
   result += (3) * tmp_2853;
   std::complex<double> tmp_2857;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2858;
      std::complex<double> tmp_2859;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2859 += B0(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPL(gO2,
            gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2858 += tmp_2859;
      tmp_2857 += (MFv(gI1)) * tmp_2858;
   }
   result += tmp_2857;
   std::complex<double> tmp_2860;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2861;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2861 += B0(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PL(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_2860 += tmp_2861;
   }
   result += tmp_2860;
   std::complex<double> tmp_2862;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2863;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2863 += B0(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PL(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PR(gO1,gI1,
            gI2)*MCha2(gI2);
      }
      tmp_2862 += tmp_2863;
   }
   result += tmp_2862;
   std::complex<double> tmp_2864;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2865;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2865 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPL(gO2,
            gI1,gI2))*CpbarUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2864 += tmp_2865;
   }
   result += tmp_2864;
   std::complex<double> tmp_2866;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2867;
      std::complex<double> tmp_2868;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2868 += B0(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPL(
            gO2,gI1,gI2))*CpbarUChibarChiRhPR(gO1,gI1,gI2);
      }
      tmp_2867 += tmp_2868;
      tmp_2866 += (MChi(gI1)) * tmp_2867;
   }
   result += tmp_2866;
   std::complex<double> tmp_2869;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2870;
      std::complex<double> tmp_2871;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2871 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPL(gO2,
            gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2870 += tmp_2871;
      tmp_2869 += (MChi(gI1)) * tmp_2870;
   }
   result += tmp_2869;
   std::complex<double> tmp_2872;
   std::complex<double> tmp_2873;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2874;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2874 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPL(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2873 += tmp_2874;
   }
   tmp_2872 += tmp_2873;
   result += (3) * tmp_2872;
   std::complex<double> tmp_2875;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2876;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2876 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePL(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2875 += tmp_2876;
   }
   result += tmp_2875;
   std::complex<double> tmp_2877;
   std::complex<double> tmp_2878;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2879;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2879 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPL(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2878 += tmp_2879;
   }
   tmp_2877 += tmp_2878;
   result += (3) * tmp_2877;
   std::complex<double> tmp_2880;
   std::complex<double> tmp_2881;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2881 += B0(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PR(gO2,gI2))*
         CpbarUChiVWmCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_2880 += tmp_2881;
   result += (-4) * tmp_2880;
   std::complex<double> tmp_2882;
   std::complex<double> tmp_2883;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2883 += B0(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PR(gO2,gI2)
         )*CpbarUChiconjVWmCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_2882 += tmp_2883;
   result += (-4) * tmp_2882;
   std::complex<double> tmp_2884;
   std::complex<double> tmp_2885;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2885 += B0(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPR(gO2,gI2))*
         CpbarUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2884 += tmp_2885;
   result += (-4) * tmp_2884;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2886;
   std::complex<double> tmp_2887;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2887 += B1(p,MCha1(gI1),MSRdp)*Conj(CpbarUChibarCha1SRdpPR(gO2,gI1
         ))*CpbarUChibarCha1SRdpPR(gO1,gI1);
   }
   tmp_2886 += tmp_2887;
   result += (-0.5) * tmp_2886;
   std::complex<double> tmp_2888;
   std::complex<double> tmp_2889;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2889 += B1(p,MCha2(gI1),MSRum)*Conj(CpbarUChibarCha2SRumPR(gO2,gI1
         ))*CpbarUChibarCha2SRumPR(gO1,gI1);
   }
   tmp_2888 += tmp_2889;
   result += (-0.5) * tmp_2888;
   std::complex<double> tmp_2890;
   std::complex<double> tmp_2891;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2892;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2892 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPR(gO2,
            gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2891 += tmp_2892;
   }
   tmp_2890 += tmp_2891;
   result += (-0.5) * tmp_2890;
   std::complex<double> tmp_2893;
   std::complex<double> tmp_2894;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2895;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2895 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPR(gO2,
            gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2894 += tmp_2895;
   }
   tmp_2893 += tmp_2894;
   result += (-1.5) * tmp_2893;
   std::complex<double> tmp_2896;
   std::complex<double> tmp_2897;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2898;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2898 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePR(gO2,
            gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2897 += tmp_2898;
   }
   tmp_2896 += tmp_2897;
   result += (-0.5) * tmp_2896;
   std::complex<double> tmp_2899;
   std::complex<double> tmp_2900;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2901;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2901 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPR(gO2,
            gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2900 += tmp_2901;
   }
   tmp_2899 += tmp_2900;
   result += (-1.5) * tmp_2899;
   std::complex<double> tmp_2902;
   std::complex<double> tmp_2903;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2904;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2904 += B1(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPR(
            gO2,gI1,gI2))*CpbarUChibarChiRhPR(gO1,gI1,gI2);
      }
      tmp_2903 += tmp_2904;
   }
   tmp_2902 += tmp_2903;
   result += (-0.5) * tmp_2902;
   std::complex<double> tmp_2905;
   std::complex<double> tmp_2906;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2907;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2907 += B1(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PR(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PR(gO1,gI1,
            gI2);
      }
      tmp_2906 += tmp_2907;
   }
   tmp_2905 += tmp_2906;
   result += (-0.5) * tmp_2905;
   std::complex<double> tmp_2908;
   std::complex<double> tmp_2909;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2910;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2910 += B1(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PR(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PR(gO1,gI1,gI2);
      }
      tmp_2909 += tmp_2910;
   }
   tmp_2908 += tmp_2909;
   result += (-0.5) * tmp_2908;
   std::complex<double> tmp_2911;
   std::complex<double> tmp_2912;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2913;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2913 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPR(gO2,
            gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2912 += tmp_2913;
   }
   tmp_2911 += tmp_2912;
   result += (-0.5) * tmp_2911;
   std::complex<double> tmp_2914;
   std::complex<double> tmp_2915;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2916;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2916 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPR(gO2,
            gI1,gI2))*CpbarUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_2915 += tmp_2916;
   }
   tmp_2914 += tmp_2915;
   result += (-0.5) * tmp_2914;
   std::complex<double> tmp_2917;
   std::complex<double> tmp_2918;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2919;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2919 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPR(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_2918 += tmp_2919;
   }
   tmp_2917 += tmp_2918;
   result += (-1.5) * tmp_2917;
   std::complex<double> tmp_2920;
   std::complex<double> tmp_2921;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2922;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2922 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePR(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_2921 += tmp_2922;
   }
   tmp_2920 += tmp_2921;
   result += (-0.5) * tmp_2920;
   std::complex<double> tmp_2923;
   std::complex<double> tmp_2924;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2925;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2925 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPR(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_2924 += tmp_2925;
   }
   tmp_2923 += tmp_2924;
   result += (-1.5) * tmp_2923;
   std::complex<double> tmp_2926;
   std::complex<double> tmp_2927;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2927 += B1(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PL(gO2,gI2)
         )*CpbarUChiconjVWmCha2PL(gO1,gI2);
   }
   tmp_2926 += tmp_2927;
   result += (-1) * tmp_2926;
   std::complex<double> tmp_2928;
   std::complex<double> tmp_2929;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2929 += B1(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PL(gO2,gI2))*
         CpbarUChiVWmCha1PL(gO1,gI2);
   }
   tmp_2928 += tmp_2929;
   result += (-1) * tmp_2928;
   std::complex<double> tmp_2930;
   std::complex<double> tmp_2931;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2931 += B1(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPL(gO2,gI2))*
         CpbarUChiVZChiPL(gO1,gI2);
   }
   tmp_2930 += tmp_2931;
   result += (-1) * tmp_2930;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2932;
   std::complex<double> tmp_2933;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2933 += B1(p,MCha1(gI1),MSRdp)*Conj(CpbarUChibarCha1SRdpPL(gO2,gI1
         ))*CpbarUChibarCha1SRdpPL(gO1,gI1);
   }
   tmp_2932 += tmp_2933;
   result += (-0.5) * tmp_2932;
   std::complex<double> tmp_2934;
   std::complex<double> tmp_2935;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2935 += B1(p,MCha2(gI1),MSRum)*Conj(CpbarUChibarCha2SRumPL(gO2,gI1
         ))*CpbarUChibarCha2SRumPL(gO1,gI1);
   }
   tmp_2934 += tmp_2935;
   result += (-0.5) * tmp_2934;
   std::complex<double> tmp_2936;
   std::complex<double> tmp_2937;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2938;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2938 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPL(gO2,
            gI1,gI2))*CpbarUChibarFvSvPL(gO1,gI1,gI2);
      }
      tmp_2937 += tmp_2938;
   }
   tmp_2936 += tmp_2937;
   result += (-0.5) * tmp_2936;
   std::complex<double> tmp_2939;
   std::complex<double> tmp_2940;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2941;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2941 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPL(gO2,
            gI1,gI2))*CpbarUChibarFdSdPL(gO1,gI1,gI2);
      }
      tmp_2940 += tmp_2941;
   }
   tmp_2939 += tmp_2940;
   result += (-1.5) * tmp_2939;
   std::complex<double> tmp_2942;
   std::complex<double> tmp_2943;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2944;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2944 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePL(gO2,
            gI1,gI2))*CpbarUChibarFeSePL(gO1,gI1,gI2);
      }
      tmp_2943 += tmp_2944;
   }
   tmp_2942 += tmp_2943;
   result += (-0.5) * tmp_2942;
   std::complex<double> tmp_2945;
   std::complex<double> tmp_2946;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2947;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2947 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPL(gO2,
            gI1,gI2))*CpbarUChibarFuSuPL(gO1,gI1,gI2);
      }
      tmp_2946 += tmp_2947;
   }
   tmp_2945 += tmp_2946;
   result += (-1.5) * tmp_2945;
   std::complex<double> tmp_2948;
   std::complex<double> tmp_2949;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2950;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2950 += B1(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPL(
            gO2,gI1,gI2))*CpbarUChibarChiRhPL(gO1,gI1,gI2);
      }
      tmp_2949 += tmp_2950;
   }
   tmp_2948 += tmp_2949;
   result += (-0.5) * tmp_2948;
   std::complex<double> tmp_2951;
   std::complex<double> tmp_2952;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2953;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2953 += B1(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PL(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PL(gO1,gI1,
            gI2);
      }
      tmp_2952 += tmp_2953;
   }
   tmp_2951 += tmp_2952;
   result += (-0.5) * tmp_2951;
   std::complex<double> tmp_2954;
   std::complex<double> tmp_2955;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2956;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2956 += B1(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PL(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PL(gO1,gI1,gI2);
      }
      tmp_2955 += tmp_2956;
   }
   tmp_2954 += tmp_2955;
   result += (-0.5) * tmp_2954;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2959;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2959 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPL(gO2,
            gI1,gI2))*CpbarUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_2958 += tmp_2959;
   }
   tmp_2957 += tmp_2958;
   result += (-0.5) * tmp_2957;
   std::complex<double> tmp_2960;
   std::complex<double> tmp_2961;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2962;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2962 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPL(gO2,
            gI1,gI2))*CpbarUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_2961 += tmp_2962;
   }
   tmp_2960 += tmp_2961;
   result += (-0.5) * tmp_2960;
   std::complex<double> tmp_2963;
   std::complex<double> tmp_2964;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2965;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2965 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPL(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_2964 += tmp_2965;
   }
   tmp_2963 += tmp_2964;
   result += (-1.5) * tmp_2963;
   std::complex<double> tmp_2966;
   std::complex<double> tmp_2967;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2968;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2968 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePL(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_2967 += tmp_2968;
   }
   tmp_2966 += tmp_2967;
   result += (-0.5) * tmp_2966;
   std::complex<double> tmp_2969;
   std::complex<double> tmp_2970;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2971;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2971 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPL(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_2970 += tmp_2971;
   }
   tmp_2969 += tmp_2970;
   result += (-1.5) * tmp_2969;
   std::complex<double> tmp_2972;
   std::complex<double> tmp_2973;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2973 += B1(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PR(gO2,gI2)
         )*CpbarUChiconjVWmCha2PR(gO1,gI2);
   }
   tmp_2972 += tmp_2973;
   result += (-1) * tmp_2972;
   std::complex<double> tmp_2974;
   std::complex<double> tmp_2975;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2975 += B1(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PR(gO2,gI2))*
         CpbarUChiVWmCha1PR(gO1,gI2);
   }
   tmp_2974 += tmp_2975;
   result += (-1) * tmp_2974;
   std::complex<double> tmp_2976;
   std::complex<double> tmp_2977;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2977 += B1(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPR(gO2,gI2))*
         CpbarUChiVZChiPR(gO1,gI2);
   }
   tmp_2976 += tmp_2977;
   result += (-1) * tmp_2976;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2978;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2979;
      std::complex<double> tmp_2980;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2980 += B0(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPL(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,gI2);
      }
      tmp_2979 += tmp_2980;
      tmp_2978 += (MCha1(gI1)) * tmp_2979;
   }
   result += tmp_2978;
   std::complex<double> tmp_2981;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2982;
      std::complex<double> tmp_2983;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2983 += B0(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPL
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(gO1,gI1,gI2);
      }
      tmp_2982 += tmp_2983;
      tmp_2981 += (MCha2(gI1)) * tmp_2982;
   }
   result += tmp_2981;
   std::complex<double> tmp_2984;
   std::complex<double> tmp_2985;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2986;
      std::complex<double> tmp_2987;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2987 += B0(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPL(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,gI2);
      }
      tmp_2986 += tmp_2987;
      tmp_2985 += (MFd(gI1)) * tmp_2986;
   }
   tmp_2984 += tmp_2985;
   result += (3) * tmp_2984;
   std::complex<double> tmp_2988;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2989;
      std::complex<double> tmp_2990;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2990 += B0(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPL(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,gI2);
      }
      tmp_2989 += tmp_2990;
      tmp_2988 += (MFe(gI1)) * tmp_2989;
   }
   result += tmp_2988;
   std::complex<double> tmp_2991;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2992;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2992 += B0(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PL(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_2991 += tmp_2992;
   }
   result += tmp_2991;
   std::complex<double> tmp_2993;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2994;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2994 += B0(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPL(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPR(gO1,gI1,
            gI2)*MChi(gI2);
      }
      tmp_2993 += tmp_2994;
   }
   result += tmp_2993;
   std::complex<double> tmp_2995;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2995 += B0(p,MChi(gI1),MSRdp)*Conj(CpbarUCha1barChiSRdpPL(gO2,gI1)
         )*CpbarUCha1barChiSRdpPR(gO1,gI1)*MChi(gI1);
   }
   result += tmp_2995;
   std::complex<double> tmp_2996;
   std::complex<double> tmp_2997;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2998;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2998 += B0(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPL(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2997 += tmp_2998;
   }
   tmp_2996 += tmp_2997;
   result += (3) * tmp_2996;
   std::complex<double> tmp_2999;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3000;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3000 += B0(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPL(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2999 += tmp_3000;
   }
   result += tmp_2999;
   std::complex<double> tmp_3001;
   std::complex<double> tmp_3002;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3002 += B0(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PR(gO2,gI2))*
         CpbarUCha1VPCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_3001 += tmp_3002;
   result += (-4) * tmp_3001;
   std::complex<double> tmp_3003;
   std::complex<double> tmp_3004;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3004 += B0(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PR(gO2,gI2))*
         CpbarUCha1VZCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_3003 += tmp_3004;
   result += (-4) * tmp_3003;
   std::complex<double> tmp_3005;
   std::complex<double> tmp_3006;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3006 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPR(gO2,gI2))
         *CpbarUCha1conjVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3005 += tmp_3006;
   result += (-4) * tmp_3005;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3007;
   std::complex<double> tmp_3008;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3009;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3009 += B1(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPR
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(gO1,gI1,gI2);
      }
      tmp_3008 += tmp_3009;
   }
   tmp_3007 += tmp_3008;
   result += (-0.5) * tmp_3007;
   std::complex<double> tmp_3010;
   std::complex<double> tmp_3011;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3012;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3012 += B1(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPR(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,gI2);
      }
      tmp_3011 += tmp_3012;
   }
   tmp_3010 += tmp_3011;
   result += (-0.5) * tmp_3010;
   std::complex<double> tmp_3013;
   std::complex<double> tmp_3014;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3015;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3015 += B1(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPR(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,gI2);
      }
      tmp_3014 += tmp_3015;
   }
   tmp_3013 += tmp_3014;
   result += (-0.5) * tmp_3013;
   std::complex<double> tmp_3016;
   std::complex<double> tmp_3017;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3018;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3018 += B1(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPR(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,gI2);
      }
      tmp_3017 += tmp_3018;
   }
   tmp_3016 += tmp_3017;
   result += (-1.5) * tmp_3016;
   std::complex<double> tmp_3019;
   std::complex<double> tmp_3020;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3020 += B1(p,MChi(gI1),MSRdp)*Conj(CpbarUCha1barChiSRdpPR(gO2,gI1)
         )*CpbarUCha1barChiSRdpPR(gO1,gI1);
   }
   tmp_3019 += tmp_3020;
   result += (-0.5) * tmp_3019;
   std::complex<double> tmp_3021;
   std::complex<double> tmp_3022;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3023;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3023 += B1(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PR(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PR(gO1,gI1,gI2);
      }
      tmp_3022 += tmp_3023;
   }
   tmp_3021 += tmp_3022;
   result += (-0.5) * tmp_3021;
   std::complex<double> tmp_3024;
   std::complex<double> tmp_3025;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3026;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3026 += B1(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPR(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPR(gO1,gI1,
            gI2);
      }
      tmp_3025 += tmp_3026;
   }
   tmp_3024 += tmp_3025;
   result += (-0.5) * tmp_3024;
   std::complex<double> tmp_3027;
   std::complex<double> tmp_3028;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3029;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3029 += B1(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPR(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPR(gO1,gI1,gI2);
      }
      tmp_3028 += tmp_3029;
   }
   tmp_3027 += tmp_3028;
   result += (-1.5) * tmp_3027;
   std::complex<double> tmp_3030;
   std::complex<double> tmp_3031;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3032;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3032 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPR(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPR(gO1,gI1,gI2);
      }
      tmp_3031 += tmp_3032;
   }
   tmp_3030 += tmp_3031;
   result += (-0.5) * tmp_3030;
   std::complex<double> tmp_3033;
   std::complex<double> tmp_3034;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3034 += B1(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PL(gO2,gI2))*
         CpbarUCha1VPCha1PL(gO1,gI2);
   }
   tmp_3033 += tmp_3034;
   result += (-1) * tmp_3033;
   std::complex<double> tmp_3035;
   std::complex<double> tmp_3036;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3036 += B1(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PL(gO2,gI2))*
         CpbarUCha1VZCha1PL(gO1,gI2);
   }
   tmp_3035 += tmp_3036;
   result += (-1) * tmp_3035;
   std::complex<double> tmp_3037;
   std::complex<double> tmp_3038;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3038 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPL(gO2,gI2))
         *CpbarUCha1conjVWmChiPL(gO1,gI2);
   }
   tmp_3037 += tmp_3038;
   result += (-1) * tmp_3037;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3039;
   std::complex<double> tmp_3040;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3041;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3041 += B1(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPL
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPL(gO1,gI1,gI2);
      }
      tmp_3040 += tmp_3041;
   }
   tmp_3039 += tmp_3040;
   result += (-0.5) * tmp_3039;
   std::complex<double> tmp_3042;
   std::complex<double> tmp_3043;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3044;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3044 += B1(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPL(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPL(gO1,gI1,gI2);
      }
      tmp_3043 += tmp_3044;
   }
   tmp_3042 += tmp_3043;
   result += (-0.5) * tmp_3042;
   std::complex<double> tmp_3045;
   std::complex<double> tmp_3046;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3047;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3047 += B1(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPL(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPL(gO1,gI1,gI2);
      }
      tmp_3046 += tmp_3047;
   }
   tmp_3045 += tmp_3046;
   result += (-0.5) * tmp_3045;
   std::complex<double> tmp_3048;
   std::complex<double> tmp_3049;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3050;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3050 += B1(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPL(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPL(gO1,gI1,gI2);
      }
      tmp_3049 += tmp_3050;
   }
   tmp_3048 += tmp_3049;
   result += (-1.5) * tmp_3048;
   std::complex<double> tmp_3051;
   std::complex<double> tmp_3052;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3052 += B1(p,MChi(gI1),MSRdp)*Conj(CpbarUCha1barChiSRdpPL(gO2,gI1)
         )*CpbarUCha1barChiSRdpPL(gO1,gI1);
   }
   tmp_3051 += tmp_3052;
   result += (-0.5) * tmp_3051;
   std::complex<double> tmp_3053;
   std::complex<double> tmp_3054;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3055;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3055 += B1(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PL(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PL(gO1,gI1,gI2);
      }
      tmp_3054 += tmp_3055;
   }
   tmp_3053 += tmp_3054;
   result += (-0.5) * tmp_3053;
   std::complex<double> tmp_3056;
   std::complex<double> tmp_3057;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3058;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3058 += B1(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPL(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPL(gO1,gI1,
            gI2);
      }
      tmp_3057 += tmp_3058;
   }
   tmp_3056 += tmp_3057;
   result += (-0.5) * tmp_3056;
   std::complex<double> tmp_3059;
   std::complex<double> tmp_3060;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3061;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3061 += B1(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPL(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPL(gO1,gI1,gI2);
      }
      tmp_3060 += tmp_3061;
   }
   tmp_3059 += tmp_3060;
   result += (-1.5) * tmp_3059;
   std::complex<double> tmp_3062;
   std::complex<double> tmp_3063;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3064;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3064 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPL(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPL(gO1,gI1,gI2);
      }
      tmp_3063 += tmp_3064;
   }
   tmp_3062 += tmp_3063;
   result += (-0.5) * tmp_3062;
   std::complex<double> tmp_3065;
   std::complex<double> tmp_3066;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3066 += B1(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PR(gO2,gI2))*
         CpbarUCha1VPCha1PR(gO1,gI2);
   }
   tmp_3065 += tmp_3066;
   result += (-1) * tmp_3065;
   std::complex<double> tmp_3067;
   std::complex<double> tmp_3068;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3068 += B1(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PR(gO2,gI2))*
         CpbarUCha1VZCha1PR(gO1,gI2);
   }
   tmp_3067 += tmp_3068;
   result += (-1) * tmp_3067;
   std::complex<double> tmp_3069;
   std::complex<double> tmp_3070;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3070 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPR(gO2,gI2))
         *CpbarUCha1conjVWmChiPR(gO1,gI2);
   }
   tmp_3069 += tmp_3070;
   result += (-1) * tmp_3069;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3071;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3072;
      std::complex<double> tmp_3073;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3073 += B0(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPL
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPR(gO1,gI1,gI2);
      }
      tmp_3072 += tmp_3073;
      tmp_3071 += (MCha1(gI1)) * tmp_3072;
   }
   result += tmp_3071;
   std::complex<double> tmp_3074;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3075;
      std::complex<double> tmp_3076;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3076 += B0(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPL(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,gI2);
      }
      tmp_3075 += tmp_3076;
      tmp_3074 += (MCha2(gI1)) * tmp_3075;
   }
   result += tmp_3074;
   std::complex<double> tmp_3077;
   std::complex<double> tmp_3078;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3079;
      std::complex<double> tmp_3080;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3080 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPL(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,gI2);
      }
      tmp_3079 += tmp_3080;
      tmp_3078 += (MFu(gI1)) * tmp_3079;
   }
   tmp_3077 += tmp_3078;
   result += (3) * tmp_3077;
   std::complex<double> tmp_3081;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3082;
      std::complex<double> tmp_3083;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3083 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePL(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,gI2);
      }
      tmp_3082 += tmp_3083;
      tmp_3081 += (MFv(gI1)) * tmp_3082;
   }
   result += tmp_3081;
   std::complex<double> tmp_3084;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3085;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3085 += B0(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PL(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3084 += tmp_3085;
   }
   result += tmp_3084;
   std::complex<double> tmp_3086;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3087;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3087 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPL(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3086 += tmp_3087;
   }
   result += tmp_3086;
   std::complex<double> tmp_3088;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3088 += B0(p,MChi(gI1),MSRum)*Conj(CpbarUCha2barChiSRumPL(gO2,gI1)
         )*CpbarUCha2barChiSRumPR(gO1,gI1)*MChi(gI1);
   }
   result += tmp_3088;
   std::complex<double> tmp_3089;
   std::complex<double> tmp_3090;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3091;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3091 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPL(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3090 += tmp_3091;
   }
   tmp_3089 += tmp_3090;
   result += (3) * tmp_3089;
   std::complex<double> tmp_3092;
   std::complex<double> tmp_3093;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3093 += B0(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PR(gO2,gI2))*
         CpbarUCha2VPCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_3092 += tmp_3093;
   result += (-4) * tmp_3092;
   std::complex<double> tmp_3094;
   std::complex<double> tmp_3095;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3095 += B0(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PR(gO2,gI2))*
         CpbarUCha2VZCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_3094 += tmp_3095;
   result += (-4) * tmp_3094;
   std::complex<double> tmp_3096;
   std::complex<double> tmp_3097;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3097 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPR(gO2,gI2))*
         CpbarUCha2VWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3096 += tmp_3097;
   result += (-4) * tmp_3096;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3098;
   std::complex<double> tmp_3099;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3100;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3100 += B1(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPR
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPR(gO1,gI1,gI2);
      }
      tmp_3099 += tmp_3100;
   }
   tmp_3098 += tmp_3099;
   result += (-0.5) * tmp_3098;
   std::complex<double> tmp_3101;
   std::complex<double> tmp_3102;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3103;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3103 += B1(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPR(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,gI2);
      }
      tmp_3102 += tmp_3103;
   }
   tmp_3101 += tmp_3102;
   result += (-0.5) * tmp_3101;
   std::complex<double> tmp_3104;
   std::complex<double> tmp_3105;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3106;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3106 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPR(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,gI2);
      }
      tmp_3105 += tmp_3106;
   }
   tmp_3104 += tmp_3105;
   result += (-1.5) * tmp_3104;
   std::complex<double> tmp_3107;
   std::complex<double> tmp_3108;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3109;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3109 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePR(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,gI2);
      }
      tmp_3108 += tmp_3109;
   }
   tmp_3107 += tmp_3108;
   result += (-0.5) * tmp_3107;
   std::complex<double> tmp_3110;
   std::complex<double> tmp_3111;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3111 += B1(p,MChi(gI1),MSRum)*Conj(CpbarUCha2barChiSRumPR(gO2,gI1)
         )*CpbarUCha2barChiSRumPR(gO1,gI1);
   }
   tmp_3110 += tmp_3111;
   result += (-0.5) * tmp_3110;
   std::complex<double> tmp_3112;
   std::complex<double> tmp_3113;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3114;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3114 += B1(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PR(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PR(gO1,gI1,gI2);
      }
      tmp_3113 += tmp_3114;
   }
   tmp_3112 += tmp_3113;
   result += (-0.5) * tmp_3112;
   std::complex<double> tmp_3115;
   std::complex<double> tmp_3116;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3117;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3117 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPR(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPR(gO1,gI1,gI2);
      }
      tmp_3116 += tmp_3117;
   }
   tmp_3115 += tmp_3116;
   result += (-0.5) * tmp_3115;
   std::complex<double> tmp_3118;
   std::complex<double> tmp_3119;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3120;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3120 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPR(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPR(gO1,gI1,gI2);
      }
      tmp_3119 += tmp_3120;
   }
   tmp_3118 += tmp_3119;
   result += (-1.5) * tmp_3118;
   std::complex<double> tmp_3121;
   std::complex<double> tmp_3122;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3122 += B1(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PL(gO2,gI2))*
         CpbarUCha2VPCha2PL(gO1,gI2);
   }
   tmp_3121 += tmp_3122;
   result += (-1) * tmp_3121;
   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3124 += B1(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PL(gO2,gI2))*
         CpbarUCha2VZCha2PL(gO1,gI2);
   }
   tmp_3123 += tmp_3124;
   result += (-1) * tmp_3123;
   std::complex<double> tmp_3125;
   std::complex<double> tmp_3126;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3126 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPL(gO2,gI2))*
         CpbarUCha2VWmChiPL(gO1,gI2);
   }
   tmp_3125 += tmp_3126;
   result += (-1) * tmp_3125;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3127;
   std::complex<double> tmp_3128;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3129;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3129 += B1(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPL
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPL(gO1,gI1,gI2);
      }
      tmp_3128 += tmp_3129;
   }
   tmp_3127 += tmp_3128;
   result += (-0.5) * tmp_3127;
   std::complex<double> tmp_3130;
   std::complex<double> tmp_3131;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3132;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3132 += B1(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPL(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPL(gO1,gI1,gI2);
      }
      tmp_3131 += tmp_3132;
   }
   tmp_3130 += tmp_3131;
   result += (-0.5) * tmp_3130;
   std::complex<double> tmp_3133;
   std::complex<double> tmp_3134;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3135;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3135 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPL(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPL(gO1,gI1,gI2);
      }
      tmp_3134 += tmp_3135;
   }
   tmp_3133 += tmp_3134;
   result += (-1.5) * tmp_3133;
   std::complex<double> tmp_3136;
   std::complex<double> tmp_3137;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3138;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3138 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePL(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePL(gO1,gI1,gI2);
      }
      tmp_3137 += tmp_3138;
   }
   tmp_3136 += tmp_3137;
   result += (-0.5) * tmp_3136;
   std::complex<double> tmp_3139;
   std::complex<double> tmp_3140;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3140 += B1(p,MChi(gI1),MSRum)*Conj(CpbarUCha2barChiSRumPL(gO2,gI1)
         )*CpbarUCha2barChiSRumPL(gO1,gI1);
   }
   tmp_3139 += tmp_3140;
   result += (-0.5) * tmp_3139;
   std::complex<double> tmp_3141;
   std::complex<double> tmp_3142;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3143;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3143 += B1(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PL(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PL(gO1,gI1,gI2);
      }
      tmp_3142 += tmp_3143;
   }
   tmp_3141 += tmp_3142;
   result += (-0.5) * tmp_3141;
   std::complex<double> tmp_3144;
   std::complex<double> tmp_3145;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3146;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3146 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPL(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPL(gO1,gI1,gI2);
      }
      tmp_3145 += tmp_3146;
   }
   tmp_3144 += tmp_3145;
   result += (-0.5) * tmp_3144;
   std::complex<double> tmp_3147;
   std::complex<double> tmp_3148;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3149;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3149 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPL(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPL(gO1,gI1,gI2);
      }
      tmp_3148 += tmp_3149;
   }
   tmp_3147 += tmp_3148;
   result += (-1.5) * tmp_3147;
   std::complex<double> tmp_3150;
   std::complex<double> tmp_3151;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3151 += B1(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PR(gO2,gI2))*
         CpbarUCha2VPCha2PR(gO1,gI2);
   }
   tmp_3150 += tmp_3151;
   result += (-1) * tmp_3150;
   std::complex<double> tmp_3152;
   std::complex<double> tmp_3153;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3153 += B1(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PR(gO2,gI2))*
         CpbarUCha2VZCha2PR(gO1,gI2);
   }
   tmp_3152 += tmp_3153;
   result += (-1) * tmp_3152;
   std::complex<double> tmp_3154;
   std::complex<double> tmp_3155;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3155 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPR(gO2,gI2))*
         CpbarUCha2VWmChiPR(gO1,gI2);
   }
   tmp_3154 += tmp_3155;
   result += (-1) * tmp_3154;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3156;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3157;
      std::complex<double> tmp_3158;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3158 += B0(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3157 += tmp_3158;
      tmp_3156 += (MCha1(gI1)) * tmp_3157;
   }
   result += tmp_3156;
   std::complex<double> tmp_3159;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3160;
      std::complex<double> tmp_3161;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3161 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3160 += tmp_3161;
      tmp_3159 += (MFe(gI1)) * tmp_3160;
   }
   result += tmp_3159;
   std::complex<double> tmp_3162;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3163;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3163 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3162 += tmp_3163;
   }
   result += tmp_3162;
   std::complex<double> tmp_3164;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3165;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3165 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3164 += tmp_3165;
   }
   result += tmp_3164;
   std::complex<double> tmp_3166;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3167;
      std::complex<double> tmp_3168;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3168 += B0(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePL(gO2
            ,gI1,gI2))*CpbarUFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3167 += tmp_3168;
      tmp_3166 += (MChi(gI1)) * tmp_3167;
   }
   result += tmp_3166;
   std::complex<double> tmp_3169;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3170;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3170 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3169 += tmp_3170;
   }
   result += tmp_3169;
   std::complex<double> tmp_3171;
   std::complex<double> tmp_3172;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3172 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3171 += tmp_3172;
   result += (-4) * tmp_3171;
   std::complex<double> tmp_3173;
   std::complex<double> tmp_3174;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3174 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3173 += tmp_3174;
   result += (-4) * tmp_3173;
   std::complex<double> tmp_3175;
   std::complex<double> tmp_3176;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3176 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3175 += tmp_3176;
   result += (-4) * tmp_3175;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3177;
   std::complex<double> tmp_3178;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3179;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3179 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPR(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3178 += tmp_3179;
   }
   tmp_3177 += tmp_3178;
   result += (-0.5) * tmp_3177;
   std::complex<double> tmp_3180;
   std::complex<double> tmp_3181;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3182;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3182 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3181 += tmp_3182;
   }
   tmp_3180 += tmp_3181;
   result += (-0.5) * tmp_3180;
   std::complex<double> tmp_3183;
   std::complex<double> tmp_3184;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3185;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3185 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_3184 += tmp_3185;
   }
   tmp_3183 += tmp_3184;
   result += (-0.5) * tmp_3183;
   std::complex<double> tmp_3186;
   std::complex<double> tmp_3187;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3188;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3188 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3187 += tmp_3188;
   }
   tmp_3186 += tmp_3187;
   result += (-0.5) * tmp_3186;
   std::complex<double> tmp_3189;
   std::complex<double> tmp_3190;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3191;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3191 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePR(gO2
            ,gI1,gI2))*CpbarUFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3190 += tmp_3191;
   }
   tmp_3189 += tmp_3190;
   result += (-0.5) * tmp_3189;
   std::complex<double> tmp_3192;
   std::complex<double> tmp_3193;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3194;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3194 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3193 += tmp_3194;
   }
   tmp_3192 += tmp_3193;
   result += (-0.5) * tmp_3192;
   std::complex<double> tmp_3195;
   std::complex<double> tmp_3196;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3196 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_3195 += tmp_3196;
   result += (-1) * tmp_3195;
   std::complex<double> tmp_3197;
   std::complex<double> tmp_3198;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3198 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_3197 += tmp_3198;
   result += (-1) * tmp_3197;
   std::complex<double> tmp_3199;
   std::complex<double> tmp_3200;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3200 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_3199 += tmp_3200;
   result += (-1) * tmp_3199;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3201;
   std::complex<double> tmp_3202;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3203;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3203 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPL(gO1,gI1,gI2);
      }
      tmp_3202 += tmp_3203;
   }
   tmp_3201 += tmp_3202;
   result += (-0.5) * tmp_3201;
   std::complex<double> tmp_3204;
   std::complex<double> tmp_3205;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3206;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3206 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3205 += tmp_3206;
   }
   tmp_3204 += tmp_3205;
   result += (-0.5) * tmp_3204;
   std::complex<double> tmp_3207;
   std::complex<double> tmp_3208;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3209;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3209 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
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
         tmp_3212 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3211 += tmp_3212;
   }
   tmp_3210 += tmp_3211;
   result += (-0.5) * tmp_3210;
   std::complex<double> tmp_3213;
   std::complex<double> tmp_3214;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3215;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3215 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePL(gO2
            ,gI1,gI2))*CpbarUFebarChiSePL(gO1,gI1,gI2);
      }
      tmp_3214 += tmp_3215;
   }
   tmp_3213 += tmp_3214;
   result += (-0.5) * tmp_3213;
   std::complex<double> tmp_3216;
   std::complex<double> tmp_3217;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3218;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3218 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3217 += tmp_3218;
   }
   tmp_3216 += tmp_3217;
   result += (-0.5) * tmp_3216;
   std::complex<double> tmp_3219;
   std::complex<double> tmp_3220;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3220 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_3219 += tmp_3220;
   result += (-1) * tmp_3219;
   std::complex<double> tmp_3221;
   std::complex<double> tmp_3222;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3222 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_3221 += tmp_3222;
   result += (-1) * tmp_3221;
   std::complex<double> tmp_3223;
   std::complex<double> tmp_3224;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3224 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_3223 += tmp_3224;
   result += (-1) * tmp_3223;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3225;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3226;
      std::complex<double> tmp_3227;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3227 += B0(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3226 += tmp_3227;
      tmp_3225 += (MCha1(gI1)) * tmp_3226;
   }
   result += tmp_3225;
   std::complex<double> tmp_3228;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3229;
      std::complex<double> tmp_3230;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3230 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3229 += tmp_3230;
      tmp_3228 += (MFd(gI1)) * tmp_3229;
   }
   result += tmp_3228;
   std::complex<double> tmp_3231;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3232;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3232 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3231 += tmp_3232;
   }
   result += tmp_3231;
   std::complex<double> tmp_3233;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3234;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3234 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3233 += tmp_3234;
   }
   result += tmp_3233;
   std::complex<double> tmp_3235;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3236;
      std::complex<double> tmp_3237;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3237 += B0(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPL(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3236 += tmp_3237;
      tmp_3235 += (MChi(gI1)) * tmp_3236;
   }
   result += tmp_3235;
   std::complex<double> tmp_3238;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3239;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3239 += B0(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PL(gO2,
            gI1,gI2))*CpbarUFdSuCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3238 += tmp_3239;
   }
   result += tmp_3238;
   std::complex<double> tmp_3240;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3241;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3241 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3240 += tmp_3241;
   }
   result += tmp_3240;
   std::complex<double> tmp_3242;
   std::complex<double> tmp_3243;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3243 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3242 += tmp_3243;
   result += (-5.333333333333333) * tmp_3242;
   std::complex<double> tmp_3244;
   std::complex<double> tmp_3245;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3245 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3244 += tmp_3245;
   result += (-4) * tmp_3244;
   std::complex<double> tmp_3246;
   std::complex<double> tmp_3247;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3247 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3246 += tmp_3247;
   result += (-4) * tmp_3246;
   std::complex<double> tmp_3248;
   std::complex<double> tmp_3249;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3249 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3248 += tmp_3249;
   result += (-4) * tmp_3248;
   std::complex<double> tmp_3250;
   std::complex<double> tmp_3251;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3251 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,0))*
         CpbarUFdSdGluPR(gO1,gI1,0);
   }
   tmp_3250 += tmp_3251;
   result += (1.3333333333333333*MGlu) * tmp_3250;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3253 += B0(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPL(gO2,0,gI2))*
         CpbarUFdbarGluSdPR(gO1,0,gI2);
   }
   tmp_3252 += tmp_3253;
   result += (1.3333333333333333*MGlu) * tmp_3252;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3254;
   std::complex<double> tmp_3255;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3256;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3256 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPR(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3255 += tmp_3256;
   }
   tmp_3254 += tmp_3255;
   result += (-0.5) * tmp_3254;
   std::complex<double> tmp_3257;
   std::complex<double> tmp_3258;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3259;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3259 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3258 += tmp_3259;
   }
   tmp_3257 += tmp_3258;
   result += (-0.5) * tmp_3257;
   std::complex<double> tmp_3260;
   std::complex<double> tmp_3261;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3262;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3262 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3261 += tmp_3262;
   }
   tmp_3260 += tmp_3261;
   result += (-0.5) * tmp_3260;
   std::complex<double> tmp_3263;
   std::complex<double> tmp_3264;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3265;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3265 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3264 += tmp_3265;
   }
   tmp_3263 += tmp_3264;
   result += (-0.5) * tmp_3263;
   std::complex<double> tmp_3266;
   std::complex<double> tmp_3267;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3268;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3268 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPR(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3267 += tmp_3268;
   }
   tmp_3266 += tmp_3267;
   result += (-0.5) * tmp_3266;
   std::complex<double> tmp_3269;
   std::complex<double> tmp_3270;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3270 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,0))*
         CpbarUFdSdGluPR(gO1,gI1,0);
   }
   tmp_3269 += tmp_3270;
   result += (-0.6666666666666666) * tmp_3269;
   std::complex<double> tmp_3271;
   std::complex<double> tmp_3272;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3273;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3273 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PR(gO2,
            gI1,gI2))*CpbarUFdSuCha2PR(gO1,gI1,gI2);
      }
      tmp_3272 += tmp_3273;
   }
   tmp_3271 += tmp_3272;
   result += (-0.5) * tmp_3271;
   std::complex<double> tmp_3274;
   std::complex<double> tmp_3275;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3276;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3276 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3275 += tmp_3276;
   }
   tmp_3274 += tmp_3275;
   result += (-0.5) * tmp_3274;
   std::complex<double> tmp_3277;
   std::complex<double> tmp_3278;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3278 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_3277 += tmp_3278;
   result += (-1.3333333333333333) * tmp_3277;
   std::complex<double> tmp_3279;
   std::complex<double> tmp_3280;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3280 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_3279 += tmp_3280;
   result += (-1) * tmp_3279;
   std::complex<double> tmp_3281;
   std::complex<double> tmp_3282;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3282 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_3281 += tmp_3282;
   result += (-1) * tmp_3281;
   std::complex<double> tmp_3283;
   std::complex<double> tmp_3284;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3284 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_3283 += tmp_3284;
   result += (-1) * tmp_3283;
   std::complex<double> tmp_3285;
   std::complex<double> tmp_3286;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3286 += B1(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPR(gO2,0,gI2))*
         CpbarUFdbarGluSdPR(gO1,0,gI2);
   }
   tmp_3285 += tmp_3286;
   result += (-0.6666666666666666) * tmp_3285;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3287;
   std::complex<double> tmp_3288;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3289;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3289 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPL(gO1,gI1,gI2);
      }
      tmp_3288 += tmp_3289;
   }
   tmp_3287 += tmp_3288;
   result += (-0.5) * tmp_3287;
   std::complex<double> tmp_3290;
   std::complex<double> tmp_3291;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3292;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3292 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3291 += tmp_3292;
   }
   tmp_3290 += tmp_3291;
   result += (-0.5) * tmp_3290;
   std::complex<double> tmp_3293;
   std::complex<double> tmp_3294;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3295;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3295 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3294 += tmp_3295;
   }
   tmp_3293 += tmp_3294;
   result += (-0.5) * tmp_3293;
   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3298;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3298 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3297 += tmp_3298;
   }
   tmp_3296 += tmp_3297;
   result += (-0.5) * tmp_3296;
   std::complex<double> tmp_3299;
   std::complex<double> tmp_3300;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3301;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3301 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPL(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPL(gO1,gI1,gI2);
      }
      tmp_3300 += tmp_3301;
   }
   tmp_3299 += tmp_3300;
   result += (-0.5) * tmp_3299;
   std::complex<double> tmp_3302;
   std::complex<double> tmp_3303;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3303 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,0))*
         CpbarUFdSdGluPL(gO1,gI1,0);
   }
   tmp_3302 += tmp_3303;
   result += (-0.6666666666666666) * tmp_3302;
   std::complex<double> tmp_3304;
   std::complex<double> tmp_3305;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3306;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3306 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PL(gO2,
            gI1,gI2))*CpbarUFdSuCha2PL(gO1,gI1,gI2);
      }
      tmp_3305 += tmp_3306;
   }
   tmp_3304 += tmp_3305;
   result += (-0.5) * tmp_3304;
   std::complex<double> tmp_3307;
   std::complex<double> tmp_3308;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3309;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3309 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3308 += tmp_3309;
   }
   tmp_3307 += tmp_3308;
   result += (-0.5) * tmp_3307;
   std::complex<double> tmp_3310;
   std::complex<double> tmp_3311;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3311 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_3310 += tmp_3311;
   result += (-1.3333333333333333) * tmp_3310;
   std::complex<double> tmp_3312;
   std::complex<double> tmp_3313;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3313 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_3312 += tmp_3313;
   result += (-1) * tmp_3312;
   std::complex<double> tmp_3314;
   std::complex<double> tmp_3315;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3315 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_3314 += tmp_3315;
   result += (-1) * tmp_3314;
   std::complex<double> tmp_3316;
   std::complex<double> tmp_3317;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3317 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_3316 += tmp_3317;
   result += (-1) * tmp_3316;
   std::complex<double> tmp_3318;
   std::complex<double> tmp_3319;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3319 += B1(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPL(gO2,0,gI2))*
         CpbarUFdbarGluSdPL(gO1,0,gI2);
   }
   tmp_3318 += tmp_3319;
   result += (-0.6666666666666666) * tmp_3318;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3320;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3321;
      std::complex<double> tmp_3322;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3322 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3321 += tmp_3322;
      tmp_3320 += (MCha2(gI1)) * tmp_3321;
   }
   result += tmp_3320;
   std::complex<double> tmp_3323;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3324;
      std::complex<double> tmp_3325;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3325 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3324 += tmp_3325;
      tmp_3323 += (MFu(gI1)) * tmp_3324;
   }
   result += tmp_3323;
   std::complex<double> tmp_3326;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3327;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3327 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3326 += tmp_3327;
   }
   result += tmp_3326;
   std::complex<double> tmp_3328;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3329;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3329 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3328 += tmp_3329;
   }
   result += tmp_3328;
   std::complex<double> tmp_3330;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3331;
      std::complex<double> tmp_3332;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3332 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3331 += tmp_3332;
      tmp_3330 += (MChi(gI1)) * tmp_3331;
   }
   result += tmp_3330;
   std::complex<double> tmp_3333;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3334;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3334 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3333 += tmp_3334;
   }
   result += tmp_3333;
   std::complex<double> tmp_3335;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3336;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3336 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3335 += tmp_3336;
   }
   result += tmp_3335;
   std::complex<double> tmp_3337;
   std::complex<double> tmp_3338;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3338 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3337 += tmp_3338;
   result += (-4) * tmp_3337;
   std::complex<double> tmp_3339;
   std::complex<double> tmp_3340;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3340 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3339 += tmp_3340;
   result += (-5.333333333333333) * tmp_3339;
   std::complex<double> tmp_3341;
   std::complex<double> tmp_3342;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3342 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3341 += tmp_3342;
   result += (-4) * tmp_3341;
   std::complex<double> tmp_3343;
   std::complex<double> tmp_3344;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3344 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3343 += tmp_3344;
   result += (-4) * tmp_3343;
   std::complex<double> tmp_3345;
   std::complex<double> tmp_3346;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3346 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3345 += tmp_3346;
   result += (1.3333333333333333*MGlu) * tmp_3345;
   std::complex<double> tmp_3347;
   std::complex<double> tmp_3348;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3348 += B0(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,0,gI2))*
         CpbarUFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3347 += tmp_3348;
   result += (1.3333333333333333*MGlu) * tmp_3347;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3349;
   std::complex<double> tmp_3350;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3351;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3351 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3350 += tmp_3351;
   }
   tmp_3349 += tmp_3350;
   result += (-0.5) * tmp_3349;
   std::complex<double> tmp_3352;
   std::complex<double> tmp_3353;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3354;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3354 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3353 += tmp_3354;
   }
   tmp_3352 += tmp_3353;
   result += (-0.5) * tmp_3352;
   std::complex<double> tmp_3355;
   std::complex<double> tmp_3356;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3357;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3357 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3356 += tmp_3357;
   }
   tmp_3355 += tmp_3356;
   result += (-0.5) * tmp_3355;
   std::complex<double> tmp_3358;
   std::complex<double> tmp_3359;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3360;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3360 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3359 += tmp_3360;
   }
   tmp_3358 += tmp_3359;
   result += (-0.5) * tmp_3358;
   std::complex<double> tmp_3361;
   std::complex<double> tmp_3362;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3363;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3363 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPR(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3362 += tmp_3363;
   }
   tmp_3361 += tmp_3362;
   result += (-0.5) * tmp_3361;
   std::complex<double> tmp_3364;
   std::complex<double> tmp_3365;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3365 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3364 += tmp_3365;
   result += (-0.6666666666666666) * tmp_3364;
   std::complex<double> tmp_3366;
   std::complex<double> tmp_3367;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3368;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3368 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PR(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2);
      }
      tmp_3367 += tmp_3368;
   }
   tmp_3366 += tmp_3367;
   result += (-0.5) * tmp_3366;
   std::complex<double> tmp_3369;
   std::complex<double> tmp_3370;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3371;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3371 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3370 += tmp_3371;
   }
   tmp_3369 += tmp_3370;
   result += (-0.5) * tmp_3369;
   std::complex<double> tmp_3372;
   std::complex<double> tmp_3373;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3373 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3372 += tmp_3373;
   result += (-1) * tmp_3372;
   std::complex<double> tmp_3374;
   std::complex<double> tmp_3375;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3375 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_3374 += tmp_3375;
   result += (-1.3333333333333333) * tmp_3374;
   std::complex<double> tmp_3376;
   std::complex<double> tmp_3377;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3377 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3376 += tmp_3377;
   result += (-1) * tmp_3376;
   std::complex<double> tmp_3378;
   std::complex<double> tmp_3379;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3379 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3378 += tmp_3379;
   result += (-1) * tmp_3378;
   std::complex<double> tmp_3380;
   std::complex<double> tmp_3381;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3381 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPR(gO2,0,gI2))*
         CpbarUFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3380 += tmp_3381;
   result += (-0.6666666666666666) * tmp_3380;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3382;
   std::complex<double> tmp_3383;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3384;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3384 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPL(gO1,gI1,gI2);
      }
      tmp_3383 += tmp_3384;
   }
   tmp_3382 += tmp_3383;
   result += (-0.5) * tmp_3382;
   std::complex<double> tmp_3385;
   std::complex<double> tmp_3386;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3387;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3387 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3386 += tmp_3387;
   }
   tmp_3385 += tmp_3386;
   result += (-0.5) * tmp_3385;
   std::complex<double> tmp_3388;
   std::complex<double> tmp_3389;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3390;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3390 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3389 += tmp_3390;
   }
   tmp_3388 += tmp_3389;
   result += (-0.5) * tmp_3388;
   std::complex<double> tmp_3391;
   std::complex<double> tmp_3392;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3393;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3393 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3392 += tmp_3393;
   }
   tmp_3391 += tmp_3392;
   result += (-0.5) * tmp_3391;
   std::complex<double> tmp_3394;
   std::complex<double> tmp_3395;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3396;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3396 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3395 += tmp_3396;
   }
   tmp_3394 += tmp_3395;
   result += (-0.5) * tmp_3394;
   std::complex<double> tmp_3397;
   std::complex<double> tmp_3398;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3398 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPL(gO1,gI1,0);
   }
   tmp_3397 += tmp_3398;
   result += (-0.6666666666666666) * tmp_3397;
   std::complex<double> tmp_3399;
   std::complex<double> tmp_3400;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3401;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3401 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PL(gO1,gI1,gI2);
      }
      tmp_3400 += tmp_3401;
   }
   tmp_3399 += tmp_3400;
   result += (-0.5) * tmp_3399;
   std::complex<double> tmp_3402;
   std::complex<double> tmp_3403;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3404;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3404 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3403 += tmp_3404;
   }
   tmp_3402 += tmp_3403;
   result += (-0.5) * tmp_3402;
   std::complex<double> tmp_3405;
   std::complex<double> tmp_3406;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3406 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3405 += tmp_3406;
   result += (-1) * tmp_3405;
   std::complex<double> tmp_3407;
   std::complex<double> tmp_3408;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3408 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_3407 += tmp_3408;
   result += (-1.3333333333333333) * tmp_3407;
   std::complex<double> tmp_3409;
   std::complex<double> tmp_3410;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3410 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3409 += tmp_3410;
   result += (-1) * tmp_3409;
   std::complex<double> tmp_3411;
   std::complex<double> tmp_3412;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3412 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3411 += tmp_3412;
   result += (-1) * tmp_3411;
   std::complex<double> tmp_3413;
   std::complex<double> tmp_3414;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3414 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,0,gI2))*
         CpbarUFubarGluSuPL(gO1,0,gI2);
   }
   tmp_3413 += tmp_3414;
   result += (-0.6666666666666666) * tmp_3413;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3415;
   std::complex<double> tmp_3416;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3417;
      std::complex<double> tmp_3418;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3418 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpbarGlubarFdSdPL(gI1,
            gI2))*CpbarGlubarFdSdPR(gI1,gI2);
      }
      tmp_3417 += tmp_3418;
      tmp_3416 += (MFd(gI1)) * tmp_3417;
   }
   tmp_3415 += tmp_3416;
   result += (0.5) * tmp_3415;
   std::complex<double> tmp_3419;
   std::complex<double> tmp_3420;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3421;
      std::complex<double> tmp_3422;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3422 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpbarGlubarFuSuPL(gI1,
            gI2))*CpbarGlubarFuSuPR(gI1,gI2);
      }
      tmp_3421 += tmp_3422;
      tmp_3420 += (MFu(gI1)) * tmp_3421;
   }
   tmp_3419 += tmp_3420;
   result += (0.5) * tmp_3419;
   std::complex<double> tmp_3423;
   std::complex<double> tmp_3424;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3425;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3425 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpbarGluconjSdFdPL(gI1,
            gI2))*CpbarGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_3424 += tmp_3425;
   }
   tmp_3423 += tmp_3424;
   result += (0.5) * tmp_3423;
   std::complex<double> tmp_3426;
   std::complex<double> tmp_3427;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3428;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3428 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpbarGluconjSuFuPL(gI1,
            gI2))*CpbarGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_3427 += tmp_3428;
   }
   tmp_3426 += tmp_3427;
   result += (0.5) * tmp_3426;
   result += 3*MGlu*B0(p,MGlu,MphiO)*Conj(CpbarGluphiOGluPL())*
      CpbarGluphiOGluPR();
   result += 3*MGlu*B0(p,MGlu,MsigmaO)*Conj(CpbarGlusigmaOGluPL())*
      CpbarGlusigmaOGluPR();
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpbarGluVGGluPR())*CpbarGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -1.5*AbsSqr(CpbarGluphiOGluPR())*B1(p,MGlu,MphiO);
   result += -1.5*AbsSqr(CpbarGlusigmaOGluPR())*B1(p,MGlu,MsigmaO);
   result += -3*AbsSqr(CpbarGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_3429;
   std::complex<double> tmp_3430;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3431;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3431 += AbsSqr(CpbarGlubarFdSdPR(gI1,gI2))*B1(p,MFd(gI1),MSd
            (gI2));
      }
      tmp_3430 += tmp_3431;
   }
   tmp_3429 += tmp_3430;
   result += (-0.25) * tmp_3429;
   std::complex<double> tmp_3432;
   std::complex<double> tmp_3433;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3434;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3434 += AbsSqr(CpbarGlubarFuSuPR(gI1,gI2))*B1(p,MFu(gI1),MSu
            (gI2));
      }
      tmp_3433 += tmp_3434;
   }
   tmp_3432 += tmp_3433;
   result += (-0.25) * tmp_3432;
   std::complex<double> tmp_3435;
   std::complex<double> tmp_3436;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3437;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3437 += AbsSqr(CpbarGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),
            MSd(gI1));
      }
      tmp_3436 += tmp_3437;
   }
   tmp_3435 += tmp_3436;
   result += (-0.25) * tmp_3435;
   std::complex<double> tmp_3438;
   std::complex<double> tmp_3439;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3440;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3440 += AbsSqr(CpbarGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),
            MSu(gI1));
      }
      tmp_3439 += tmp_3440;
   }
   tmp_3438 += tmp_3439;
   result += (-0.25) * tmp_3438;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -1.5*AbsSqr(CpbarGluphiOGluPL())*B1(p,MGlu,MphiO);
   result += -1.5*AbsSqr(CpbarGlusigmaOGluPL())*B1(p,MGlu,MsigmaO);
   result += -3*AbsSqr(CpbarGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_3441;
   std::complex<double> tmp_3442;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3443;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3443 += AbsSqr(CpbarGlubarFdSdPL(gI1,gI2))*B1(p,MFd(gI1),MSd
            (gI2));
      }
      tmp_3442 += tmp_3443;
   }
   tmp_3441 += tmp_3442;
   result += (-0.25) * tmp_3441;
   std::complex<double> tmp_3444;
   std::complex<double> tmp_3445;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3446;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3446 += AbsSqr(CpbarGlubarFuSuPL(gI1,gI2))*B1(p,MFu(gI1),MSu
            (gI2));
      }
      tmp_3445 += tmp_3446;
   }
   tmp_3444 += tmp_3445;
   result += (-0.25) * tmp_3444;
   std::complex<double> tmp_3447;
   std::complex<double> tmp_3448;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3449;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3449 += AbsSqr(CpbarGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),
            MSd(gI1));
      }
      tmp_3448 += tmp_3449;
   }
   tmp_3447 += tmp_3448;
   result += (-0.25) * tmp_3447;
   std::complex<double> tmp_3450;
   std::complex<double> tmp_3451;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3452;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3452 += AbsSqr(CpbarGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),
            MSu(gI1));
      }
      tmp_3451 += tmp_3452;
   }
   tmp_3450 += tmp_3451;
   result += (-0.25) * tmp_3450;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3453;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3454;
      std::complex<double> tmp_3455;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3455 += B0(p,MCha2(gI1),MSe(gI2))*Conj(CpbarFvbarCha2SePL(
            gO2,gI1,gI2))*CpbarFvbarCha2SePR(gO1,gI1,gI2);
      }
      tmp_3454 += tmp_3455;
      tmp_3453 += (MCha2(gI1)) * tmp_3454;
   }
   result += tmp_3453;
   std::complex<double> tmp_3456;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3457;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3457 += B0(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePL(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3456 += tmp_3457;
   }
   result += tmp_3456;
   std::complex<double> tmp_3458;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3459;
      std::complex<double> tmp_3460;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3460 += B0(p,MChi(gI1),MSv(gI2))*Conj(CpbarFvbarChiSvPL(gO2,
            gI1,gI2))*CpbarFvbarChiSvPR(gO1,gI1,gI2);
      }
      tmp_3459 += tmp_3460;
      tmp_3458 += (MChi(gI1)) * tmp_3459;
   }
   result += tmp_3458;
   std::complex<double> tmp_3461;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3462;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3462 += B0(p,MCha1(gI2),MSe(gI1))*Conj(CpbarFvSeCha1PL(gO2,
            gI1,gI2))*CpbarFvSeCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3461 += tmp_3462;
   }
   result += tmp_3461;
   std::complex<double> tmp_3463;
   std::complex<double> tmp_3464;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3464 += B0(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePR(gO2,gI2))*
         CpbarFvconjVWmFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3463 += tmp_3464;
   result += (-4) * tmp_3463;
   std::complex<double> tmp_3465;
   std::complex<double> tmp_3466;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3466 += B0(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3465 += tmp_3466;
   result += (-4) * tmp_3465;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3467;
   std::complex<double> tmp_3468;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3469;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3469 += B1(p,MCha2(gI1),MSe(gI2))*Conj(CpbarFvbarCha2SePR(
            gO2,gI1,gI2))*CpbarFvbarCha2SePR(gO1,gI1,gI2);
      }
      tmp_3468 += tmp_3469;
   }
   tmp_3467 += tmp_3468;
   result += (-0.5) * tmp_3467;
   std::complex<double> tmp_3470;
   std::complex<double> tmp_3471;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3472;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3472 += B1(p,MChi(gI1),MSv(gI2))*Conj(CpbarFvbarChiSvPR(gO2,
            gI1,gI2))*CpbarFvbarChiSvPR(gO1,gI1,gI2);
      }
      tmp_3471 += tmp_3472;
   }
   tmp_3470 += tmp_3471;
   result += (-0.5) * tmp_3470;
   std::complex<double> tmp_3473;
   std::complex<double> tmp_3474;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3475;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3475 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePR(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePR(gO1,gI1,gI2);
      }
      tmp_3474 += tmp_3475;
   }
   tmp_3473 += tmp_3474;
   result += (-0.5) * tmp_3473;
   std::complex<double> tmp_3476;
   std::complex<double> tmp_3477;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3478;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3478 += B1(p,MCha1(gI2),MSe(gI1))*Conj(CpbarFvSeCha1PR(gO2,
            gI1,gI2))*CpbarFvSeCha1PR(gO1,gI1,gI2);
      }
      tmp_3477 += tmp_3478;
   }
   tmp_3476 += tmp_3477;
   result += (-0.5) * tmp_3476;
   std::complex<double> tmp_3479;
   std::complex<double> tmp_3480;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3480 += B1(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePL(gO2,gI2))*
         CpbarFvconjVWmFePL(gO1,gI2);
   }
   tmp_3479 += tmp_3480;
   result += (-1) * tmp_3479;
   std::complex<double> tmp_3481;
   std::complex<double> tmp_3482;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3482 += B1(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPL(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2);
   }
   tmp_3481 += tmp_3482;
   result += (-1) * tmp_3481;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3483;
   std::complex<double> tmp_3484;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3485;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3485 += B1(p,MCha2(gI1),MSe(gI2))*Conj(CpbarFvbarCha2SePL(
            gO2,gI1,gI2))*CpbarFvbarCha2SePL(gO1,gI1,gI2);
      }
      tmp_3484 += tmp_3485;
   }
   tmp_3483 += tmp_3484;
   result += (-0.5) * tmp_3483;
   std::complex<double> tmp_3486;
   std::complex<double> tmp_3487;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3488;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3488 += B1(p,MChi(gI1),MSv(gI2))*Conj(CpbarFvbarChiSvPL(gO2,
            gI1,gI2))*CpbarFvbarChiSvPL(gO1,gI1,gI2);
      }
      tmp_3487 += tmp_3488;
   }
   tmp_3486 += tmp_3487;
   result += (-0.5) * tmp_3486;
   std::complex<double> tmp_3489;
   std::complex<double> tmp_3490;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3491;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3491 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePL(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePL(gO1,gI1,gI2);
      }
      tmp_3490 += tmp_3491;
   }
   tmp_3489 += tmp_3490;
   result += (-0.5) * tmp_3489;
   std::complex<double> tmp_3492;
   std::complex<double> tmp_3493;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3494;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3494 += B1(p,MCha1(gI2),MSe(gI1))*Conj(CpbarFvSeCha1PL(gO2,
            gI1,gI2))*CpbarFvSeCha1PL(gO1,gI1,gI2);
      }
      tmp_3493 += tmp_3494;
   }
   tmp_3492 += tmp_3493;
   result += (-0.5) * tmp_3492;
   std::complex<double> tmp_3495;
   std::complex<double> tmp_3496;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3496 += B1(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePR(gO2,gI2))*
         CpbarFvconjVWmFePR(gO1,gI2);
   }
   tmp_3495 += tmp_3496;
   result += (-1) * tmp_3495;
   std::complex<double> tmp_3497;
   std::complex<double> tmp_3498;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3498 += B1(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPR(gO1,gI2);
   }
   tmp_3497 += tmp_3498;
   result += (-1) * tmp_3497;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   result += -4*AbsSqr(CpVZconjSRdpSRdp())*B00(p,MSRdp,MSRdp);
   result += -4*AbsSqr(CpVZconjSRumSRum())*B00(p,MSRum,MSRum);
   result += A0(MSRdp)*CpVZVZconjSRdpSRdp();
   result += A0(MSRum)*CpVZVZconjSRumSRum();
   std::complex<double> tmp_3499;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3499 += A0(MRh(gI1))*CpVZVZconjRhRh(gI1,gI1);
   }
   result += tmp_3499;
   std::complex<double> tmp_3500;
   std::complex<double> tmp_3501;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3502;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3502 += AbsSqr(CpVZconjRhRh(gI1,gI2))*B00(p,MRh(gI1),MRh(gI2
            ));
      }
      tmp_3501 += tmp_3502;
   }
   tmp_3500 += tmp_3501;
   result += (-4) * tmp_3500;
   std::complex<double> tmp_3503;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3504;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3504 += (AbsSqr(CpVZbarCha1Cha1PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha1Cha1PR(gI1,gI2)))*H0(p,MCha1(gI1),MCha1(gI2));
         tmp_3504 += 4*B0(p,MCha1(gI1),MCha1(gI2))*MCha1(gI1)*MCha1(gI2)*
            Re(Conj(CpVZbarCha1Cha1PL(gI1,gI2))*CpVZbarCha1Cha1PR(gI1,gI2));
      }
      tmp_3503 += tmp_3504;
   }
   result += tmp_3503;
   std::complex<double> tmp_3505;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3506;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3506 += (AbsSqr(CpVZbarCha2Cha2PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha2Cha2PR(gI1,gI2)))*H0(p,MCha2(gI1),MCha2(gI2));
         tmp_3506 += 4*B0(p,MCha2(gI1),MCha2(gI2))*MCha2(gI1)*MCha2(gI2)*
            Re(Conj(CpVZbarCha2Cha2PL(gI1,gI2))*CpVZbarCha2Cha2PR(gI1,gI2));
      }
      tmp_3505 += tmp_3506;
   }
   result += tmp_3505;
   std::complex<double> tmp_3507;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3507 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_3507;
   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3510;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3510 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_3509 += tmp_3510;
   }
   tmp_3508 += tmp_3509;
   result += (-4) * tmp_3508;
   std::complex<double> tmp_3511;
   std::complex<double> tmp_3512;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3512 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_3511 += tmp_3512;
   result += (0.5) * tmp_3511;
   std::complex<double> tmp_3513;
   std::complex<double> tmp_3514;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3515;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3515 += AbsSqr(CpVZhhAh(gI1,1 + gI2))*B00(p,MAh(1 + gI2),Mhh
            (gI1));
      }
      tmp_3514 += tmp_3515;
   }
   tmp_3513 += tmp_3514;
   result += (-4) * tmp_3513;
   std::complex<double> tmp_3516;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3517;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3517 += (AbsSqr(CpVZbarChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZbarChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3517 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZbarChiChiPL(gI1,gI2))*CpVZbarChiChiPR(gI1,gI2));
      }
      tmp_3516 += tmp_3517;
   }
   result += tmp_3516;
   std::complex<double> tmp_3518;
   std::complex<double> tmp_3519;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3519 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3518 += tmp_3519;
   result += (3) * tmp_3518;
   std::complex<double> tmp_3520;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3520 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3520;
   std::complex<double> tmp_3521;
   std::complex<double> tmp_3522;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3522 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3521 += tmp_3522;
   result += (3) * tmp_3521;
   std::complex<double> tmp_3523;
   std::complex<double> tmp_3524;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3525;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3525 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_3524 += tmp_3525;
   }
   tmp_3523 += tmp_3524;
   result += (-12) * tmp_3523;
   std::complex<double> tmp_3526;
   std::complex<double> tmp_3527;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3528;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3528 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_3527 += tmp_3528;
   }
   tmp_3526 += tmp_3527;
   result += (-4) * tmp_3526;
   std::complex<double> tmp_3529;
   std::complex<double> tmp_3530;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3531;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3531 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_3530 += tmp_3531;
   }
   tmp_3529 += tmp_3530;
   result += (-12) * tmp_3529;
   std::complex<double> tmp_3532;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3532 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3532;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   result += A0(MSRdp)*CpVWmconjVWmconjSRdpSRdp();
   result += A0(MSRum)*CpVWmconjVWmconjSRumSRum();
   std::complex<double> tmp_3533;
   std::complex<double> tmp_3534;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3534 += AbsSqr(CpconjVWmconjRhSRum(gI1))*B00(p,MSRum,MRh(gI1));
   }
   tmp_3533 += tmp_3534;
   result += (-4) * tmp_3533;
   std::complex<double> tmp_3535;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3535 += A0(MRh(gI1))*CpVWmconjVWmconjRhRh(gI1,gI1);
   }
   result += tmp_3535;
   std::complex<double> tmp_3536;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3537;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3537 += (AbsSqr(CpconjVWmbarCha1ChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarCha1ChiPR(gI1,gI2)))*H0(p,MCha1(gI1),MChi(gI2));
         tmp_3537 += 4*B0(p,MCha1(gI1),MChi(gI2))*MCha1(gI1)*MChi(gI2)*Re
            (Conj(CpconjVWmbarCha1ChiPL(gI1,gI2))*CpconjVWmbarCha1ChiPR(gI1,gI2));
      }
      tmp_3536 += tmp_3537;
   }
   result += tmp_3536;
   std::complex<double> tmp_3538;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3538 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3538;
   std::complex<double> tmp_3539;
   std::complex<double> tmp_3540;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3541;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3541 += AbsSqr(CpconjVWmHpmhh(1 + gI1,gI2))*B00(p,Mhh(gI2),
            MHpm(1 + gI1));
      }
      tmp_3540 += tmp_3541;
   }
   tmp_3539 += tmp_3540;
   result += (-4) * tmp_3539;
   std::complex<double> tmp_3542;
   std::complex<double> tmp_3543;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3544;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3544 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3543 += tmp_3544;
   }
   tmp_3542 += tmp_3543;
   result += (-4) * tmp_3542;
   std::complex<double> tmp_3545;
   std::complex<double> tmp_3546;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3546 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3545 += tmp_3546;
   result += (0.5) * tmp_3545;
   std::complex<double> tmp_3547;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3548;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3548 += (AbsSqr(CpconjVWmbarChiCha2PL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarChiCha2PR(gI1,gI2)))*H0(p,MChi(gI1),MCha2(gI2));
         tmp_3548 += 4*B0(p,MChi(gI1),MCha2(gI2))*MCha2(gI2)*MChi(gI1)*Re
            (Conj(CpconjVWmbarChiCha2PL(gI1,gI2))*CpconjVWmbarChiCha2PR(gI1,gI2));
      }
      tmp_3547 += tmp_3548;
   }
   result += tmp_3547;
   std::complex<double> tmp_3549;
   std::complex<double> tmp_3550;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3550 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3549 += tmp_3550;
   result += (3) * tmp_3549;
   std::complex<double> tmp_3551;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3551 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3551;
   std::complex<double> tmp_3552;
   std::complex<double> tmp_3553;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3553 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3552 += tmp_3553;
   result += (3) * tmp_3552;
   std::complex<double> tmp_3554;
   std::complex<double> tmp_3555;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3556;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3556 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3555 += tmp_3556;
   }
   tmp_3554 += tmp_3555;
   result += (-12) * tmp_3554;
   std::complex<double> tmp_3557;
   std::complex<double> tmp_3558;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3558 += AbsSqr(CpconjVWmconjSRdpRh(gI2))*B00(p,MRh(gI2),MSRdp);
   }
   tmp_3557 += tmp_3558;
   result += (-4) * tmp_3557;
   std::complex<double> tmp_3559;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3559 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3559;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3560;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3561;
      std::complex<double> tmp_3562;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3562 += B0(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3561 += tmp_3562;
      tmp_3560 += (MCha1(gI1)) * tmp_3561;
   }
   result += tmp_3560;
   std::complex<double> tmp_3563;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3564;
      std::complex<double> tmp_3565;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3565 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3564 += tmp_3565;
      tmp_3563 += (MFe(gI1)) * tmp_3564;
   }
   result += tmp_3563;
   std::complex<double> tmp_3566;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3567;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3567 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3566 += tmp_3567;
   }
   result += tmp_3566;
   std::complex<double> tmp_3568;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3569;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3569 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3568 += tmp_3569;
   }
   result += tmp_3568;
   std::complex<double> tmp_3570;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3571;
      std::complex<double> tmp_3572;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3572 += B0(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePL(gO2,
            gI1,gI2))*CpbarFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3571 += tmp_3572;
      tmp_3570 += (MChi(gI1)) * tmp_3571;
   }
   result += tmp_3570;
   std::complex<double> tmp_3573;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3574;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3574 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3573 += tmp_3574;
   }
   result += tmp_3573;
   std::complex<double> tmp_3575;
   std::complex<double> tmp_3576;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3576 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3575 += tmp_3576;
   result += (-4) * tmp_3575;
   std::complex<double> tmp_3577;
   std::complex<double> tmp_3578;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3578 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3577 += tmp_3578;
   result += (-4) * tmp_3577;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3579;
   std::complex<double> tmp_3580;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3581;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3581 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPR(
            gO2,gI1,gI2))*CpbarFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3580 += tmp_3581;
   }
   tmp_3579 += tmp_3580;
   result += (-0.5) * tmp_3579;
   std::complex<double> tmp_3582;
   std::complex<double> tmp_3583;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3584;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3584 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3583 += tmp_3584;
   }
   tmp_3582 += tmp_3583;
   result += (-0.5) * tmp_3582;
   std::complex<double> tmp_3585;
   std::complex<double> tmp_3586;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3587;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3587 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_3586 += tmp_3587;
   }
   tmp_3585 += tmp_3586;
   result += (-0.5) * tmp_3585;
   std::complex<double> tmp_3588;
   std::complex<double> tmp_3589;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3590;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3590 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3589 += tmp_3590;
   }
   tmp_3588 += tmp_3589;
   result += (-0.5) * tmp_3588;
   std::complex<double> tmp_3591;
   std::complex<double> tmp_3592;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3593;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3593 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePR(gO2,
            gI1,gI2))*CpbarFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3592 += tmp_3593;
   }
   tmp_3591 += tmp_3592;
   result += (-0.5) * tmp_3591;
   std::complex<double> tmp_3594;
   std::complex<double> tmp_3595;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3596;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3596 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3595 += tmp_3596;
   }
   tmp_3594 += tmp_3595;
   result += (-0.5) * tmp_3594;
   std::complex<double> tmp_3597;
   std::complex<double> tmp_3598;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3598 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3597 += tmp_3598;
   result += (-1) * tmp_3597;
   std::complex<double> tmp_3599;
   std::complex<double> tmp_3600;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3600 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3599 += tmp_3600;
   result += (-1) * tmp_3599;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3601;
   std::complex<double> tmp_3602;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3603;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3603 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarFebarCha1SvPL(gO1,gI1,gI2);
      }
      tmp_3602 += tmp_3603;
   }
   tmp_3601 += tmp_3602;
   result += (-0.5) * tmp_3601;
   std::complex<double> tmp_3604;
   std::complex<double> tmp_3605;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3606;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3606 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3605 += tmp_3606;
   }
   tmp_3604 += tmp_3605;
   result += (-0.5) * tmp_3604;
   std::complex<double> tmp_3607;
   std::complex<double> tmp_3608;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3609;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3609 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3608 += tmp_3609;
   }
   tmp_3607 += tmp_3608;
   result += (-0.5) * tmp_3607;
   std::complex<double> tmp_3610;
   std::complex<double> tmp_3611;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3612;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3612 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3611 += tmp_3612;
   }
   tmp_3610 += tmp_3611;
   result += (-0.5) * tmp_3610;
   std::complex<double> tmp_3613;
   std::complex<double> tmp_3614;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3615;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3615 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePL(gO2,
            gI1,gI2))*CpbarFebarChiSePL(gO1,gI1,gI2);
      }
      tmp_3614 += tmp_3615;
   }
   tmp_3613 += tmp_3614;
   result += (-0.5) * tmp_3613;
   std::complex<double> tmp_3616;
   std::complex<double> tmp_3617;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3618;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3618 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3617 += tmp_3618;
   }
   tmp_3616 += tmp_3617;
   result += (-0.5) * tmp_3616;
   std::complex<double> tmp_3619;
   std::complex<double> tmp_3620;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3620 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3619 += tmp_3620;
   result += (-1) * tmp_3619;
   std::complex<double> tmp_3621;
   std::complex<double> tmp_3622;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3622 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3621 += tmp_3622;
   result += (-1) * tmp_3621;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3623;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3624;
      std::complex<double> tmp_3625;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3625 += B0(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3624 += tmp_3625;
      tmp_3623 += (MCha1(gI1)) * tmp_3624;
   }
   result += tmp_3623;
   std::complex<double> tmp_3626;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3627;
      std::complex<double> tmp_3628;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3628 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3627 += tmp_3628;
      tmp_3626 += (MFd(gI1)) * tmp_3627;
   }
   result += tmp_3626;
   std::complex<double> tmp_3629;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3630;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3630 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3629 += tmp_3630;
   }
   result += tmp_3629;
   std::complex<double> tmp_3631;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3632;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3632 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3631 += tmp_3632;
   }
   result += tmp_3631;
   std::complex<double> tmp_3633;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3634;
      std::complex<double> tmp_3635;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3635 += B0(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPL(gO2,
            gI1,gI2))*CpbarFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3634 += tmp_3635;
      tmp_3633 += (MChi(gI1)) * tmp_3634;
   }
   result += tmp_3633;
   std::complex<double> tmp_3636;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3637;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3637 += B0(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PL(gO2,
            gI1,gI2))*CpbarFdSuCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3636 += tmp_3637;
   }
   result += tmp_3636;
   std::complex<double> tmp_3638;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3639;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3639 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3638 += tmp_3639;
   }
   result += tmp_3638;
   std::complex<double> tmp_3640;
   std::complex<double> tmp_3641;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3641 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3640 += tmp_3641;
   result += (-4) * tmp_3640;
   std::complex<double> tmp_3642;
   std::complex<double> tmp_3643;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3643 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3642 += tmp_3643;
   result += (-4) * tmp_3642;
   std::complex<double> tmp_3644;
   std::complex<double> tmp_3645;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3645 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,0))*
         CpbarFdSdGluPR(gO1,gI1,0);
   }
   tmp_3644 += tmp_3645;
   result += (1.3333333333333333*MGlu) * tmp_3644;
   std::complex<double> tmp_3646;
   std::complex<double> tmp_3647;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3647 += B0(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPL(gO2,0,gI2))*
         CpbarFdbarGluSdPR(gO1,0,gI2);
   }
   tmp_3646 += tmp_3647;
   result += (1.3333333333333333*MGlu) * tmp_3646;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3648;
   std::complex<double> tmp_3649;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3650;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3650 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPR(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3649 += tmp_3650;
   }
   tmp_3648 += tmp_3649;
   result += (-0.5) * tmp_3648;
   std::complex<double> tmp_3651;
   std::complex<double> tmp_3652;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3653;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3653 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3652 += tmp_3653;
   }
   tmp_3651 += tmp_3652;
   result += (-0.5) * tmp_3651;
   std::complex<double> tmp_3654;
   std::complex<double> tmp_3655;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3656;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3656 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3655 += tmp_3656;
   }
   tmp_3654 += tmp_3655;
   result += (-0.5) * tmp_3654;
   std::complex<double> tmp_3657;
   std::complex<double> tmp_3658;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3659;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3659 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3658 += tmp_3659;
   }
   tmp_3657 += tmp_3658;
   result += (-0.5) * tmp_3657;
   std::complex<double> tmp_3660;
   std::complex<double> tmp_3661;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3662;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3662 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPR(gO2,
            gI1,gI2))*CpbarFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3661 += tmp_3662;
   }
   tmp_3660 += tmp_3661;
   result += (-0.5) * tmp_3660;
   std::complex<double> tmp_3663;
   std::complex<double> tmp_3664;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3664 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,0))*
         CpbarFdSdGluPR(gO1,gI1,0);
   }
   tmp_3663 += tmp_3664;
   result += (-0.6666666666666666) * tmp_3663;
   std::complex<double> tmp_3665;
   std::complex<double> tmp_3666;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3667;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3667 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PR(gO2,
            gI1,gI2))*CpbarFdSuCha2PR(gO1,gI1,gI2);
      }
      tmp_3666 += tmp_3667;
   }
   tmp_3665 += tmp_3666;
   result += (-0.5) * tmp_3665;
   std::complex<double> tmp_3668;
   std::complex<double> tmp_3669;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3670;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3670 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3669 += tmp_3670;
   }
   tmp_3668 += tmp_3669;
   result += (-0.5) * tmp_3668;
   std::complex<double> tmp_3671;
   std::complex<double> tmp_3672;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3672 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3671 += tmp_3672;
   result += (-1) * tmp_3671;
   std::complex<double> tmp_3673;
   std::complex<double> tmp_3674;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3674 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3673 += tmp_3674;
   result += (-1) * tmp_3673;
   std::complex<double> tmp_3675;
   std::complex<double> tmp_3676;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3676 += B1(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPR(gO2,0,gI2))*
         CpbarFdbarGluSdPR(gO1,0,gI2);
   }
   tmp_3675 += tmp_3676;
   result += (-0.6666666666666666) * tmp_3675;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3677;
   std::complex<double> tmp_3678;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3679;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3679 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPL(gO1,gI1,gI2);
      }
      tmp_3678 += tmp_3679;
   }
   tmp_3677 += tmp_3678;
   result += (-0.5) * tmp_3677;
   std::complex<double> tmp_3680;
   std::complex<double> tmp_3681;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3682;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3682 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3681 += tmp_3682;
   }
   tmp_3680 += tmp_3681;
   result += (-0.5) * tmp_3680;
   std::complex<double> tmp_3683;
   std::complex<double> tmp_3684;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3685;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3685 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3684 += tmp_3685;
   }
   tmp_3683 += tmp_3684;
   result += (-0.5) * tmp_3683;
   std::complex<double> tmp_3686;
   std::complex<double> tmp_3687;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3688;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3688 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3687 += tmp_3688;
   }
   tmp_3686 += tmp_3687;
   result += (-0.5) * tmp_3686;
   std::complex<double> tmp_3689;
   std::complex<double> tmp_3690;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3691;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3691 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPL(gO2,
            gI1,gI2))*CpbarFdbarChiSdPL(gO1,gI1,gI2);
      }
      tmp_3690 += tmp_3691;
   }
   tmp_3689 += tmp_3690;
   result += (-0.5) * tmp_3689;
   std::complex<double> tmp_3692;
   std::complex<double> tmp_3693;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3693 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,0))*
         CpbarFdSdGluPL(gO1,gI1,0);
   }
   tmp_3692 += tmp_3693;
   result += (-0.6666666666666666) * tmp_3692;
   std::complex<double> tmp_3694;
   std::complex<double> tmp_3695;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3696;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3696 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PL(gO2,
            gI1,gI2))*CpbarFdSuCha2PL(gO1,gI1,gI2);
      }
      tmp_3695 += tmp_3696;
   }
   tmp_3694 += tmp_3695;
   result += (-0.5) * tmp_3694;
   std::complex<double> tmp_3697;
   std::complex<double> tmp_3698;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3699;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3699 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3698 += tmp_3699;
   }
   tmp_3697 += tmp_3698;
   result += (-0.5) * tmp_3697;
   std::complex<double> tmp_3700;
   std::complex<double> tmp_3701;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3701 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3700 += tmp_3701;
   result += (-1) * tmp_3700;
   std::complex<double> tmp_3702;
   std::complex<double> tmp_3703;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3703 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3702 += tmp_3703;
   result += (-1) * tmp_3702;
   std::complex<double> tmp_3704;
   std::complex<double> tmp_3705;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3705 += B1(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPL(gO2,0,gI2))*
         CpbarFdbarGluSdPL(gO1,0,gI2);
   }
   tmp_3704 += tmp_3705;
   result += (-0.6666666666666666) * tmp_3704;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3706;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3707;
      std::complex<double> tmp_3708;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3708 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarFubarCha2SdPR(gO1,gI1,gI2);
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
         tmp_3711 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3710 += tmp_3711;
      tmp_3709 += (MFu(gI1)) * tmp_3710;
   }
   result += tmp_3709;
   std::complex<double> tmp_3712;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3713;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3713 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3712 += tmp_3713;
   }
   result += tmp_3712;
   std::complex<double> tmp_3714;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3715;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3715 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3714 += tmp_3715;
   }
   result += tmp_3714;
   std::complex<double> tmp_3716;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3717;
      std::complex<double> tmp_3718;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3718 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPL(gO2,
            gI1,gI2))*CpbarFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3717 += tmp_3718;
      tmp_3716 += (MChi(gI1)) * tmp_3717;
   }
   result += tmp_3716;
   std::complex<double> tmp_3719;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3720;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3720 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PL(gO2,
            gI1,gI2))*CpbarFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3719 += tmp_3720;
   }
   result += tmp_3719;
   std::complex<double> tmp_3721;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3722;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3722 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3721 += tmp_3722;
   }
   result += tmp_3721;
   std::complex<double> tmp_3723;
   std::complex<double> tmp_3724;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3724 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3723 += tmp_3724;
   result += (-4) * tmp_3723;
   std::complex<double> tmp_3725;
   std::complex<double> tmp_3726;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3726 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3725 += tmp_3726;
   result += (-4) * tmp_3725;
   std::complex<double> tmp_3727;
   std::complex<double> tmp_3728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3728 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3727 += tmp_3728;
   result += (-4) * tmp_3727;
   std::complex<double> tmp_3729;
   std::complex<double> tmp_3730;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3730 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,0))*
         CpbarFuSuGluPR(gO1,gI1,0);
   }
   tmp_3729 += tmp_3730;
   result += (1.3333333333333333*MGlu) * tmp_3729;
   std::complex<double> tmp_3731;
   std::complex<double> tmp_3732;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3732 += B0(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPL(gO2,0,gI2))*
         CpbarFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3731 += tmp_3732;
   result += (1.3333333333333333*MGlu) * tmp_3731;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3733;
   std::complex<double> tmp_3734;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3735;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3735 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarFubarCha2SdPR(gO1,gI1,gI2);
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
         tmp_3738 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
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
         tmp_3741 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
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
         tmp_3744 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
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
         tmp_3747 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPR(gO2,
            gI1,gI2))*CpbarFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3746 += tmp_3747;
   }
   tmp_3745 += tmp_3746;
   result += (-0.5) * tmp_3745;
   std::complex<double> tmp_3748;
   std::complex<double> tmp_3749;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3749 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,0))*
         CpbarFuSuGluPR(gO1,gI1,0);
   }
   tmp_3748 += tmp_3749;
   result += (-0.6666666666666666) * tmp_3748;
   std::complex<double> tmp_3750;
   std::complex<double> tmp_3751;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3752;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3752 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PR(gO2,
            gI1,gI2))*CpbarFuSdCha1PR(gO1,gI1,gI2);
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
         tmp_3755 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3754 += tmp_3755;
   }
   tmp_3753 += tmp_3754;
   result += (-0.5) * tmp_3753;
   std::complex<double> tmp_3756;
   std::complex<double> tmp_3757;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3757 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3756 += tmp_3757;
   result += (-1) * tmp_3756;
   std::complex<double> tmp_3758;
   std::complex<double> tmp_3759;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3759 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3758 += tmp_3759;
   result += (-1) * tmp_3758;
   std::complex<double> tmp_3760;
   std::complex<double> tmp_3761;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3761 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3760 += tmp_3761;
   result += (-1) * tmp_3760;
   std::complex<double> tmp_3762;
   std::complex<double> tmp_3763;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3763 += B1(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPR(gO2,0,gI2))*
         CpbarFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3762 += tmp_3763;
   result += (-0.6666666666666666) * tmp_3762;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3764;
   std::complex<double> tmp_3765;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3766;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3766 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarFubarCha2SdPL(gO1,gI1,gI2);
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
         tmp_3769 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
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
         tmp_3772 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
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
         tmp_3775 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
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
         tmp_3778 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPL(gO2,
            gI1,gI2))*CpbarFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3777 += tmp_3778;
   }
   tmp_3776 += tmp_3777;
   result += (-0.5) * tmp_3776;
   std::complex<double> tmp_3779;
   std::complex<double> tmp_3780;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3780 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,0))*
         CpbarFuSuGluPL(gO1,gI1,0);
   }
   tmp_3779 += tmp_3780;
   result += (-0.6666666666666666) * tmp_3779;
   std::complex<double> tmp_3781;
   std::complex<double> tmp_3782;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3783;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3783 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PL(gO2,
            gI1,gI2))*CpbarFuSdCha1PL(gO1,gI1,gI2);
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
         tmp_3786 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3785 += tmp_3786;
   }
   tmp_3784 += tmp_3785;
   result += (-0.5) * tmp_3784;
   std::complex<double> tmp_3787;
   std::complex<double> tmp_3788;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3788 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3787 += tmp_3788;
   result += (-1) * tmp_3787;
   std::complex<double> tmp_3789;
   std::complex<double> tmp_3790;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3790 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3789 += tmp_3790;
   result += (-1) * tmp_3789;
   std::complex<double> tmp_3791;
   std::complex<double> tmp_3792;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3792 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3791 += tmp_3792;
   result += (-1) * tmp_3791;
   std::complex<double> tmp_3793;
   std::complex<double> tmp_3794;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3794 += B1(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPL(gO2,0,gI2))*
         CpbarFubarGluSuPL(gO1,0,gI2);
   }
   tmp_3793 += tmp_3794;
   result += (-0.6666666666666666) * tmp_3793;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3795;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3796;
      std::complex<double> tmp_3797;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3797 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3796 += tmp_3797;
      tmp_3795 += (MCha2(gI1)) * tmp_3796;
   }
   result += tmp_3795;
   std::complex<double> tmp_3798;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3799;
      std::complex<double> tmp_3800;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3800 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3799 += tmp_3800;
      tmp_3798 += (MFu(gI1)) * tmp_3799;
   }
   result += tmp_3798;
   std::complex<double> tmp_3801;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3802;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3802 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3801 += tmp_3802;
   }
   result += tmp_3801;
   std::complex<double> tmp_3803;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3804;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3804 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3803 += tmp_3804;
   }
   result += tmp_3803;
   std::complex<double> tmp_3805;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3806;
      std::complex<double> tmp_3807;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3807 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3806 += tmp_3807;
      tmp_3805 += (MChi(gI1)) * tmp_3806;
   }
   result += tmp_3805;
   std::complex<double> tmp_3808;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3809;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3809 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3808 += tmp_3809;
   }
   result += tmp_3808;
   std::complex<double> tmp_3810;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3811;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3811 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3810 += tmp_3811;
   }
   result += tmp_3810;
   std::complex<double> tmp_3812;
   std::complex<double> tmp_3813;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3813 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3812 += tmp_3813;
   result += (-4) * tmp_3812;
   std::complex<double> tmp_3814;
   std::complex<double> tmp_3815;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3815 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3814 += tmp_3815;
   result += (-4) * tmp_3814;
   std::complex<double> tmp_3816;
   std::complex<double> tmp_3817;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3817 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3816 += tmp_3817;
   result += (-4) * tmp_3816;
   std::complex<double> tmp_3818;
   std::complex<double> tmp_3819;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3819 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3818 += tmp_3819;
   result += (1.3333333333333333*MGlu) * tmp_3818;
   std::complex<double> tmp_3820;
   std::complex<double> tmp_3821;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3821 += B0(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,0,gI2))*
         CpbarUFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3820 += tmp_3821;
   result += (1.3333333333333333*MGlu) * tmp_3820;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3822;
   std::complex<double> tmp_3823;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3824;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3824 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3823 += tmp_3824;
   }
   tmp_3822 += tmp_3823;
   result += (-0.5) * tmp_3822;
   std::complex<double> tmp_3825;
   std::complex<double> tmp_3826;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3827;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3827 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3826 += tmp_3827;
   }
   tmp_3825 += tmp_3826;
   result += (-0.5) * tmp_3825;
   std::complex<double> tmp_3828;
   std::complex<double> tmp_3829;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3830;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3830 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3829 += tmp_3830;
   }
   tmp_3828 += tmp_3829;
   result += (-0.5) * tmp_3828;
   std::complex<double> tmp_3831;
   std::complex<double> tmp_3832;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3833;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3833 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3832 += tmp_3833;
   }
   tmp_3831 += tmp_3832;
   result += (-0.5) * tmp_3831;
   std::complex<double> tmp_3834;
   std::complex<double> tmp_3835;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3836;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3836 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPR(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3835 += tmp_3836;
   }
   tmp_3834 += tmp_3835;
   result += (-0.5) * tmp_3834;
   std::complex<double> tmp_3837;
   std::complex<double> tmp_3838;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3838 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3837 += tmp_3838;
   result += (-0.6666666666666666) * tmp_3837;
   std::complex<double> tmp_3839;
   std::complex<double> tmp_3840;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3841;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3841 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PR(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2);
      }
      tmp_3840 += tmp_3841;
   }
   tmp_3839 += tmp_3840;
   result += (-0.5) * tmp_3839;
   std::complex<double> tmp_3842;
   std::complex<double> tmp_3843;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3844;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3844 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3843 += tmp_3844;
   }
   tmp_3842 += tmp_3843;
   result += (-0.5) * tmp_3842;
   std::complex<double> tmp_3845;
   std::complex<double> tmp_3846;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3846 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3845 += tmp_3846;
   result += (-1) * tmp_3845;
   std::complex<double> tmp_3847;
   std::complex<double> tmp_3848;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3848 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3847 += tmp_3848;
   result += (-1) * tmp_3847;
   std::complex<double> tmp_3849;
   std::complex<double> tmp_3850;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3850 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3849 += tmp_3850;
   result += (-1) * tmp_3849;
   std::complex<double> tmp_3851;
   std::complex<double> tmp_3852;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3852 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPR(gO2,0,gI2))*
         CpbarUFubarGluSuPR(gO1,0,gI2);
   }
   tmp_3851 += tmp_3852;
   result += (-0.6666666666666666) * tmp_3851;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3853;
   std::complex<double> tmp_3854;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3855;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3855 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPL(gO1,gI1,gI2);
      }
      tmp_3854 += tmp_3855;
   }
   tmp_3853 += tmp_3854;
   result += (-0.5) * tmp_3853;
   std::complex<double> tmp_3856;
   std::complex<double> tmp_3857;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3858;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3858 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3857 += tmp_3858;
   }
   tmp_3856 += tmp_3857;
   result += (-0.5) * tmp_3856;
   std::complex<double> tmp_3859;
   std::complex<double> tmp_3860;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3861;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3861 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3860 += tmp_3861;
   }
   tmp_3859 += tmp_3860;
   result += (-0.5) * tmp_3859;
   std::complex<double> tmp_3862;
   std::complex<double> tmp_3863;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3864;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3864 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3863 += tmp_3864;
   }
   tmp_3862 += tmp_3863;
   result += (-0.5) * tmp_3862;
   std::complex<double> tmp_3865;
   std::complex<double> tmp_3866;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3867;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3867 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3866 += tmp_3867;
   }
   tmp_3865 += tmp_3866;
   result += (-0.5) * tmp_3865;
   std::complex<double> tmp_3868;
   std::complex<double> tmp_3869;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3869 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPL(gO1,gI1,0);
   }
   tmp_3868 += tmp_3869;
   result += (-0.6666666666666666) * tmp_3868;
   std::complex<double> tmp_3870;
   std::complex<double> tmp_3871;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3872;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3872 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PL(gO1,gI1,gI2);
      }
      tmp_3871 += tmp_3872;
   }
   tmp_3870 += tmp_3871;
   result += (-0.5) * tmp_3870;
   std::complex<double> tmp_3873;
   std::complex<double> tmp_3874;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3875;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3875 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3874 += tmp_3875;
   }
   tmp_3873 += tmp_3874;
   result += (-0.5) * tmp_3873;
   std::complex<double> tmp_3876;
   std::complex<double> tmp_3877;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3877 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3876 += tmp_3877;
   result += (-1) * tmp_3876;
   std::complex<double> tmp_3878;
   std::complex<double> tmp_3879;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3879 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3878 += tmp_3879;
   result += (-1) * tmp_3878;
   std::complex<double> tmp_3880;
   std::complex<double> tmp_3881;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3881 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3880 += tmp_3881;
   result += (-1) * tmp_3880;
   std::complex<double> tmp_3882;
   std::complex<double> tmp_3883;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3883 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,0,gI2))*
         CpbarUFubarGluSuPL(gO1,0,gI2);
   }
   tmp_3882 += tmp_3883;
   result += (-0.6666666666666666) * tmp_3882;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   result += -(A0(MSRdp)*CpUhhconjSRdpSRdp(gO1));
   result += -(A0(MSRum)*CpUhhconjSRumSRum(gO1));
   result += 4*A0(MVWm)*CpUhhconjVWmVWm(gO1);
   result += 2*A0(MVZ)*CpUhhVZVZ(gO1);
   std::complex<double> tmp_3884;
   std::complex<double> tmp_3885;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3885 += A0(MRh(gI1))*CpUhhconjRhRh(gO1,gI1,gI1);
   }
   tmp_3884 += tmp_3885;
   result += (-1) * tmp_3884;
   std::complex<double> tmp_3886;
   std::complex<double> tmp_3887;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3887 += A0(MCha1(gI1))*(CpUhhbarCha1Cha1PL(gO1,gI1,gI1) +
         CpUhhbarCha1Cha1PR(gO1,gI1,gI1))*MCha1(gI1);
   }
   tmp_3886 += tmp_3887;
   result += (2) * tmp_3886;
   std::complex<double> tmp_3888;
   std::complex<double> tmp_3889;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3889 += A0(MCha2(gI1))*(CpUhhbarCha2Cha2PL(gO1,gI1,gI1) +
         CpUhhbarCha2Cha2PR(gO1,gI1,gI1))*MCha2(gI1);
   }
   tmp_3888 += tmp_3889;
   result += (2) * tmp_3888;
   std::complex<double> tmp_3890;
   std::complex<double> tmp_3891;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3891 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3890 += tmp_3891;
   result += (-1) * tmp_3890;
   std::complex<double> tmp_3892;
   std::complex<double> tmp_3893;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3893 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3892 += tmp_3893;
   result += (6) * tmp_3892;
   std::complex<double> tmp_3894;
   std::complex<double> tmp_3895;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3895 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3894 += tmp_3895;
   result += (2) * tmp_3894;
   std::complex<double> tmp_3896;
   std::complex<double> tmp_3897;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3897 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3896 += tmp_3897;
   result += (6) * tmp_3896;
   std::complex<double> tmp_3898;
   std::complex<double> tmp_3899;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3899 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3898 += tmp_3899;
   result += (-0.5) * tmp_3898;
   std::complex<double> tmp_3900;
   std::complex<double> tmp_3901;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3901 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3900 += tmp_3901;
   result += (-1) * tmp_3900;
   std::complex<double> tmp_3902;
   std::complex<double> tmp_3903;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3903 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3902 += tmp_3903;
   result += (-0.5) * tmp_3902;
   std::complex<double> tmp_3904;
   std::complex<double> tmp_3905;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3905 += A0(MChi(gI1))*(CpUhhbarChiChiPL(gO1,gI1,gI1) +
         CpUhhbarChiChiPR(gO1,gI1,gI1))*MChi(gI1);
   }
   tmp_3904 += tmp_3905;
   result += (2) * tmp_3904;
   std::complex<double> tmp_3906;
   std::complex<double> tmp_3907;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3907 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3906 += tmp_3907;
   result += (-3) * tmp_3906;
   std::complex<double> tmp_3908;
   std::complex<double> tmp_3909;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3909 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3908 += tmp_3909;
   result += (-1) * tmp_3908;
   std::complex<double> tmp_3910;
   std::complex<double> tmp_3911;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3911 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3910 += tmp_3911;
   result += (-3) * tmp_3910;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_phiO() const
{
   std::complex<double> result;

   std::complex<double> tmp_3912;
   std::complex<double> tmp_3913;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3913 += A0(MSd(gI1))*CpphiOconjSdSd(gI1,gI1);
   }
   tmp_3912 += tmp_3913;
   result += (-0.5) * tmp_3912;
   std::complex<double> tmp_3914;
   std::complex<double> tmp_3915;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3915 += A0(MSu(gI1))*CpphiOconjSuSu(gI1,gI1);
   }
   tmp_3914 += tmp_3915;
   result += (-0.5) * tmp_3914;
   result += 6*MGlu*A0(MGlu)*(CpphiObarGluGluPL(0,0) + CpphiObarGluGluPR(0,0));

   return result * oneOver16PiSqr;

}










void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVG) old_MVG(MVG), new_MVG(MVG);

   do {
      PHYSICAL(MVG) = 0.;

      new_MVG = PHYSICAL(MVG);
      diff = MaxRelDiff(new_MVG, old_MVG);
      old_MVG = new_MVG;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::VG);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::VG);
}

void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MGlu) old_MGlu(MGlu), new_MGlu(MGlu);

   do {
      const double M_tree(MGlu);
      const double p = old_MGlu;
      const double self_energy_1  = Re(self_energy_Glu_1(p));
      const double self_energy_PL = Re(self_energy_Glu_PL(p));
      const double self_energy_PR = Re(self_energy_Glu_PR(p));
      const auto M_loop = M_tree - self_energy_1 - M_tree * (
         self_energy_PL + self_energy_PR);
      PHYSICAL(MGlu) = calculate_singlet_mass(M_loop);

      new_MGlu = PHYSICAL(MGlu);
      diff = MaxRelDiff(new_MGlu, old_MGlu);
      old_MGlu = new_MGlu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Glu);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Glu);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFv) old_MFv(MFv), new_MFv(MFv);

   do {
      PHYSICAL(MFv).setConstant(0.);

      new_MFv = PHYSICAL(MFv);
      diff = MaxRelDiff(new_MFv, old_MFv);
      old_MFv = new_MFv;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Fv);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Fv);
}

void CLASSNAME::calculate_MSRdp_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::SRdp))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSRdp) old_MSRdp(MSRdp), new_MSRdp(MSRdp);

   do {
      const double M_tree(get_mass_matrix_SRdp());
      const double p = old_MSRdp;
      double self_energy = Re(self_energy_SRdp(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(MSRdp) = SignedAbsSqrt(mass_sqr);

      new_MSRdp = PHYSICAL(MSRdp);
      diff = MaxRelDiff(new_MSRdp, old_MSRdp);
      old_MSRdp = new_MSRdp;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::SRdp);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::SRdp);
}

void CLASSNAME::calculate_MSRum_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::SRum))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSRum) old_MSRum(MSRum), new_MSRum(MSRum);

   do {
      const double M_tree(get_mass_matrix_SRum());
      const double p = old_MSRum;
      double self_energy = Re(self_energy_SRum(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(MSRum) = SignedAbsSqrt(mass_sqr);

      new_MSRum = PHYSICAL(MSRum);
      diff = MaxRelDiff(new_MSRum, old_MSRum);
      old_MSRum = new_MSRum;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::SRum);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::SRum);
}

void CLASSNAME::calculate_MsigmaO_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::sigmaO))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MsigmaO) old_MsigmaO(MsigmaO), new_MsigmaO(MsigmaO);

   do {
      const double M_tree(get_mass_matrix_sigmaO());
      const double p = old_MsigmaO;
      double self_energy = Re(self_energy_sigmaO(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(MsigmaO) = SignedAbsSqrt(mass_sqr);

      new_MsigmaO = PHYSICAL(MsigmaO);
      diff = MaxRelDiff(new_MsigmaO, old_MsigmaO);
      old_MsigmaO = new_MsigmaO;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::sigmaO);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::sigmaO
         );
}

void CLASSNAME::calculate_MphiO_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::phiO))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MphiO) old_MphiO(MphiO), new_MphiO(MphiO);

   do {
      const double M_tree(get_mass_matrix_phiO());
      const double p = old_MphiO;
      double self_energy = Re(self_energy_phiO(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(MphiO) = SignedAbsSqrt(mass_sqr);

      new_MphiO = PHYSICAL(MphiO);
      diff = MaxRelDiff(new_MphiO, old_MphiO);
      old_MphiO = new_MphiO;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::phiO);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::phiO);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVP) old_MVP(MVP), new_MVP(MVP);

   do {
      PHYSICAL(MVP) = 0.;

      new_MVP = PHYSICAL(MVP);
      diff = MaxRelDiff(new_MVP, old_MVP);
      old_MVP = new_MVP;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::VP);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::VP);
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::VZ))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVZ) old_MVZ(MVZ), new_MVZ(MVZ);

   do {
      const double M_tree(Sqr(MVZ));
      const double p = old_MVZ;
      const double self_energy = Re(self_energy_VZ(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(MRSSMtower_info::VZ);

      PHYSICAL(MVZ) = AbsSqrt(mass_sqr);

      new_MVZ = PHYSICAL(MVZ);
      diff = MaxRelDiff(new_MVZ, old_MVZ);
      old_MVZ = new_MVZ;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::VZ);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::VZ);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Sd))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSd) old_MSd(MSd), new_MSd(MSd);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sd());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSd(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Sd(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZD;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZD, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Sd,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZD);
         #endif

         PHYSICAL(MSd(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZD) = mix_ZD;
      }

      new_MSd = PHYSICAL(MSd);
      diff = MaxRelDiff(new_MSd, old_MSd);
      old_MSd = new_MSd;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Sd);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Sd);
}

void CLASSNAME::calculate_MSv_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Sv))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSv) old_MSv(MSv), new_MSv(MSv);

   do {
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Sv());

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MSv(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Sv(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZV;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZV, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Sv,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZV);
         #endif

         PHYSICAL(MSv(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZV) = mix_ZV;
      }

      new_MSv = PHYSICAL(MSv);
      diff = MaxRelDiff(new_MSv, old_MSv);
      old_MSv = new_MSv;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Sv);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Sv);
}

void CLASSNAME::calculate_MSu_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Su))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSu) old_MSu(MSu), new_MSu(MSu);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Su());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSu(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Su(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZU;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZU, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Su,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZU);
         #endif

         PHYSICAL(MSu(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZU) = mix_ZU;
      }

      new_MSu = PHYSICAL(MSu);
      diff = MaxRelDiff(new_MSu, old_MSu);
      old_MSu = new_MSu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Su);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Su);
}

void CLASSNAME::calculate_MSe_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Se))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSe) old_MSe(MSe), new_MSe(MSe);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Se());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSe(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Se(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZE;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZE, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Se,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZE);
         #endif

         PHYSICAL(MSe(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZE) = mix_ZE;
      }

      new_MSe = PHYSICAL(MSe);
      diff = MaxRelDiff(new_MSe, old_MSe);
      old_MSe = new_MSe;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Se);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Se);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::hh))
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
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::hh,
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
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Ah))
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
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA);
         #endif

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
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Ah);
}

void CLASSNAME::calculate_MRh_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Rh))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MRh) old_MRh(MRh), new_MRh(MRh);

   do {
      Eigen::Matrix<double,2,2> self_energy;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Rh());

      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MRh(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = i1; i2 < 2; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Rh(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,2,2> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZHR;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZHR, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Rh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZHR);
         #endif

         PHYSICAL(MRh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZHR) = mix_ZHR;
      }

      new_MRh = PHYSICAL(MRh);
      diff = MaxRelDiff(new_MRh, old_MRh);
      old_MRh = new_MRh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Rh);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Rh);
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::Hpm))
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
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(MRSSMtower_info::Hpm,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
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
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Hpm);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MChi) old_MChi(MChi), new_MChi(MChi);

   do {
      Eigen::Matrix<double,4,4> self_energy_1;
      Eigen::Matrix<double,4,4> self_energy_PL;
      Eigen::Matrix<double,4,4> self_energy_PR;
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Chi());
      for (unsigned es = 0; es < 4; ++es) {
         const double p = Abs(old_MChi(es));
         for (unsigned i1 = 0; i1 < 4; ++i1) {
            for (unsigned i2 = 0; i2 < 4; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Chi_1(p
                  ,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Chi_PL(
                  p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Chi_PR(
                  p,i1,i2));
            }
         }
         const Eigen::Matrix<double,4,4> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,4,4> M_loop(M_tree + delta_M);
         Eigen::Array<double,4,1> eigen_values;
         decltype(ZN1) mix_ZN1;
         decltype(ZN2) mix_ZN2;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZN1, mix_ZN2,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Chi,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZN1, mix_ZN2);
      #endif
         if (es == 0) {
            PHYSICAL(ZN1) = mix_ZN1;
            PHYSICAL(ZN2) = mix_ZN2;
         }
         PHYSICAL(MChi(es)) = Abs(eigen_values(es));
      }

      new_MChi = PHYSICAL(MChi);
      diff = MaxRelDiff(new_MChi, old_MChi);
      old_MChi = new_MChi;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Chi);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Chi);
}

void CLASSNAME::calculate_MCha1_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MCha1) old_MCha1(MCha1), new_MCha1(MCha1);

   do {
      Eigen::Matrix<double,2,2> self_energy_1;
      Eigen::Matrix<double,2,2> self_energy_PL;
      Eigen::Matrix<double,2,2> self_energy_PR;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha1());
      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MCha1(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = 0; i2 < 2; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Cha1_1(
                  p,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Cha1_PL
                  (p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Cha1_PR
                  (p,i1,i2));
            }
         }
         const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
         Eigen::Array<double,2,1> eigen_values;
         decltype(UM1) mix_UM1;
         decltype(UP1) mix_UP1;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_UM1, mix_UP1,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Cha1,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_UM1, mix_UP1);
      #endif
         if (es == 0) {
            PHYSICAL(UM1) = mix_UM1;
            PHYSICAL(UP1) = mix_UP1;
         }
         PHYSICAL(MCha1(es)) = Abs(eigen_values(es));
      }

      new_MCha1 = PHYSICAL(MCha1);
      diff = MaxRelDiff(new_MCha1, old_MCha1);
      old_MCha1 = new_MCha1;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Cha1);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Cha1);
}

void CLASSNAME::calculate_MCha2_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MCha2) old_MCha2(MCha2), new_MCha2(MCha2);

   do {
      Eigen::Matrix<double,2,2> self_energy_1;
      Eigen::Matrix<double,2,2> self_energy_PL;
      Eigen::Matrix<double,2,2> self_energy_PR;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha2());
      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MCha2(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = 0; i2 < 2; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Cha2_1(
                  p,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Cha2_PL
                  (p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Cha2_PR
                  (p,i1,i2));
            }
         }
         const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
         Eigen::Array<double,2,1> eigen_values;
         decltype(UM2) mix_UM2;
         decltype(UP2) mix_UP2;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_UM2, mix_UP2,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Cha2,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_UM2, mix_UP2);
      #endif
         if (es == 0) {
            PHYSICAL(UM2) = mix_UM2;
            PHYSICAL(UP2) = mix_UP2;
         }
         PHYSICAL(MCha2(es)) = Abs(eigen_values(es));
      }

      new_MCha2 = PHYSICAL(MCha2);
      diff = MaxRelDiff(new_MCha2, old_MCha2);
      old_MCha2 = new_MCha2;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Cha2);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Cha2);
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFe) old_MFe(MFe), new_MFe(MFe);

   do {
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFe(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = 0; i2 < 3; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fe_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fe_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fe_PR(p
                  ,i1,i2));
            }
         }
         const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZEL) mix_ZEL;
         decltype(ZER) mix_ZER;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Fe,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER);
      #endif
         if (es == 0) {
            PHYSICAL(ZEL) = mix_ZEL;
            PHYSICAL(ZER) = mix_ZER;
         }
         PHYSICAL(MFe(es)) = Abs(eigen_values(es));
      }

      new_MFe = PHYSICAL(MFe);
      diff = MaxRelDiff(new_MFe, old_MFe);
      old_MFe = new_MFe;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Fe);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Fe);
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFd) old_MFd(MFd), new_MFd(MFd);

   do {
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFd(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = 0; i2 < 3; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fd_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fd_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fd_PR(p
                  ,i1,i2));
            }
         }
         const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZDL) mix_ZDL;
         decltype(ZDR) mix_ZDR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Fd,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR);
      #endif
         if (es == 0) {
            PHYSICAL(ZDL) = mix_ZDL;
            PHYSICAL(ZDR) = mix_ZDR;
         }
         PHYSICAL(MFd(es)) = Abs(eigen_values(es));
      }

      new_MFd = PHYSICAL(MFd);
      diff = MaxRelDiff(new_MFd, old_MFd);
      old_MFd = new_MFd;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Fd);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Fd);
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFu) old_MFu(MFu), new_MFu(MFu);

   do {
      double qcd_1l = 0.;

      {
         const double currentScale = get_scale();
         qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))
            /Sqr(currentScale)))*Sqr(g3);
      }

      double qcd_2l = 0.;

      if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
         const double currentScale = get_scale();
         qcd_2l = -0.005191204615668296*Power(g3,4) -
            0.0032883224409535764*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(2))
            ) - 0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(currentScale,2)
            /Sqr(MFu(2))));
      }

      double qcd_3l = 0.;

      if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
         const double currentScale = get_scale();
         qcd_3l = 0;
      }

      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFu(es));
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
                  self_energy_1(i1,i2)  = Re(
                     self_energy_Fu_1(p,i1,i2));
                  self_energy_PL(i1,i2) = Re(
                     self_energy_Fu_PL(p,i1,i2));
                  self_energy_PR(i1,i2) = Re(
                     self_energy_Fu_PR(p,i1,i2));
               }
            }
         }
         Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZUL) mix_ZUL;
         decltype(ZUR) mix_ZUR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSMtower_info::Fu,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR);
      #endif
         if (es == 0) {
            PHYSICAL(ZUL) = mix_ZUL;
            PHYSICAL(ZUR) = mix_ZUR;
         }
         PHYSICAL(MFu(es)) = Abs(eigen_values(es));
      }

      new_MFu = PHYSICAL(MFu);
      diff = MaxRelDiff(new_MFu, old_MFu);
      old_MFu = new_MFu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::Fu);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::Fu);
}

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::VWm))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVWm) old_MVWm(MVWm), new_MVWm(MVWm);

   do {
      const double M_tree(Sqr(MVWm));
      const double p = old_MVWm;
      const double self_energy = Re(self_energy_VWm(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(MRSSMtower_info::VWm);

      PHYSICAL(MVWm) = AbsSqrt(mass_sqr);

      new_MVWm = PHYSICAL(MVWm);
      diff = MaxRelDiff(new_MVWm, old_MVWm);
      old_MVWm = new_MVWm;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(MRSSMtower_info::VWm);
   else
      problems.unflag_no_pole_mass_convergence(MRSSMtower_info::VWm);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MRSSMtower_info::VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(MRSSMtower_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(MRSSMtower_info::VZ);

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
   const double qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(idx))
      /Sqr(currentScale)))*Sqr(g3);
   double qcd_2l = 0., qcd_3l = 0.;

   if (get_thresholds() > 1) {
      qcd_2l = -0.003408916029785599*Power(g3,4) -
         0.0011495761378943394*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(currentScale,2)
         /Sqr(MFu(idx))));
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
      problems.flag_tachyon(MRSSMtower_info::VZ);
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
      problems.flag_tachyon(MRSSMtower_info::VWm);
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


std::ostream& operator<<(std::ostream& ostr, const MRSSMtower_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
