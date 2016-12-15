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

// File generated at Thu 15 Dec 2016 13:01:15

/**
 * @file MSSMRHN_mass_eigenstates.cpp
 * @brief implementation of the MSSMRHN model class
 *
 * Contains the definition of the MSSMRHN model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Thu 15 Dec 2016 13:01:15 with FlexibleSUSY
 * 1.7.2 (git commit: 0d19299fef514160cb7541a03abb9b2c3365f927) and SARAH 4.9.1 .
 */

#include "MSSMRHN_mass_eigenstates.hpp"
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

using namespace MSSMRHN_info;

#define CLASSNAME MSSMRHN_mass_eigenstates

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

CLASSNAME::MSSMRHN_mass_eigenstates(const MSSMRHN_input_parameters& input_)
   : MSSMRHN_soft_parameters(input_)
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
   , problems(MSSMRHN_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MGlu(0), MSd(Eigen::Array<double,6,1>::Zero()), MSu(Eigen::Array<
      double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()), MSv(
      Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,2,1>::Zero()),
      MAh(Eigen::Array<double,2,1>::Zero()), MHpm(Eigen::Array<double,2,1>::Zero()
      ), MChi(Eigen::Array<double,4,1>::Zero()), MFv(Eigen::Array<double,6,1>
      ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,3
      ,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double
      ,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZU(Eigen::Matrix<double,6,6>::Zero(
      )), ZE(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,2,2>::Zero()), ZA(Eigen::Matrix<double,2,
      2>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,4,4>::Zero()), UV(Eigen::Matrix<std::complex<double>,6,
      6>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~MSSMRHN_mass_eigenstates()
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

const MSSMRHN_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MSSMRHN_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MSSMRHN_physical& physical_)
{
   physical = physical_;
}

const Problems<MSSMRHN_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<MSSMRHN_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
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
   MSSMRHN_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_BMu(gsl_vector_get(x, 0));
   model->set_Mu(INPUT(SignMu) * Abs(gsl_vector_get(x, 1)));


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

   BMu = solver->get_solution(0);
   Mu = LOCALINPUT(SignMu)*Abs(solver->get_solution(1));


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

   BMu = Re((0.05*(-20*mHd2*vd*vu + 20*mHu2*vd*vu - 3*Power(vd,3)*vu*Sqr(g1) +
      3*vd*Power(vu,3)*Sqr(g1) - 5*Power(vd,3)*vu*Sqr(g2) + 5*vd*Power(vu,3)*Sqr(
      g2)))/(Sqr(vd) - Sqr(vu)));
   Mu = Re(0.15811388300841897*LOCALINPUT(SignMu)*Sqrt((-40*mHd2*vd + 40*vu*BMu
      - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*
      vd*Sqr(g2)*Sqr(vu))/vd));

   const bool is_finite = IsFinite(BMu) && IsFinite(Mu);

   if (!is_finite) {
      BMu = old_BMu;
      Mu = old_Mu;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = 0;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) - 3*Power(
      vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*
      Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) - 3*Power(
      vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*
      Sqr(vd)))/vu);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);

   if (!is_finite) {
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

   x_init[0] = BMu;
   x_init[1] = Mu;


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

   double BMu;
   double Mu;

   BMu = Re((0.05*(-20*mHd2*vd*vu + 20*mHu2*vd*vu + 20*vu*tadpole[0] - 20*vd*
      tadpole[1] - 3*Power(vd,3)*vu*Sqr(g1) + 3*vd*Power(vu,3)*Sqr(g1) - 5*Power(
      vd,3)*vu*Sqr(g2) + 5*vd*Power(vu,3)*Sqr(g2)))/(Sqr(vd) - Sqr(vu)));
   Mu = Re(0.15811388300841897*LOCALINPUT(SignMu)*Sqrt((-40*mHd2*vd + 40*vu*BMu
      + 40*tadpole[0] - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) + 3*vd*Sqr(
      g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu))/vd));

   const bool is_finite = IsFinite(BMu) && IsFinite(Mu);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = BMu;
   ewsb_parameters[1] = Mu;


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
   MSSMRHN_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double BMu = gsl_vector_get(x, 0);
   const double Mu = INPUT(SignMu) * Abs(gsl_vector_get(x, 1));

   model->set_BMu(BMu);
   model->set_Mu(Mu);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = BMu;
   ewsb_parameters[1] = Mu;


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
           "MSSMRHN\n"
           "========================================\n";
   MSSMRHN_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
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
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UV = " << UV << '\n';
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

   solve_ewsb_tree_level_custom();

   calculate_MVPVZ();
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
   calculate_MSv();
   calculate_MSe();
   calculate_MSu();
   calculate_MSd();
   calculate_MGlu();
   calculate_MVG();

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;

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
   std::future<void> fut_MVP;
   std::future<void> fut_MVZ;
   std::future<void> fut_MSd;
   std::future<void> fut_MSu;
   std::future<void> fut_MSe;
   std::future<void> fut_MSv;
   std::future<void> fut_Mhh;
   std::future<void> fut_MAh;
   std::future<void> fut_MHpm;
   std::future<void> fut_MChi;
   std::future<void> fut_MFv;
   std::future<void> fut_MCha;
   std::future<void> fut_MFe;
   std::future<void> fut_MFd;
   std::future<void> fut_MFu;
   std::future<void> fut_MVWm;

   if (calculate_bsm_pole_masses) {
      fut_MAh = run_async([obj_ptr] () { obj_ptr->calculate_MAh_pole(); });
      fut_MCha = run_async([obj_ptr] () { obj_ptr->calculate_MCha_pole(); });
      fut_MChi = run_async([obj_ptr] () { obj_ptr->calculate_MChi_pole(); });
      fut_MGlu = run_async([obj_ptr] () { obj_ptr->calculate_MGlu_pole(); });
      fut_Mhh = run_async([obj_ptr] () { obj_ptr->calculate_Mhh_pole(); });
      fut_MHpm = run_async([obj_ptr] () { obj_ptr->calculate_MHpm_pole(); });
      fut_MSd = run_async([obj_ptr] () { obj_ptr->calculate_MSd_pole(); });
      fut_MSe = run_async([obj_ptr] () { obj_ptr->calculate_MSe_pole(); });
      fut_MSu = run_async([obj_ptr] () { obj_ptr->calculate_MSu_pole(); });
      fut_MSv = run_async([obj_ptr] () { obj_ptr->calculate_MSv_pole(); });
   }

   if (calculate_sm_pole_masses) {
      fut_MVG = run_async([obj_ptr] () { obj_ptr->calculate_MVG_pole(); });
      fut_MVP = run_async([obj_ptr] () { obj_ptr->calculate_MVP_pole(); });
      fut_MVZ = run_async([obj_ptr] () { obj_ptr->calculate_MVZ_pole(); });
      fut_MFv = run_async([obj_ptr] () { obj_ptr->calculate_MFv_pole(); });
      fut_MFe = run_async([obj_ptr] () { obj_ptr->calculate_MFe_pole(); });
      fut_MFd = run_async([obj_ptr] () { obj_ptr->calculate_MFd_pole(); });
      fut_MFu = run_async([obj_ptr] () { obj_ptr->calculate_MFu_pole(); });
      fut_MVWm = run_async([obj_ptr] () { obj_ptr->calculate_MVWm_pole(); });
   }

   if (fut_MAh.valid()) fut_MAh.get();
   if (fut_MCha.valid()) fut_MCha.get();
   if (fut_MChi.valid()) fut_MChi.get();
   if (fut_MGlu.valid()) fut_MGlu.get();
   if (fut_Mhh.valid()) fut_Mhh.get();
   if (fut_MHpm.valid()) fut_MHpm.get();
   if (fut_MSd.valid()) fut_MSd.get();
   if (fut_MSe.valid()) fut_MSe.get();
   if (fut_MSu.valid()) fut_MSu.get();
   if (fut_MSv.valid()) fut_MSv.get();
   if (fut_MVG.valid()) fut_MVG.get();
   if (fut_MVP.valid()) fut_MVP.get();
   if (fut_MVZ.valid()) fut_MVZ.get();
   if (fut_MFv.valid()) fut_MFv.get();
   if (fut_MFe.valid()) fut_MFe.get();
   if (fut_MFd.valid()) fut_MFd.get();
   if (fut_MFu.valid()) fut_MFu.get();
   if (fut_MVWm.valid()) fut_MVWm.get();

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
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSv) = MSv;
   PHYSICAL(ZV) = ZV;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(UV) = UV;
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Sd);
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Su);
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Se);
   if (PHYSICAL(MSv).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Sv);
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) problems.flag_tachyon(hh);
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
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,6,1>::Zero();
   ZV = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MFv = Eigen::Matrix<double,6,1>::Zero();
   UV = Eigen::Matrix<std::complex<double>,6,6>::Zero();
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

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   MSSMRHN_soft_parameters::clear();
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
   MSu(0) = pars(8);
   MSu(1) = pars(9);
   MSu(2) = pars(10);
   MSu(3) = pars(11);
   MSu(4) = pars(12);
   MSu(5) = pars(13);
   MSe(0) = pars(14);
   MSe(1) = pars(15);
   MSe(2) = pars(16);
   MSe(3) = pars(17);
   MSe(4) = pars(18);
   MSe(5) = pars(19);
   MSv(0) = pars(20);
   MSv(1) = pars(21);
   MSv(2) = pars(22);
   MSv(3) = pars(23);
   MSv(4) = pars(24);
   MSv(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   MAh(0) = pars(28);
   MAh(1) = pars(29);
   MHpm(0) = pars(30);
   MHpm(1) = pars(31);
   MChi(0) = pars(32);
   MChi(1) = pars(33);
   MChi(2) = pars(34);
   MChi(3) = pars(35);
   MFv(0) = pars(36);
   MFv(1) = pars(37);
   MFv(2) = pars(38);
   MFv(3) = pars(39);
   MFv(4) = pars(40);
   MFv(5) = pars(41);
   MCha(0) = pars(42);
   MCha(1) = pars(43);
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
   pars(2) = MSd(0);
   pars(3) = MSd(1);
   pars(4) = MSd(2);
   pars(5) = MSd(3);
   pars(6) = MSd(4);
   pars(7) = MSd(5);
   pars(8) = MSu(0);
   pars(9) = MSu(1);
   pars(10) = MSu(2);
   pars(11) = MSu(3);
   pars(12) = MSu(4);
   pars(13) = MSu(5);
   pars(14) = MSe(0);
   pars(15) = MSe(1);
   pars(16) = MSe(2);
   pars(17) = MSe(3);
   pars(18) = MSe(4);
   pars(19) = MSe(5);
   pars(20) = MSv(0);
   pars(21) = MSv(1);
   pars(22) = MSv(2);
   pars(23) = MSv(3);
   pars(24) = MSv(4);
   pars(25) = MSv(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = MAh(0);
   pars(29) = MAh(1);
   pars(30) = MHpm(0);
   pars(31) = MHpm(1);
   pars(32) = MChi(0);
   pars(33) = MChi(1);
   pars(34) = MChi(2);
   pars(35) = MChi(3);
   pars(36) = MFv(0);
   pars(37) = MFv(1);
   pars(38) = MFv(2);
   pars(39) = MFv(3);
   pars(40) = MFv(4);
   pars(41) = MFv(5);
   pars(42) = MCha(0);
   pars(43) = MCha(1);
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

std::string CLASSNAME::name() const
{
   return "MSSMRHN";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   MSSMRHN_soft_parameters::run_to(scale, eps);
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
   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)
      ) + AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) +
      Conj(Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) +
      Conj(Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) -
      0.7071067811865475*vu*Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) -
      0.7071067811865475*vu*Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) -
      0.7071067811865475*vu*Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) +
      Conj(Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) -
      0.7071067811865475*vu*Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) -
      0.7071067811865475*vu*Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) -
      0.7071067811865475*vu*Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) -
      0.7071067811865475*vu*Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) -
      0.7071067811865475*vu*Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) -
      0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)
      ) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Sd, eigenvalue_error > precision
      * Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif


   if (MSd.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Sd);
   }

   MSd = AbsSqrt(MSd);
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
      0.7071067811865475*vd*Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) -
      0.7071067811865475*vd*Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) -
      0.7071067811865475*vd*Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) +
      Conj(Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) -
      0.7071067811865475*vd*Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) -
      0.7071067811865475*vd*Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) -
      0.7071067811865475*vd*Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) -
      0.7071067811865475*vd*Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) -
      0.7071067811865475*vd*Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) -
      0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      0,0)) + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      1,0)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      2,0)) + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Su, eigenvalue_error > precision
      * Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif


   if (MSu.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Su);
   }

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
      0.7071067811865475*vu*Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) -
      0.7071067811865475*vu*Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) -
      0.7071067811865475*vu*Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) +
      Conj(Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) -
      0.7071067811865475*vu*Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) -
      0.7071067811865475*vu*Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) -
      0.7071067811865475*vu*Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) -
      0.7071067811865475*vu*Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) -
      0.7071067811865475*vu*Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) -
      0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)
      ) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Se, eigenvalue_error > precision
      * Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif


   if (MSe.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sv() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yv(0,0)) + AbsSqr(Yv(1,0)) + AbsSqr(Yv(2,0)))*Sqr(
      vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(0,1) +
      Conj(Yv(1,0))*Yv(1,1) + Conj(Yv(2,0))*Yv(2,1));
   mass_matrix_Sv(0,2) = ml2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(0,2) +
      Conj(Yv(1,0))*Yv(1,2) + Conj(Yv(2,0))*Yv(2,2));
   mass_matrix_Sv(0,3) = 0.7071067811865475*vu*Conj(TYv(0,0)) -
      0.7071067811865475*vd*Conj(Yv(0,0))*Mu;
   mass_matrix_Sv(0,4) = 0.7071067811865475*vu*Conj(TYv(1,0)) -
      0.7071067811865475*vd*Conj(Yv(1,0))*Mu;
   mass_matrix_Sv(0,5) = 0.7071067811865475*vu*Conj(TYv(2,0)) -
      0.7071067811865475*vd*Conj(Yv(2,0))*Mu;
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yv(0,1)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(2,1)))*Sqr(
      vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(0,1))*Yv(0,2) +
      Conj(Yv(1,1))*Yv(1,2) + Conj(Yv(2,1))*Yv(2,2));
   mass_matrix_Sv(1,3) = 0.7071067811865475*vu*Conj(TYv(0,1)) -
      0.7071067811865475*vd*Conj(Yv(0,1))*Mu;
   mass_matrix_Sv(1,4) = 0.7071067811865475*vu*Conj(TYv(1,1)) -
      0.7071067811865475*vd*Conj(Yv(1,1))*Mu;
   mass_matrix_Sv(1,5) = 0.7071067811865475*vu*Conj(TYv(2,1)) -
      0.7071067811865475*vd*Conj(Yv(2,1))*Mu;
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yv(0,2)) + AbsSqr(Yv(1,2)) + AbsSqr(Yv(2,2)))*Sqr(
      vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(2,3) = 0.7071067811865475*vu*Conj(TYv(0,2)) -
      0.7071067811865475*vd*Conj(Yv(0,2))*Mu;
   mass_matrix_Sv(2,4) = 0.7071067811865475*vu*Conj(TYv(1,2)) -
      0.7071067811865475*vd*Conj(Yv(1,2))*Mu;
   mass_matrix_Sv(2,5) = 0.7071067811865475*vu*Conj(TYv(2,2)) -
      0.7071067811865475*vd*Conj(Yv(2,2))*Mu;
   mass_matrix_Sv(3,3) = 0.25*(AbsSqr(Mv(0,0)) + AbsSqr(Mv(0,1)) + AbsSqr
      (Mv(0,2))) + 0.25*(AbsSqr(Mv(0,0)) + AbsSqr(Mv(1,0)) + AbsSqr(Mv(2,0))) +
      0.25*(AbsSqr(Mv(0,0)) + Conj(Mv(1,0))*Mv(0,1) + Conj(Mv(2,0))*Mv(0,2)) +
      0.25*(AbsSqr(Mv(0,0)) + Conj(Mv(0,1))*Mv(1,0) + Conj(Mv(0,2))*Mv(2,0)) +
      mv2(0,0) + 0.5*(AbsSqr(Yv(0,0)) + AbsSqr(Yv(0,1)) + AbsSqr(Yv(0,2)))*Sqr
      (vu);
   mass_matrix_Sv(3,4) = 0.25*(Conj(Mv(1,0))*Mv(0,0) + Conj(Mv(1,1))*Mv(0
      ,1) + Conj(Mv(1,2))*Mv(0,2)) + 0.25*(Conj(Mv(0,1))*Mv(0,0) + Conj(Mv(1,1)
      )*Mv(0,1) + Conj(Mv(2,1))*Mv(0,2)) + 0.25*(Conj(Mv(1,0))*Mv(0,0) + Conj(
      Mv(1,1))*Mv(1,0) + Conj(Mv(1,2))*Mv(2,0)) + 0.25*(Conj(Mv(0,1))*Mv(0,0) +
      Conj(Mv(1,1))*Mv(1,0) + Conj(Mv(2,1))*Mv(2,0)) + mv2(0,1) + 0.5*Sqr(vu)*
      (Conj(Yv(1,0))*Yv(0,0) + Conj(Yv(1,1))*Yv(0,1) + Conj(Yv(1,2))*Yv(0,2));
   mass_matrix_Sv(3,5) = 0.25*(Conj(Mv(0,2))*Mv(0,0) + Conj(Mv(1,2))*Mv(0
      ,1) + Conj(Mv(2,2))*Mv(0,2)) + 0.25*(Conj(Mv(2,0))*Mv(0,0) + Conj(Mv(2,1)
      )*Mv(0,1) + Conj(Mv(2,2))*Mv(0,2)) + 0.25*(Conj(Mv(0,2))*Mv(0,0) + Conj(
      Mv(1,2))*Mv(1,0) + Conj(Mv(2,2))*Mv(2,0)) + 0.25*(Conj(Mv(2,0))*Mv(0,0) +
      Conj(Mv(2,1))*Mv(1,0) + Conj(Mv(2,2))*Mv(2,0)) + mv2(0,2) + 0.5*Sqr(vu)*
      (Conj(Yv(2,0))*Yv(0,0) + Conj(Yv(2,1))*Yv(0,1) + Conj(Yv(2,2))*Yv(0,2));
   mass_matrix_Sv(4,4) = 0.25*(AbsSqr(Mv(1,0)) + AbsSqr(Mv(1,1)) + AbsSqr
      (Mv(1,2))) + 0.25*(AbsSqr(Mv(0,1)) + AbsSqr(Mv(1,1)) + AbsSqr(Mv(2,1))) +
      0.25*(AbsSqr(Mv(1,1)) + Conj(Mv(0,1))*Mv(1,0) + Conj(Mv(2,1))*Mv(1,2)) +
      0.25*(AbsSqr(Mv(1,1)) + Conj(Mv(1,0))*Mv(0,1) + Conj(Mv(1,2))*Mv(2,1)) +
      mv2(1,1) + 0.5*(AbsSqr(Yv(1,0)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(1,2)))*Sqr
      (vu);
   mass_matrix_Sv(4,5) = 0.25*(Conj(Mv(0,2))*Mv(1,0) + Conj(Mv(1,2))*Mv(1
      ,1) + Conj(Mv(2,2))*Mv(1,2)) + 0.25*(Conj(Mv(2,0))*Mv(1,0) + Conj(Mv(2,1)
      )*Mv(1,1) + Conj(Mv(2,2))*Mv(1,2)) + 0.25*(Conj(Mv(0,2))*Mv(0,1) + Conj(
      Mv(1,2))*Mv(1,1) + Conj(Mv(2,2))*Mv(2,1)) + 0.25*(Conj(Mv(2,0))*Mv(0,1) +
      Conj(Mv(2,1))*Mv(1,1) + Conj(Mv(2,2))*Mv(2,1)) + mv2(1,2) + 0.5*Sqr(vu)*
      (Conj(Yv(2,0))*Yv(1,0) + Conj(Yv(2,1))*Yv(1,1) + Conj(Yv(2,2))*Yv(1,2));
   mass_matrix_Sv(5,5) = 0.25*(AbsSqr(Mv(0,2)) + AbsSqr(Mv(1,2)) + AbsSqr
      (Mv(2,2))) + 0.25*(AbsSqr(Mv(2,0)) + AbsSqr(Mv(2,1)) + AbsSqr(Mv(2,2))) +
      0.25*(AbsSqr(Mv(2,2)) + Conj(Mv(2,0))*Mv(0,2) + Conj(Mv(2,1))*Mv(1,2)) +
      0.25*(AbsSqr(Mv(2,2)) + Conj(Mv(0,2))*Mv(2,0) + Conj(Mv(1,2))*Mv(2,1)) +
      mv2(2,2) + 0.5*(AbsSqr(Yv(2,0)) + AbsSqr(Yv(2,1)) + AbsSqr(Yv(2,2)))*Sqr
      (vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Sv, eigenvalue_error > precision
      * Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif


   if (MSv.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + AbsSqr(Mu) + 0.225*Sqr(g1)*Sqr(vd) +
      0.375*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*Conj(BMu) - 0.15*vd*vu*Sqr(g1) -
      0.25*vd*vu*Sqr(g2);
   mass_matrix_hh(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.225*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::hh, eigenvalue_error > precision
      * Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif


   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + AbsSqr(Mu) + 0.3872983346207417*g1*g2*Cos
      (ThetaW())*Sin(ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*
      Sqr(vd)*Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) - 0.3872983346207417*g1*
      g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW(
      ))) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW
      ())*Sqr(vu) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2
      )*Sqr(vu)*Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Ah, eigenvalue_error > precision
      * Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif


   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + AbsSqr(Mu) + 0.075*Sqr(g1)*Sqr(vd) +
      0.375*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = Conj(BMu);
   mass_matrix_Hpm(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) +
      0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Hpm, eigenvalue_error > precision
      * Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif


   if (MHpm.minCoeff() < 0.) {
      problems.flag_tachyon(MSSMRHN_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
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
   problems.flag_bad_mass(MSSMRHN_info::Chi, eigenvalue_error > precision
      * Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif

}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(0,3) = 0.7071067811865475*vu*Yv(0,0);
   mass_matrix_Fv(0,4) = 0.7071067811865475*vu*Yv(1,0);
   mass_matrix_Fv(0,5) = 0.7071067811865475*vu*Yv(2,0);
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(1,3) = 0.7071067811865475*vu*Yv(0,1);
   mass_matrix_Fv(1,4) = 0.7071067811865475*vu*Yv(1,1);
   mass_matrix_Fv(1,5) = 0.7071067811865475*vu*Yv(2,1);
   mass_matrix_Fv(2,2) = 0;
   mass_matrix_Fv(2,3) = 0.7071067811865475*vu*Yv(0,2);
   mass_matrix_Fv(2,4) = 0.7071067811865475*vu*Yv(1,2);
   mass_matrix_Fv(2,5) = 0.7071067811865475*vu*Yv(2,2);
   mass_matrix_Fv(3,3) = Mv(0,0);
   mass_matrix_Fv(3,4) = 0.5*Mv(0,1) + 0.5*Mv(1,0);
   mass_matrix_Fv(3,5) = 0.5*Mv(0,2) + 0.5*Mv(2,0);
   mass_matrix_Fv(4,4) = Mv(1,1);
   mass_matrix_Fv(4,5) = 0.5*Mv(1,2) + 0.5*Mv(2,1);
   mass_matrix_Fv(5,5) = Mv(2,2);

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   const auto mass_matrix_Fv(get_mass_matrix_Fv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Fv, eigenvalue_error > precision
      * Abs(MFv(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(MSSMRHN_info::Cha, eigenvalue_error > precision
      * Abs(MCha(0)));
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
   problems.flag_bad_mass(MSSMRHN_info::Fe, eigenvalue_error > precision
      * Abs(MFe(0)));
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
   problems.flag_bad_mass(MSSMRHN_info::Fd, eigenvalue_error > precision
      * Abs(MFd(0)));
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
   problems.flag_bad_mass(MSSMRHN_info::Fu, eigenvalue_error > precision
      * Abs(MFu(0)));
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
      problems.flag_tachyon(MSSMRHN_info::VWm);
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
   double result = Re(mHd2*vd + vd*AbsSqr(Mu) - 0.5*vu*BMu - 0.5*vu*Conj(BMu) +
      0.075*Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) - 0.075*vd*Sqr(g1)*
      Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = Re(mHu2*vu + vu*AbsSqr(Mu) - 0.5*vd*BMu - 0.5*vd*Conj(BMu) +
      0.075*Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) - 0.075*vu*Sqr(g1)*
      Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd));

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

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_14;
   std::complex<double> tmp_15;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_15 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_14 += tmp_15;
   result += (0.1*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_14;
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
   std::complex<double> tmp_22;
   std::complex<double> tmp_23;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_23 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_22 += tmp_23;
   result += (-0.1*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_22;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_24;
      std::complex<double> tmp_25;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_25 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_24 += tmp_25;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_24;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_26;
   std::complex<double> tmp_27;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_27 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_26 += tmp_27;
   result += (0.1*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_26;
   std::complex<double> tmp_28;
   std::complex<double> tmp_29;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_30;
      std::complex<double> tmp_31;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_32;
         std::complex<double> tmp_33;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_33 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_32 += tmp_33;
         tmp_31 += (KroneckerDelta(gO2,3 + j2)) * tmp_32;
      }
      tmp_30 += tmp_31;
      tmp_29 += (KroneckerDelta(gO1,3 + j3)) * tmp_30;
   }
   tmp_28 += tmp_29;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_28;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_34;
      std::complex<double> tmp_35;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_35 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_34 += tmp_35;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_34;
   }
   std::complex<double> tmp_36;
   std::complex<double> tmp_37;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_37 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_36 += tmp_37;
   result += (-0.1*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_36;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_38;
      std::complex<double> tmp_39;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_39 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_38 += tmp_39;
      result += (UP(gI2,1)) * tmp_38;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_40;
   std::complex<double> tmp_41;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_42;
      std::complex<double> tmp_43;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_43 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_42 += tmp_43;
      tmp_41 += (Conj(ZUL(gI1,j2))) * tmp_42;
   }
   tmp_40 += tmp_41;
   result += (Conj(UM(gI2,1))) * tmp_40;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_44;
   std::complex<double> tmp_45;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_45 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_44 += tmp_45;
   result += (-0.3651483716701107*g1*ZN(gI2,0)) * tmp_44;
   if (gO2 < 3) {
      std::complex<double> tmp_46;
      std::complex<double> tmp_47;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_47 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_46 += tmp_47;
      result += (-ZN(gI2,2)) * tmp_46;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_48;
   std::complex<double> tmp_49;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_50;
      std::complex<double> tmp_51;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_51 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_50 += tmp_51;
      tmp_49 += (Conj(ZDL(gI1,j2))) * tmp_50;
   }
   tmp_48 += tmp_49;
   result += (-Conj(ZN(gI2,2))) * tmp_48;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_52;
   std::complex<double> tmp_54;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_54 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_52 += tmp_54;
   std::complex<double> tmp_53;
   std::complex<double> tmp_55;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_55 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_53 += tmp_55;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_52 * tmp_53;
   std::complex<double> tmp_56;
   std::complex<double> tmp_58;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_58 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_56 += tmp_58;
   std::complex<double> tmp_57;
   std::complex<double> tmp_59;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_59 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_57 += tmp_59;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_56 * tmp_57;
   std::complex<double> tmp_60;
   std::complex<double> tmp_62;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_62 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_60 += tmp_62;
   std::complex<double> tmp_61;
   std::complex<double> tmp_63;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_63 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_61 += tmp_63;
   result += (-0.05*Sqr(g1)) * tmp_60 * tmp_61;
   std::complex<double> tmp_64;
   std::complex<double> tmp_66;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_66 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_64 += tmp_66;
   std::complex<double> tmp_65;
   std::complex<double> tmp_67;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_67 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_65 += tmp_67;
   result += (-0.1*Sqr(g1)) * tmp_64 * tmp_65;
   std::complex<double> tmp_68;
   std::complex<double> tmp_70;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_70 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_68 += tmp_70;
   std::complex<double> tmp_69;
   std::complex<double> tmp_71;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_71 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_69 += tmp_71;
   result += (-0.05*Sqr(g1)) * tmp_68 * tmp_69;
   std::complex<double> tmp_72;
   std::complex<double> tmp_74;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_74 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_72 += tmp_74;
   std::complex<double> tmp_73;
   std::complex<double> tmp_75;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_75 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_73 += tmp_75;
   result += (-0.1*Sqr(g1)) * tmp_72 * tmp_73;
   std::complex<double> tmp_76;
   std::complex<double> tmp_78;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_78 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_76 += tmp_78;
   std::complex<double> tmp_77;
   std::complex<double> tmp_79;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_79 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_77 += tmp_79;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_76 * tmp_77;
   std::complex<double> tmp_80;
   std::complex<double> tmp_82;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_82 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_80 += tmp_82;
   std::complex<double> tmp_81;
   std::complex<double> tmp_83;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_83 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_81 += tmp_83;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_80 * tmp_81;
   std::complex<double> tmp_84;
   std::complex<double> tmp_86;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_87;
      std::complex<double> tmp_88;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_88 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_87 += tmp_88;
      tmp_86 += (Conj(ZD(gI2,j2))) * tmp_87;
   }
   tmp_84 += tmp_86;
   std::complex<double> tmp_85;
   std::complex<double> tmp_89;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_90;
      std::complex<double> tmp_91;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_91 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_90 += tmp_91;
      tmp_89 += (ZD(gI1,j4)) * tmp_90;
   }
   tmp_85 += tmp_89;
   result += (-1) * tmp_84 * tmp_85;
   if (gO1 < 3) {
      std::complex<double> tmp_92;
      std::complex<double> tmp_93;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_93 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_92 += tmp_93;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_92;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_94;
      std::complex<double> tmp_95;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_95 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_94 += tmp_95;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_94;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_96;
      std::complex<double> tmp_97;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_97 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_96 += tmp_97;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_96;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_98;
      std::complex<double> tmp_99;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_99 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_98 += tmp_99;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_98;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_100;
      std::complex<double> tmp_101;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_101 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_100 += tmp_101;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_100;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_102;
      std::complex<double> tmp_103;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_103 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_102 += tmp_103;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_102;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_104;
      std::complex<double> tmp_106;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_106 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_104 += tmp_106;
      std::complex<double> tmp_105;
      std::complex<double> tmp_107;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_108;
         std::complex<double> tmp_109;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_109 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_108 += tmp_109;
         tmp_107 += (ZD(gI1,j4)) * tmp_108;
      }
      tmp_105 += tmp_107;
      result += (-3) * tmp_104 * tmp_105;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_111 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_110 += tmp_111;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_110;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_113 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_112 += tmp_113;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_112;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_114;
      std::complex<double> tmp_115;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_115 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_114 += tmp_115;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_114;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_117 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_116 += tmp_117;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_116;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_118;
      std::complex<double> tmp_120;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_121;
         std::complex<double> tmp_122;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_122 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_121 += tmp_122;
         tmp_120 += (Conj(ZD(gI2,j2))) * tmp_121;
      }
      tmp_118 += tmp_120;
      std::complex<double> tmp_119;
      std::complex<double> tmp_123;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_123 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_119 += tmp_123;
      result += (-3) * tmp_118 * tmp_119;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_126;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_126 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_124 += tmp_126;
      std::complex<double> tmp_125;
      std::complex<double> tmp_127;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_127 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_125 += tmp_127;
      result += (-1) * tmp_124 * tmp_125;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_129;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_129 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_128 += tmp_129;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_128;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_130;
      std::complex<double> tmp_131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_131 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_130 += tmp_131;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_130;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_133 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_132 += tmp_133;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_132;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_135 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_134 += tmp_135;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_134;
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

   std::complex<double> tmp_136;
   std::complex<double> tmp_138;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_138 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_136 += tmp_138;
   std::complex<double> tmp_137;
   std::complex<double> tmp_139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_139 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_137 += tmp_139;
   result += (0.05*Sqr(g1)) * tmp_136 * tmp_137;
   std::complex<double> tmp_140;
   std::complex<double> tmp_142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_142 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_140 += tmp_142;
   std::complex<double> tmp_141;
   std::complex<double> tmp_143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_143 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_141 += tmp_143;
   result += (-0.1*Sqr(g1)) * tmp_140 * tmp_141;
   std::complex<double> tmp_144;
   std::complex<double> tmp_146;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_146 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_144 += tmp_146;
   std::complex<double> tmp_145;
   std::complex<double> tmp_147;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_147 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_145 += tmp_147;
   result += (0.05*Sqr(g1)) * tmp_144 * tmp_145;
   std::complex<double> tmp_148;
   std::complex<double> tmp_150;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_150 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_148 += tmp_150;
   std::complex<double> tmp_149;
   std::complex<double> tmp_151;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_151 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_149 += tmp_151;
   result += (-0.1*Sqr(g1)) * tmp_148 * tmp_149;
   if (gO1 < 3) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_153 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_152 += tmp_153;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_152;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_154;
      std::complex<double> tmp_155;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_155 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_154 += tmp_155;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_154;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_156;
      std::complex<double> tmp_157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_157 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_156 += tmp_157;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_156;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_158;
      std::complex<double> tmp_159;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_159 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_158 += tmp_159;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_158;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_160;
      std::complex<double> tmp_161;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_161 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_160 += tmp_161;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_160;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_162;
      std::complex<double> tmp_163;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_163 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_162 += tmp_163;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_162;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_164;
      std::complex<double> tmp_166;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_166 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_164 += tmp_166;
      std::complex<double> tmp_165;
      std::complex<double> tmp_167;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_168;
         std::complex<double> tmp_169;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_169 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_168 += tmp_169;
         tmp_167 += (ZE(gI1,j4)) * tmp_168;
      }
      tmp_165 += tmp_167;
      result += (-1) * tmp_164 * tmp_165;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_172;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_173;
         std::complex<double> tmp_174;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_174 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_173 += tmp_174;
         tmp_172 += (Conj(ZE(gI2,j2))) * tmp_173;
      }
      tmp_170 += tmp_172;
      std::complex<double> tmp_171;
      std::complex<double> tmp_175;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_175 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_171 += tmp_175;
      result += (-1) * tmp_170 * tmp_171;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_176;
   std::complex<double> tmp_178;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_178 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_176 += tmp_178;
   std::complex<double> tmp_177;
   std::complex<double> tmp_179;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_179 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_177 += tmp_179;
   result += (-0.05*Sqr(g1)) * tmp_176 * tmp_177;
   std::complex<double> tmp_180;
   std::complex<double> tmp_182;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_182 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_180 += tmp_182;
   std::complex<double> tmp_181;
   std::complex<double> tmp_183;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_183 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_181 += tmp_183;
   result += (0.2*Sqr(g1)) * tmp_180 * tmp_181;
   std::complex<double> tmp_184;
   std::complex<double> tmp_186;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_186 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_184 += tmp_186;
   std::complex<double> tmp_185;
   std::complex<double> tmp_187;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_187 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_185 += tmp_187;
   result += (-0.05*Sqr(g1)) * tmp_184 * tmp_185;
   std::complex<double> tmp_188;
   std::complex<double> tmp_190;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_190 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_188 += tmp_190;
   std::complex<double> tmp_189;
   std::complex<double> tmp_191;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_191 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_189 += tmp_191;
   result += (0.2*Sqr(g1)) * tmp_188 * tmp_189;
   std::complex<double> tmp_192;
   std::complex<double> tmp_194;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_195;
      std::complex<double> tmp_196;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_196 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_195 += tmp_196;
      tmp_194 += (Conj(ZU(gI2,j2))) * tmp_195;
   }
   tmp_192 += tmp_194;
   std::complex<double> tmp_193;
   std::complex<double> tmp_197;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_198;
      std::complex<double> tmp_199;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_199 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_198 += tmp_199;
      tmp_197 += (ZU(gI1,j4)) * tmp_198;
   }
   tmp_193 += tmp_197;
   result += (-1) * tmp_192 * tmp_193;
   if (gO1 < 3) {
      std::complex<double> tmp_200;
      std::complex<double> tmp_201;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_201 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_200 += tmp_201;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_200;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_202;
      std::complex<double> tmp_203;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_203 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_202 += tmp_203;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_202;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_204;
      std::complex<double> tmp_205;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_205 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_204 += tmp_205;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_204;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_206;
      std::complex<double> tmp_207;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_207 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_206 += tmp_207;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_206;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_208;
      std::complex<double> tmp_209;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_209 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_208 += tmp_209;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_208;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_210;
      std::complex<double> tmp_211;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_211 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_210 += tmp_211;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_210;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_212;
      std::complex<double> tmp_214;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_214 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_212 += tmp_214;
      std::complex<double> tmp_213;
      std::complex<double> tmp_215;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_215 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_213 += tmp_215;
      result += (-1) * tmp_212 * tmp_213;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_216;
   std::complex<double> tmp_218;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_218 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_216 += tmp_218;
   std::complex<double> tmp_217;
   std::complex<double> tmp_219;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_219 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_217 += tmp_219;
   result += (0.05*Sqr(g1)) * tmp_216 * tmp_217;
   std::complex<double> tmp_220;
   std::complex<double> tmp_222;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_222 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_220 += tmp_222;
   std::complex<double> tmp_221;
   std::complex<double> tmp_223;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_223 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_221 += tmp_223;
   result += (0.05*Sqr(g1)) * tmp_220 * tmp_221;
   if (gO1 < 3) {
      std::complex<double> tmp_224;
      std::complex<double> tmp_225;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_225 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_224 += tmp_225;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_224;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_226;
      std::complex<double> tmp_227;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_227 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_226 += tmp_227;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_226;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_228;
      std::complex<double> tmp_229;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_229 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_228 += tmp_229;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_228;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_230;
      std::complex<double> tmp_231;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_231 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_230 += tmp_231;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_230;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_232;
   std::complex<double> tmp_233;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_234;
      std::complex<double> tmp_235;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_235 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_234 += tmp_235;
      tmp_233 += (Conj(ZD(gI1,j2))) * tmp_234;
   }
   tmp_232 += tmp_233;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_232
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_236;
      std::complex<double> tmp_237;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_237 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_236 += tmp_237;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_236;
   }
   std::complex<double> tmp_238;
   std::complex<double> tmp_239;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_240;
      std::complex<double> tmp_241;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_241 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_240 += tmp_241;
      tmp_239 += (Conj(ZD(gI1,j2))) * tmp_240;
   }
   tmp_238 += tmp_239;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1))
      * tmp_238;
   if (gO2 < 3) {
      std::complex<double> tmp_242;
      std::complex<double> tmp_243;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_243 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_242 += tmp_243;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_242;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_244;
   std::complex<double> tmp_245;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_245 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_244 += tmp_245;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_244;
   std::complex<double> tmp_246;
   std::complex<double> tmp_247;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_248;
      std::complex<double> tmp_249;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_249 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_248 += tmp_249;
      tmp_247 += (Conj(ZD(gI1,j2))) * tmp_248;
   }
   tmp_246 += tmp_247;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_246;
   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_254;
         std::complex<double> tmp_255;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_255 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_254 += tmp_255;
         tmp_253 += (KroneckerDelta(gO2,3 + j2)) * tmp_254;
      }
      tmp_252 += tmp_253;
      tmp_251 += (Conj(ZD(gI1,3 + j3))) * tmp_252;
   }
   tmp_250 += tmp_251;
   result += (-(vd*ZH(gI2,0))) * tmp_250;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_257 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_256 += tmp_257;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_256;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_258;
      std::complex<double> tmp_259;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_260;
         std::complex<double> tmp_261;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_261 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_260 += tmp_261;
         tmp_259 += (Conj(ZD(gI1,j2))) * tmp_260;
      }
      tmp_258 += tmp_259;
      result += (-(vd*ZH(gI2,0))) * tmp_258;
   }
   std::complex<double> tmp_262;
   std::complex<double> tmp_263;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_263 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_262 += tmp_263;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_262;
   std::complex<double> tmp_264;
   std::complex<double> tmp_265;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_266;
      std::complex<double> tmp_267;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_267 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_266 += tmp_267;
      tmp_265 += (Conj(ZD(gI1,j2))) * tmp_266;
   }
   tmp_264 += tmp_265;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_264;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_269 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_268 += tmp_269;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_268;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_270;
   std::complex<double> tmp_271;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_272;
      std::complex<double> tmp_273;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_273 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_272 += tmp_273;
      tmp_271 += (Conj(ZU(gI1,j2))) * tmp_272;
   }
   tmp_270 += tmp_271;
   result += (ZP(gI2,0)) * tmp_270;
   std::complex<double> tmp_274;
   std::complex<double> tmp_275;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_276;
      std::complex<double> tmp_277;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_278;
         std::complex<double> tmp_279;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_279 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_278 += tmp_279;
         tmp_277 += (KroneckerDelta(gO2,3 + j2)) * tmp_278;
      }
      tmp_276 += tmp_277;
      tmp_275 += (Conj(ZU(gI1,3 + j3))) * tmp_276;
   }
   tmp_274 += tmp_275;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_274;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_280;
      std::complex<double> tmp_281;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_281 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_280 += tmp_281;
      result += (Mu*ZP(gI2,0)) * tmp_280;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_282;
      std::complex<double> tmp_283;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_284;
         std::complex<double> tmp_285;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_285 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_284 += tmp_285;
         tmp_283 += (Conj(ZU(gI1,j2))) * tmp_284;
      }
      tmp_282 += tmp_283;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_282;
   }
   std::complex<double> tmp_286;
   std::complex<double> tmp_287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_288;
      std::complex<double> tmp_289;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_289 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_288 += tmp_289;
      tmp_287 += (Conj(ZU(gI1,j2))) * tmp_288;
   }
   tmp_286 += tmp_287;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_286;
   std::complex<double> tmp_290;
   std::complex<double> tmp_291;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_292;
      std::complex<double> tmp_293;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_294;
         std::complex<double> tmp_295;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_295 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_294 += tmp_295;
         tmp_293 += (KroneckerDelta(gO2,3 + j2)) * tmp_294;
      }
      tmp_292 += tmp_293;
      tmp_291 += (Conj(ZU(gI1,3 + j3))) * tmp_292;
   }
   tmp_290 += tmp_291;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_290;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_296;
      std::complex<double> tmp_297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_297 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_296 += tmp_297;
      result += (ZP(gI2,1)) * tmp_296;
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

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_302;
   std::complex<double> tmp_303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_303 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_302 += tmp_303;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_302;

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

   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_305 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_304 += tmp_305;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_304;
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

   std::complex<double> tmp_306;
   std::complex<double> tmp_307;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_307 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_306 += tmp_307;
   result += (0.2581988897471611*g1*Sin(ThetaW())) * tmp_306;
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

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_308;
   std::complex<double> tmp_309;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_309 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_308 += tmp_309;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_308;
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

   std::complex<double> tmp_310;
   std::complex<double> tmp_311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_311 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_310 += tmp_311;
   result += (-0.2*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_310;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_313 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_312 += tmp_313;
   result += (0.2*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_312;
   std::complex<double> tmp_314;
   std::complex<double> tmp_315;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_316;
      std::complex<double> tmp_317;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_318;
         std::complex<double> tmp_319;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_319 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_318 += tmp_319;
         tmp_317 += (KroneckerDelta(gO2,3 + j2)) * tmp_318;
      }
      tmp_316 += tmp_317;
      tmp_315 += (KroneckerDelta(gO1,3 + j3)) * tmp_316;
   }
   tmp_314 += tmp_315;
   result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_314;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_320;
      std::complex<double> tmp_321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_321 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_320 += tmp_321;
      result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_320;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_322;
   std::complex<double> tmp_323;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_323 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_322 += tmp_323;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_322;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_324 += tmp_325;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_324;
   }
   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_327 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_326 += tmp_327;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_326;
   std::complex<double> tmp_328;
   std::complex<double> tmp_329;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_330;
      std::complex<double> tmp_331;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_332;
         std::complex<double> tmp_333;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_333 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_332 += tmp_333;
         tmp_331 += (KroneckerDelta(gO2,3 + j2)) * tmp_332;
      }
      tmp_330 += tmp_331;
      tmp_329 += (KroneckerDelta(gO1,3 + j3)) * tmp_330;
   }
   tmp_328 += tmp_329;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_328;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_334;
   std::complex<double> tmp_335;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_335 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_334 += tmp_335;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_334;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   std::complex<double> tmp_336;
   std::complex<double> tmp_337;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_337 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_336 += tmp_337;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_336;
   std::complex<double> tmp_338;
   std::complex<double> tmp_339;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_340;
      std::complex<double> tmp_341;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_342;
         std::complex<double> tmp_343;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_343 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_342 += tmp_343;
         tmp_341 += (KroneckerDelta(gO2,3 + j2)) * tmp_342;
      }
      tmp_340 += tmp_341;
      tmp_339 += (KroneckerDelta(gO1,3 + j3)) * tmp_340;
   }
   tmp_338 += tmp_339;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_338;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_344;
      std::complex<double> tmp_345;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_345 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_344 += tmp_345;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_344;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_346;
      std::complex<double> tmp_347;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_347 += Conj(Yd(j1,gO2))*ZDR(gI2,j1);
      }
      tmp_346 += tmp_347;
      result += (UM(gI1,1)) * tmp_346;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_348;
   std::complex<double> tmp_349;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_350;
      std::complex<double> tmp_351;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_351 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_350 += tmp_351;
      tmp_349 += (Conj(ZDL(gI2,j2))) * tmp_350;
   }
   tmp_348 += tmp_349;
   result += (Conj(UP(gI1,1))) * tmp_348;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_354;
      std::complex<double> tmp_355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_355 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_354 += tmp_355;
      tmp_353 += (Conj(ZD(gI2,j2))) * tmp_354;
   }
   tmp_352 += tmp_353;
   result += (Conj(Mu)*ZP(gI1,0)) * tmp_352;
   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_358;
      std::complex<double> tmp_359;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_360;
         std::complex<double> tmp_361;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_361 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_360 += tmp_361;
         tmp_359 += (KroneckerDelta(gO2,3 + j2)) * tmp_360;
      }
      tmp_358 += tmp_359;
      tmp_357 += (Conj(ZD(gI2,3 + j3))) * tmp_358;
   }
   tmp_356 += tmp_357;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_356;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_362;
      std::complex<double> tmp_363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_363 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_362 += tmp_363;
      result += (ZP(gI1,0)) * tmp_362;
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
         tmp_365 += (Conj(ZD(gI2,j2))) * tmp_366;
      }
      tmp_364 += tmp_365;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_364;
   }
   std::complex<double> tmp_368;
   std::complex<double> tmp_369;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_370;
      std::complex<double> tmp_371;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_371 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_370 += tmp_371;
      tmp_369 += (Conj(ZD(gI2,j2))) * tmp_370;
   }
   tmp_368 += tmp_369;
   result += (ZP(gI1,1)) * tmp_368;
   std::complex<double> tmp_372;
   std::complex<double> tmp_373;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_374;
      std::complex<double> tmp_375;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_376;
         std::complex<double> tmp_377;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_377 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_376 += tmp_377;
         tmp_375 += (KroneckerDelta(gO2,3 + j2)) * tmp_376;
      }
      tmp_374 += tmp_375;
      tmp_373 += (Conj(ZD(gI2,3 + j3))) * tmp_374;
   }
   tmp_372 += tmp_373;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_372;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_378;
      std::complex<double> tmp_379;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_379 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_378 += tmp_379;
      result += (Mu*ZP(gI1,1)) * tmp_378;
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
         tmp_381 += (Conj(ZD(gI2,j2))) * tmp_382;
      }
      tmp_380 += tmp_381;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_380;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_384;
   std::complex<double> tmp_385;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_385 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_384 += tmp_385;
   result += (0.7302967433402214*g1*ZN(gI2,0)) * tmp_384;
   if (gO2 < 3) {
      std::complex<double> tmp_386;
      std::complex<double> tmp_387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_387 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_386 += tmp_387;
      result += (-ZN(gI2,3)) * tmp_386;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_390;
      std::complex<double> tmp_391;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_391 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_390 += tmp_391;
      tmp_389 += (Conj(ZUL(gI1,j2))) * tmp_390;
   }
   tmp_388 += tmp_389;
   result += (-Conj(ZN(gI2,3))) * tmp_388;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZN(gI2,0))*Conj(ZUL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*Conj(ZUL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_392;
   std::complex<double> tmp_394;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_394 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_392 += tmp_394;
   std::complex<double> tmp_393;
   std::complex<double> tmp_395;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_395 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_393 += tmp_395;
   result += (0.1*Sqr(g1)) * tmp_392 * tmp_393;
   std::complex<double> tmp_396;
   std::complex<double> tmp_398;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_398 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_396 += tmp_398;
   std::complex<double> tmp_397;
   std::complex<double> tmp_399;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_399 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_397 += tmp_399;
   result += (0.2*Sqr(g1)) * tmp_396 * tmp_397;
   std::complex<double> tmp_400;
   std::complex<double> tmp_402;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_402 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_400 += tmp_402;
   std::complex<double> tmp_401;
   std::complex<double> tmp_403;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_403 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_401 += tmp_403;
   result += (0.1*Sqr(g1)) * tmp_400 * tmp_401;
   std::complex<double> tmp_404;
   std::complex<double> tmp_406;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_406 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_404 += tmp_406;
   std::complex<double> tmp_405;
   std::complex<double> tmp_407;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_407 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_405 += tmp_407;
   result += (0.2*Sqr(g1)) * tmp_404 * tmp_405;
   std::complex<double> tmp_408;
   std::complex<double> tmp_410;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_411;
      std::complex<double> tmp_412;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_412 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_411 += tmp_412;
      tmp_410 += (Conj(ZD(gI2,j2))) * tmp_411;
   }
   tmp_408 += tmp_410;
   std::complex<double> tmp_409;
   std::complex<double> tmp_413;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_414;
      std::complex<double> tmp_415;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_415 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_414 += tmp_415;
      tmp_413 += (ZD(gI1,j4)) * tmp_414;
   }
   tmp_409 += tmp_413;
   result += (-1) * tmp_408 * tmp_409;
   if (gO1 < 3) {
      std::complex<double> tmp_416;
      std::complex<double> tmp_417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_417 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_416 += tmp_417;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_416;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_418;
      std::complex<double> tmp_419;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_419 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_418 += tmp_419;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_418;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_420;
      std::complex<double> tmp_421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_421 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_420 += tmp_421;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_420;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_422;
      std::complex<double> tmp_423;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_423 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_422 += tmp_423;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_422;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_424;
      std::complex<double> tmp_425;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_425 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_424 += tmp_425;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_424;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_426;
      std::complex<double> tmp_427;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_427 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_426 += tmp_427;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_426;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_430;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_430 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_428 += tmp_430;
      std::complex<double> tmp_429;
      std::complex<double> tmp_431;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_431 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_429 += tmp_431;
      result += (-1) * tmp_428 * tmp_429;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_432;
   std::complex<double> tmp_434;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_434 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_432 += tmp_434;
   std::complex<double> tmp_433;
   std::complex<double> tmp_435;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_435 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_433 += tmp_435;
   result += (-0.1*Sqr(g1)) * tmp_432 * tmp_433;
   std::complex<double> tmp_436;
   std::complex<double> tmp_438;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_438 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_436 += tmp_438;
   std::complex<double> tmp_437;
   std::complex<double> tmp_439;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_439 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_437 += tmp_439;
   result += (0.2*Sqr(g1)) * tmp_436 * tmp_437;
   std::complex<double> tmp_440;
   std::complex<double> tmp_442;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_442 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_440 += tmp_442;
   std::complex<double> tmp_441;
   std::complex<double> tmp_443;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_443 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_441 += tmp_443;
   result += (-0.1*Sqr(g1)) * tmp_440 * tmp_441;
   std::complex<double> tmp_444;
   std::complex<double> tmp_446;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_446 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_444 += tmp_446;
   std::complex<double> tmp_445;
   std::complex<double> tmp_447;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_447 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_445 += tmp_447;
   result += (0.2*Sqr(g1)) * tmp_444 * tmp_445;
   if (gO1 < 3) {
      std::complex<double> tmp_448;
      std::complex<double> tmp_449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_449 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_448 += tmp_449;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_448;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_450;
      std::complex<double> tmp_451;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_451 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_450 += tmp_451;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_450;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_452;
      std::complex<double> tmp_453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_453 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_452 += tmp_453;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_452;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_454;
      std::complex<double> tmp_455;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_455 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_454 += tmp_455;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_454;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_456;
      std::complex<double> tmp_457;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_457 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_456 += tmp_457;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_456;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_458;
      std::complex<double> tmp_459;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_459 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_458 += tmp_459;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_458;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_460;
   std::complex<double> tmp_462;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_462 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_460 += tmp_462;
   std::complex<double> tmp_461;
   std::complex<double> tmp_463;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_463 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_461 += tmp_463;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_460 * tmp_461;
   std::complex<double> tmp_464;
   std::complex<double> tmp_466;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_466 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_464 += tmp_466;
   std::complex<double> tmp_465;
   std::complex<double> tmp_467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_467 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_465 += tmp_467;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_464 * tmp_465;
   std::complex<double> tmp_468;
   std::complex<double> tmp_470;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_470 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_468 += tmp_470;
   std::complex<double> tmp_469;
   std::complex<double> tmp_471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_471 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_469 += tmp_471;
   result += (0.1*Sqr(g1)) * tmp_468 * tmp_469;
   std::complex<double> tmp_472;
   std::complex<double> tmp_474;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_474 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_472 += tmp_474;
   std::complex<double> tmp_473;
   std::complex<double> tmp_475;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_475 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_473 += tmp_475;
   result += (-0.4*Sqr(g1)) * tmp_472 * tmp_473;
   std::complex<double> tmp_476;
   std::complex<double> tmp_478;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_478 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_476 += tmp_478;
   std::complex<double> tmp_477;
   std::complex<double> tmp_479;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_479 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_477 += tmp_479;
   result += (0.1*Sqr(g1)) * tmp_476 * tmp_477;
   std::complex<double> tmp_480;
   std::complex<double> tmp_482;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_482 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_480 += tmp_482;
   std::complex<double> tmp_481;
   std::complex<double> tmp_483;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_483 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_481 += tmp_483;
   result += (-0.4*Sqr(g1)) * tmp_480 * tmp_481;
   std::complex<double> tmp_484;
   std::complex<double> tmp_486;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_486 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_484 += tmp_486;
   std::complex<double> tmp_485;
   std::complex<double> tmp_487;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_487 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_485 += tmp_487;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_484 * tmp_485;
   std::complex<double> tmp_488;
   std::complex<double> tmp_490;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_490 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_488 += tmp_490;
   std::complex<double> tmp_489;
   std::complex<double> tmp_491;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_491 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_489 += tmp_491;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_488 * tmp_489;
   std::complex<double> tmp_492;
   std::complex<double> tmp_494;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_495;
      std::complex<double> tmp_496;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_496 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_495 += tmp_496;
      tmp_494 += (Conj(ZU(gI2,j2))) * tmp_495;
   }
   tmp_492 += tmp_494;
   std::complex<double> tmp_493;
   std::complex<double> tmp_497;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_498;
      std::complex<double> tmp_499;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_499 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_498 += tmp_499;
      tmp_497 += (ZU(gI1,j4)) * tmp_498;
   }
   tmp_493 += tmp_497;
   result += (-1) * tmp_492 * tmp_493;
   if (gO1 < 3) {
      std::complex<double> tmp_500;
      std::complex<double> tmp_501;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_501 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_500 += tmp_501;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_500;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_502;
      std::complex<double> tmp_503;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_503 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_502 += tmp_503;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_502;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_504;
      std::complex<double> tmp_505;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_505 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_504 += tmp_505;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_504;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_506;
      std::complex<double> tmp_507;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_507 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_506 += tmp_507;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_506;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_508;
      std::complex<double> tmp_509;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_509 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_508 += tmp_509;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_508;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_510;
      std::complex<double> tmp_511;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_511 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_510 += tmp_511;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_510;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_512;
      std::complex<double> tmp_514;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_514 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_512 += tmp_514;
      std::complex<double> tmp_513;
      std::complex<double> tmp_515;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_516;
         std::complex<double> tmp_517;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_517 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_516 += tmp_517;
         tmp_515 += (ZU(gI1,j4)) * tmp_516;
      }
      tmp_513 += tmp_515;
      result += (-3) * tmp_512 * tmp_513;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_518;
      std::complex<double> tmp_519;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_519 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_518 += tmp_519;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_518;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_520;
      std::complex<double> tmp_521;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_521 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_520 += tmp_521;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_520;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_522;
      std::complex<double> tmp_523;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_523 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_522 += tmp_523;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_522;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_525 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_524 += tmp_525;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_524;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_526;
      std::complex<double> tmp_528;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_529;
         std::complex<double> tmp_530;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_530 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_529 += tmp_530;
         tmp_528 += (Conj(ZU(gI2,j2))) * tmp_529;
      }
      tmp_526 += tmp_528;
      std::complex<double> tmp_527;
      std::complex<double> tmp_531;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_531 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_527 += tmp_531;
      result += (-3) * tmp_526 * tmp_527;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_532;
      std::complex<double> tmp_534;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_534 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_532 += tmp_534;
      std::complex<double> tmp_533;
      std::complex<double> tmp_535;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_535 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_533 += tmp_535;
      result += (-1) * tmp_532 * tmp_533;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_536;
      std::complex<double> tmp_537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_537 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_536 += tmp_537;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_536;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_538;
      std::complex<double> tmp_539;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_539 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_538 += tmp_539;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_538;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_540;
      std::complex<double> tmp_541;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_541 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_540 += tmp_541;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_540;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_542;
      std::complex<double> tmp_543;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_543 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_542 += tmp_543;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_542;
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

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_544;
   std::complex<double> tmp_546;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_546 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_544 += tmp_546;
   std::complex<double> tmp_545;
   std::complex<double> tmp_547;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_547 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_545 += tmp_547;
   result += (-0.1*Sqr(g1)) * tmp_544 * tmp_545;
   std::complex<double> tmp_548;
   std::complex<double> tmp_550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_550 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_548 += tmp_550;
   std::complex<double> tmp_549;
   std::complex<double> tmp_551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_551 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_549 += tmp_551;
   result += (-0.1*Sqr(g1)) * tmp_548 * tmp_549;
   if (gO1 < 3) {
      std::complex<double> tmp_552;
      std::complex<double> tmp_553;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_553 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_552 += tmp_553;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_552;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_554;
      std::complex<double> tmp_555;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_555 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_554 += tmp_555;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_554;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_556;
      std::complex<double> tmp_557;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_557 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_556 += tmp_557;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_556;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_558;
      std::complex<double> tmp_559;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_559 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_558 += tmp_559;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_558;
   }
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
            tmp_565 += Conj(Yv(j3,j4))*Conj(ZV(gI2,3 + j3));
         }
         tmp_564 += tmp_565;
         tmp_563 += (ZV(gI1,j4)) * tmp_564;
      }
      tmp_561 += tmp_563;
      result += (-1) * tmp_560 * tmp_561;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_566;
      std::complex<double> tmp_568;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_569;
         std::complex<double> tmp_570;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_570 += Yv(j1,j2)*ZV(gI1,3 + j1);
         }
         tmp_569 += tmp_570;
         tmp_568 += (Conj(ZV(gI2,j2))) * tmp_569;
      }
      tmp_566 += tmp_568;
      std::complex<double> tmp_567;
      std::complex<double> tmp_571;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_571 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_567 += tmp_571;
      result += (-1) * tmp_566 * tmp_567;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_572;
   std::complex<double> tmp_573;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_574;
      std::complex<double> tmp_575;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_575 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_574 += tmp_575;
      tmp_573 += (Conj(ZU(gI1,j2))) * tmp_574;
   }
   tmp_572 += tmp_573;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,0))
      * tmp_572;
   if (gO2 < 3) {
      std::complex<double> tmp_576;
      std::complex<double> tmp_577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_577 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_576 += tmp_577;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,0)) *
         tmp_576;
   }
   std::complex<double> tmp_578;
   std::complex<double> tmp_579;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_580;
      std::complex<double> tmp_581;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_581 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_580 += tmp_581;
      tmp_579 += (Conj(ZU(gI1,j2))) * tmp_580;
   }
   tmp_578 += tmp_579;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) * tmp_578
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_582;
      std::complex<double> tmp_583;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_583 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_582 += tmp_583;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_582;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_584;
   std::complex<double> tmp_585;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_585 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_584 += tmp_585;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_584;
   std::complex<double> tmp_586;
   std::complex<double> tmp_587;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_588;
      std::complex<double> tmp_589;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_589 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_588 += tmp_589;
      tmp_587 += (Conj(ZU(gI1,j2))) * tmp_588;
   }
   tmp_586 += tmp_587;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,0)) * tmp_586;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_590;
      std::complex<double> tmp_591;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_591 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_590 += tmp_591;
      result += (0.7071067811865475*Mu*ZH(gI2,0)) * tmp_590;
   }
   std::complex<double> tmp_592;
   std::complex<double> tmp_593;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_593 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_592 += tmp_593;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_592;
   std::complex<double> tmp_594;
   std::complex<double> tmp_595;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_596;
      std::complex<double> tmp_597;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_597 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_596 += tmp_597;
      tmp_595 += (Conj(ZU(gI1,j2))) * tmp_596;
   }
   tmp_594 += tmp_595;
   result += (-0.7071067811865475*ZH(gI2,1)) * tmp_594;
   std::complex<double> tmp_598;
   std::complex<double> tmp_599;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_600;
      std::complex<double> tmp_601;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_602;
         std::complex<double> tmp_603;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_603 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_602 += tmp_603;
         tmp_601 += (KroneckerDelta(gO2,3 + j2)) * tmp_602;
      }
      tmp_600 += tmp_601;
      tmp_599 += (Conj(ZU(gI1,3 + j3))) * tmp_600;
   }
   tmp_598 += tmp_599;
   result += (-(vu*ZH(gI2,1))) * tmp_598;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_604;
      std::complex<double> tmp_605;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_605 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_604 += tmp_605;
      result += (-0.7071067811865475*ZH(gI2,1)) * tmp_604;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_606;
      std::complex<double> tmp_607;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_608;
         std::complex<double> tmp_609;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_609 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_608 += tmp_609;
         tmp_607 += (Conj(ZU(gI1,j2))) * tmp_608;
      }
      tmp_606 += tmp_607;
      result += (-(vu*ZH(gI2,1))) * tmp_606;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_610;
   std::complex<double> tmp_611;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_611 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_610 += tmp_611;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_610;

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

   std::complex<double> tmp_612;
   std::complex<double> tmp_613;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_613 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_612 += tmp_613;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_612;
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

   std::complex<double> tmp_614;
   std::complex<double> tmp_615;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_615 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_614 += tmp_615;
   result += (-0.5163977794943222*g1*Sin(ThetaW())) * tmp_614;
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

   std::complex<double> tmp_616;
   std::complex<double> tmp_617;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_617 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_616 += tmp_617;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_616;
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

   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_619 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_618 += tmp_619;
   result += (0.3*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_618;
   std::complex<double> tmp_620;
   std::complex<double> tmp_621;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_622;
      std::complex<double> tmp_623;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_624;
         std::complex<double> tmp_625;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_625 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_624 += tmp_625;
         tmp_623 += (KroneckerDelta(gO2,3 + j2)) * tmp_624;
      }
      tmp_622 += tmp_623;
      tmp_621 += (KroneckerDelta(gO1,3 + j3)) * tmp_622;
   }
   tmp_620 += tmp_621;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_620;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_626 += tmp_627;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_626;
   }
   std::complex<double> tmp_628;
   std::complex<double> tmp_629;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_629 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_628 += tmp_629;
   result += (-0.3*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_628;
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

   std::complex<double> tmp_630;
   std::complex<double> tmp_631;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_631 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_630 += tmp_631;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_630;
   std::complex<double> tmp_632;
   std::complex<double> tmp_633;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_634;
      std::complex<double> tmp_635;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_636;
         std::complex<double> tmp_637;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_637 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_636 += tmp_637;
         tmp_635 += (KroneckerDelta(gO2,3 + j2)) * tmp_636;
      }
      tmp_634 += tmp_635;
      tmp_633 += (KroneckerDelta(gO1,3 + j3)) * tmp_634;
   }
   tmp_632 += tmp_633;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_632;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_638;
   std::complex<double> tmp_639;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_639 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_638 += tmp_639;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_638;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_640;
      std::complex<double> tmp_641;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_641 += Conj(Yv(j1,gO2))*Yv(j1,gO1);
      }
      tmp_640 += tmp_641;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_640;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_642;
   std::complex<double> tmp_643;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_643 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_642 += tmp_643;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_642;
   std::complex<double> tmp_644;
   std::complex<double> tmp_645;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_646;
      std::complex<double> tmp_647;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_648;
         std::complex<double> tmp_649;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_649 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_648 += tmp_649;
         tmp_647 += (KroneckerDelta(gO2,3 + j2)) * tmp_648;
      }
      tmp_646 += tmp_647;
      tmp_645 += (KroneckerDelta(gO1,3 + j3)) * tmp_646;
   }
   tmp_644 += tmp_645;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_644;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_650;
      std::complex<double> tmp_651;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_651 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_650 += tmp_651;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_650;
   }
   std::complex<double> tmp_652;
   std::complex<double> tmp_653;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_653 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_652 += tmp_653;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_652;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_654;
   std::complex<double> tmp_655;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_655 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_654 += tmp_655;
   result += (-1.0954451150103321*g1*ZN(gI2,0)) * tmp_654;
   if (gO2 < 3) {
      std::complex<double> tmp_656;
      std::complex<double> tmp_657;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_657 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_656 += tmp_657;
      result += (-ZN(gI2,2)) * tmp_656;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_658;
   std::complex<double> tmp_659;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_660;
      std::complex<double> tmp_661;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_661 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_660 += tmp_661;
      tmp_659 += (Conj(ZEL(gI1,j2))) * tmp_660;
   }
   tmp_658 += tmp_659;
   result += (-Conj(ZN(gI2,2))) * tmp_658;
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_662;
   std::complex<double> tmp_664;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_664 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_662 += tmp_664;
   std::complex<double> tmp_663;
   std::complex<double> tmp_665;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_665 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_663 += tmp_665;
   result += (-0.05*Sqr(g1)) * tmp_662 * tmp_663;
   std::complex<double> tmp_666;
   std::complex<double> tmp_668;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_668 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_666 += tmp_668;
   std::complex<double> tmp_667;
   std::complex<double> tmp_669;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_669 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_667 += tmp_669;
   result += (-0.1*Sqr(g1)) * tmp_666 * tmp_667;
   std::complex<double> tmp_670;
   std::complex<double> tmp_672;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_672 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_670 += tmp_672;
   std::complex<double> tmp_671;
   std::complex<double> tmp_673;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_673 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_671 += tmp_673;
   result += (-0.05*Sqr(g1)) * tmp_670 * tmp_671;
   std::complex<double> tmp_674;
   std::complex<double> tmp_676;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_676 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_674 += tmp_676;
   std::complex<double> tmp_675;
   std::complex<double> tmp_677;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_677 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_675 += tmp_677;
   result += (-0.1*Sqr(g1)) * tmp_674 * tmp_675;
   if (gO1 < 3) {
      std::complex<double> tmp_678;
      std::complex<double> tmp_679;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_679 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_678 += tmp_679;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_678;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_680;
      std::complex<double> tmp_681;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_681 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_680 += tmp_681;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_680;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_682;
      std::complex<double> tmp_683;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_683 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_682 += tmp_683;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_682;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_684;
      std::complex<double> tmp_685;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_685 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_684 += tmp_685;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_684;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_686;
      std::complex<double> tmp_687;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_687 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_686 += tmp_687;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_686;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_689 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_688 += tmp_689;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_688;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_690;
      std::complex<double> tmp_692;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_692 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_690 += tmp_692;
      std::complex<double> tmp_691;
      std::complex<double> tmp_693;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_694;
         std::complex<double> tmp_695;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_695 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_694 += tmp_695;
         tmp_693 += (ZD(gI1,j4)) * tmp_694;
      }
      tmp_691 += tmp_693;
      result += (-1) * tmp_690 * tmp_691;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_696;
      std::complex<double> tmp_698;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_699;
         std::complex<double> tmp_700;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_700 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_699 += tmp_700;
         tmp_698 += (Conj(ZD(gI2,j2))) * tmp_699;
      }
      tmp_696 += tmp_698;
      std::complex<double> tmp_697;
      std::complex<double> tmp_701;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_701 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_697 += tmp_701;
      result += (-1) * tmp_696 * tmp_697;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_702;
   std::complex<double> tmp_704;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_704 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_702 += tmp_704;
   std::complex<double> tmp_703;
   std::complex<double> tmp_705;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_705 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_703 += tmp_705;
   result += (-0.3*Sqr(g1)) * tmp_702 * tmp_703;
   std::complex<double> tmp_706;
   std::complex<double> tmp_708;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_708 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_706 += tmp_708;
   std::complex<double> tmp_707;
   std::complex<double> tmp_709;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_709 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_707 += tmp_709;
   result += (0.15*Sqr(g1)) * tmp_706 * tmp_707;
   std::complex<double> tmp_710;
   std::complex<double> tmp_712;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_712 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_710 += tmp_712;
   std::complex<double> tmp_711;
   std::complex<double> tmp_713;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_713 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_711 += tmp_713;
   result += (-0.3*Sqr(g1)) * tmp_710 * tmp_711;
   std::complex<double> tmp_714;
   std::complex<double> tmp_716;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_716 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_714 += tmp_716;
   std::complex<double> tmp_715;
   std::complex<double> tmp_717;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_717 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_715 += tmp_717;
   result += (0.15*Sqr(g1)) * tmp_714 * tmp_715;
   std::complex<double> tmp_718;
   std::complex<double> tmp_720;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_720 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_718 += tmp_720;
   std::complex<double> tmp_719;
   std::complex<double> tmp_721;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_721 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_719 += tmp_721;
   result += (-0.3*Sqr(g1)) * tmp_718 * tmp_719;
   std::complex<double> tmp_722;
   std::complex<double> tmp_724;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_724 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_722 += tmp_724;
   std::complex<double> tmp_723;
   std::complex<double> tmp_725;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_725 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_723 += tmp_725;
   result += (-0.3*Sqr(g1)) * tmp_722 * tmp_723;
   std::complex<double> tmp_726;
   std::complex<double> tmp_728;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_729;
      std::complex<double> tmp_730;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_730 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_729 += tmp_730;
      tmp_728 += (Conj(ZE(gI2,j2))) * tmp_729;
   }
   tmp_726 += tmp_728;
   std::complex<double> tmp_727;
   std::complex<double> tmp_731;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_732;
      std::complex<double> tmp_733;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_733 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_732 += tmp_733;
      tmp_731 += (ZE(gI1,j4)) * tmp_732;
   }
   tmp_727 += tmp_731;
   result += (-1) * tmp_726 * tmp_727;
   if (gO1 < 3) {
      std::complex<double> tmp_734;
      std::complex<double> tmp_735;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_735 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_734 += tmp_735;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_734;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_736;
      std::complex<double> tmp_737;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_737 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_736 += tmp_737;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_736;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_738;
      std::complex<double> tmp_739;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_739 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_738 += tmp_739;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_738;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_740;
      std::complex<double> tmp_741;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_741 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_740 += tmp_741;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_740;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_742;
      std::complex<double> tmp_743;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_743 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_742 += tmp_743;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_742;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_744;
      std::complex<double> tmp_745;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_745 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_744 += tmp_745;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_744;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_746;
      std::complex<double> tmp_748;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_748 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_746 += tmp_748;
      std::complex<double> tmp_747;
      std::complex<double> tmp_749;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_750;
         std::complex<double> tmp_751;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_751 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_750 += tmp_751;
         tmp_749 += (ZE(gI1,j4)) * tmp_750;
      }
      tmp_747 += tmp_749;
      result += (-1) * tmp_746 * tmp_747;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_752;
      std::complex<double> tmp_753;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_753 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_752 += tmp_753;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_752;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_754;
      std::complex<double> tmp_755;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_755 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_754 += tmp_755;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_754;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_756;
      std::complex<double> tmp_758;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_759;
         std::complex<double> tmp_760;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_760 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_759 += tmp_760;
         tmp_758 += (Conj(ZE(gI2,j2))) * tmp_759;
      }
      tmp_756 += tmp_758;
      std::complex<double> tmp_757;
      std::complex<double> tmp_761;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_761 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_757 += tmp_761;
      result += (-1) * tmp_756 * tmp_757;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_762;
      std::complex<double> tmp_764;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_764 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_762 += tmp_764;
      std::complex<double> tmp_763;
      std::complex<double> tmp_765;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_765 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_763 += tmp_765;
      result += (-1) * tmp_762 * tmp_763;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_766;
      std::complex<double> tmp_767;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_767 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_766 += tmp_767;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_766;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_768;
      std::complex<double> tmp_769;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_769 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_768 += tmp_769;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_768;
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

   std::complex<double> tmp_770;
   std::complex<double> tmp_772;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_772 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_770 += tmp_772;
   std::complex<double> tmp_771;
   std::complex<double> tmp_773;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_773 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_771 += tmp_773;
   result += (-0.05*Sqr(g1)) * tmp_770 * tmp_771;
   std::complex<double> tmp_774;
   std::complex<double> tmp_776;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_776 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_774 += tmp_776;
   std::complex<double> tmp_775;
   std::complex<double> tmp_777;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_777 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_775 += tmp_777;
   result += (0.2*Sqr(g1)) * tmp_774 * tmp_775;
   std::complex<double> tmp_778;
   std::complex<double> tmp_780;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_780 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_778 += tmp_780;
   std::complex<double> tmp_779;
   std::complex<double> tmp_781;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_781 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_779 += tmp_781;
   result += (-0.05*Sqr(g1)) * tmp_778 * tmp_779;
   std::complex<double> tmp_782;
   std::complex<double> tmp_784;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_784 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_782 += tmp_784;
   std::complex<double> tmp_783;
   std::complex<double> tmp_785;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_785 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_783 += tmp_785;
   result += (0.2*Sqr(g1)) * tmp_782 * tmp_783;
   if (gO1 < 3) {
      std::complex<double> tmp_786;
      std::complex<double> tmp_787;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_787 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_786 += tmp_787;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_786;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_788;
      std::complex<double> tmp_789;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_789 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_788 += tmp_789;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_788;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_790;
      std::complex<double> tmp_791;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_791 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_790 += tmp_791;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_790;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_792;
      std::complex<double> tmp_793;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_793 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_792 += tmp_793;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_792;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_794;
      std::complex<double> tmp_795;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_795 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_794 += tmp_795;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_794;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_796;
      std::complex<double> tmp_797;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_797 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_796 += tmp_797;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_796;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_798;
   std::complex<double> tmp_800;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_800 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_798 += tmp_800;
   std::complex<double> tmp_799;
   std::complex<double> tmp_801;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_801 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_799 += tmp_801;
   result += (0.15*Sqr(g1)) * tmp_798 * tmp_799;
   std::complex<double> tmp_802;
   std::complex<double> tmp_804;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_804 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_802 += tmp_804;
   std::complex<double> tmp_803;
   std::complex<double> tmp_805;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_805 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_803 += tmp_805;
   result += (0.15*Sqr(g1)) * tmp_802 * tmp_803;
   std::complex<double> tmp_806;
   std::complex<double> tmp_808;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_809;
      std::complex<double> tmp_810;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_810 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_809 += tmp_810;
      tmp_808 += (Conj(ZV(gI2,j2))) * tmp_809;
   }
   tmp_806 += tmp_808;
   std::complex<double> tmp_807;
   std::complex<double> tmp_811;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_812;
      std::complex<double> tmp_813;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_813 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_812 += tmp_813;
      tmp_811 += (ZV(gI1,j4)) * tmp_812;
   }
   tmp_807 += tmp_811;
   result += (-1) * tmp_806 * tmp_807;
   if (gO1 < 3) {
      std::complex<double> tmp_814;
      std::complex<double> tmp_815;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_815 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_814 += tmp_815;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_814;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_816;
      std::complex<double> tmp_817;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_817 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_816 += tmp_817;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_816;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_818;
      std::complex<double> tmp_819;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_819 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_818 += tmp_819;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_818;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_820;
      std::complex<double> tmp_821;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_821 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_820 += tmp_821;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_820;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_824;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_824 += Yv(j1,gO1)*ZV(gI1,3 + j1);
      }
      tmp_822 += tmp_824;
      std::complex<double> tmp_823;
      std::complex<double> tmp_825;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_825 += Conj(Yv(j3,gO2))*Conj(ZV(gI2,3 + j3));
      }
      tmp_823 += tmp_825;
      result += (-1) * tmp_822 * tmp_823;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeconjSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_826;
      std::complex<double> tmp_827;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_828;
         std::complex<double> tmp_829;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_829 += Conj(Yv(j1,gO2))*Mv(j1,j2);
         }
         tmp_828 += tmp_829;
         tmp_827 += (ZV(gI1,3 + j2)) * tmp_828;
      }
      tmp_826 += tmp_827;
      result += (0.5*ZP(gI2,1)) * tmp_826;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_830;
      std::complex<double> tmp_831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_832;
         std::complex<double> tmp_833;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_833 += Conj(Yv(j1,gO2))*Mv(j2,j1);
         }
         tmp_832 += tmp_833;
         tmp_831 += (ZV(gI1,3 + j2)) * tmp_832;
      }
      tmp_830 += tmp_831;
      result += (0.5*ZP(gI2,1)) * tmp_830;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_834;
   std::complex<double> tmp_835;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_836;
      std::complex<double> tmp_837;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_837 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_836 += tmp_837;
      tmp_835 += (Conj(ZE(gI1,j2))) * tmp_836;
   }
   tmp_834 += tmp_835;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_834
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_838;
      std::complex<double> tmp_839;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_839 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_838 += tmp_839;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_838;
   }
   std::complex<double> tmp_840;
   std::complex<double> tmp_841;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_842;
      std::complex<double> tmp_843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_843 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_842 += tmp_843;
      tmp_841 += (Conj(ZE(gI1,j2))) * tmp_842;
   }
   tmp_840 += tmp_841;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1))
      * tmp_840;
   if (gO2 < 3) {
      std::complex<double> tmp_844;
      std::complex<double> tmp_845;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_845 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_844 += tmp_845;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_844;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_846;
   std::complex<double> tmp_847;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_847 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_846 += tmp_847;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_846;
   std::complex<double> tmp_848;
   std::complex<double> tmp_849;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_850;
      std::complex<double> tmp_851;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_851 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_850 += tmp_851;
      tmp_849 += (Conj(ZE(gI1,j2))) * tmp_850;
   }
   tmp_848 += tmp_849;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_848;
   std::complex<double> tmp_852;
   std::complex<double> tmp_853;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_854;
      std::complex<double> tmp_855;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_856;
         std::complex<double> tmp_857;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_857 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_856 += tmp_857;
         tmp_855 += (KroneckerDelta(gO2,3 + j2)) * tmp_856;
      }
      tmp_854 += tmp_855;
      tmp_853 += (Conj(ZE(gI1,3 + j3))) * tmp_854;
   }
   tmp_852 += tmp_853;
   result += (-(vd*ZH(gI2,0))) * tmp_852;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_858;
      std::complex<double> tmp_859;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_859 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_858 += tmp_859;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_858;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_860;
      std::complex<double> tmp_861;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_862;
         std::complex<double> tmp_863;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_863 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_862 += tmp_863;
         tmp_861 += (Conj(ZE(gI1,j2))) * tmp_862;
      }
      tmp_860 += tmp_861;
      result += (-(vd*ZH(gI2,0))) * tmp_860;
   }
   std::complex<double> tmp_864;
   std::complex<double> tmp_865;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_865 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_864 += tmp_865;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_864;
   std::complex<double> tmp_866;
   std::complex<double> tmp_867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_868;
      std::complex<double> tmp_869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_869 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_868 += tmp_869;
      tmp_867 += (Conj(ZE(gI1,j2))) * tmp_868;
   }
   tmp_866 += tmp_867;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_866;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_870;
      std::complex<double> tmp_871;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_871 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_870 += tmp_871;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_870;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_872;
   std::complex<double> tmp_873;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_874;
      std::complex<double> tmp_875;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_875 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_874 += tmp_875;
      tmp_873 += (Conj(ZV(gI1,j2))) * tmp_874;
   }
   tmp_872 += tmp_873;
   result += (ZP(gI2,0)) * tmp_872;
   std::complex<double> tmp_876;
   std::complex<double> tmp_877;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_880;
         std::complex<double> tmp_881;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_881 += Conj(Yv(j3,j1))*Ye(j2,j1);
         }
         tmp_880 += tmp_881;
         tmp_879 += (KroneckerDelta(gO2,3 + j2)) * tmp_880;
      }
      tmp_878 += tmp_879;
      tmp_877 += (Conj(ZV(gI1,3 + j3))) * tmp_878;
   }
   tmp_876 += tmp_877;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_876;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_882;
      std::complex<double> tmp_883;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_883 += Conj(Yv(j1,gO2))*Conj(ZV(gI1,3 + j1));
      }
      tmp_882 += tmp_883;
      result += (Mu*ZP(gI2,0)) * tmp_882;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_884;
      std::complex<double> tmp_885;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_886;
         std::complex<double> tmp_887;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_887 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_886 += tmp_887;
         tmp_885 += (Conj(ZV(gI1,j2))) * tmp_886;
      }
      tmp_884 += tmp_885;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_884;
   }
   std::complex<double> tmp_888;
   std::complex<double> tmp_889;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_890;
      std::complex<double> tmp_891;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_891 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_890 += tmp_891;
      tmp_889 += (Conj(ZV(gI1,j2))) * tmp_890;
   }
   tmp_888 += tmp_889;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_888;
   std::complex<double> tmp_892;
   std::complex<double> tmp_893;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_894;
      std::complex<double> tmp_895;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_896;
         std::complex<double> tmp_897;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_897 += Conj(Yv(j3,j1))*Ye(j2,j1);
         }
         tmp_896 += tmp_897;
         tmp_895 += (KroneckerDelta(gO2,3 + j2)) * tmp_896;
      }
      tmp_894 += tmp_895;
      tmp_893 += (Conj(ZV(gI1,3 + j3))) * tmp_894;
   }
   tmp_892 += tmp_893;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_892;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_898;
      std::complex<double> tmp_899;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_899 += Conj(ZV(gI1,3 + j1))*Conj(TYv(j1,gO2));
      }
      tmp_898 += tmp_899;
      result += (ZP(gI2,1)) * tmp_898;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_900;
      std::complex<double> tmp_901;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_902;
         std::complex<double> tmp_903;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_903 += Conj(Yv(j1,gO2))*Yv(j1,j2);
         }
         tmp_902 += tmp_903;
         tmp_901 += (Conj(ZV(gI1,j2))) * tmp_902;
      }
      tmp_900 += tmp_901;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_900;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_904;
      std::complex<double> tmp_905;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_905 += Conj(Yv(j1,gO2))*UV(gI1,3 + j1);
      }
      tmp_904 += tmp_905;
      result += (UP(gI2,1)) * tmp_904;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_906;
   std::complex<double> tmp_907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_908;
      std::complex<double> tmp_909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_909 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_908 += tmp_909;
      tmp_907 += (Conj(UV(gI1,j2))) * tmp_908;
   }
   tmp_906 += tmp_907;
   result += (Conj(UM(gI2,1))) * tmp_906;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(UV(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVPSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_910;
   std::complex<double> tmp_911;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_911 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_910 += tmp_911;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_910;
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

   std::complex<double> tmp_912;
   std::complex<double> tmp_913;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_913 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_912 += tmp_913;
   result += (0.7745966692414834*g1*Sin(ThetaW())) * tmp_912;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(ThetaW());
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

std::complex<double> CLASSNAME::CpUSvconjUSvVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += 0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()));
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

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   std::complex<double> tmp_914;
   std::complex<double> tmp_915;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_916;
      std::complex<double> tmp_917;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_918;
         std::complex<double> tmp_919;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_919 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_918 += tmp_919;
         tmp_917 += (KroneckerDelta(gO2,3 + j2)) * tmp_918;
      }
      tmp_916 += tmp_917;
      tmp_915 += (KroneckerDelta(gO1,3 + j3)) * tmp_916;
   }
   tmp_914 += tmp_915;
   result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_914;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_920;
      std::complex<double> tmp_921;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_921 += Conj(Yv(j1,gO2))*Yv(j1,gO1);
      }
      tmp_920 += tmp_921;
      result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_920;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_922;
      std::complex<double> tmp_923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_923 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_922 += tmp_923;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_922;
   }
   std::complex<double> tmp_924;
   std::complex<double> tmp_925;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_926;
      std::complex<double> tmp_927;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_928;
         std::complex<double> tmp_929;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_929 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_928 += tmp_929;
         tmp_927 += (KroneckerDelta(gO2,3 + j2)) * tmp_928;
      }
      tmp_926 += tmp_927;
      tmp_925 += (KroneckerDelta(gO1,3 + j3)) * tmp_926;
   }
   tmp_924 += tmp_925;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_924;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   std::complex<double> tmp_930;
   std::complex<double> tmp_931;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_932;
      std::complex<double> tmp_933;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_934;
         std::complex<double> tmp_935;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_935 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_934 += tmp_935;
         tmp_933 += (KroneckerDelta(gO2,3 + j2)) * tmp_934;
      }
      tmp_932 += tmp_933;
      tmp_931 += (KroneckerDelta(gO1,3 + j3)) * tmp_932;
   }
   tmp_930 += tmp_931;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_930;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_936;
      std::complex<double> tmp_937;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_937 += Conj(Yv(j1,gO2))*Yv(j1,gO1);
      }
      tmp_936 += tmp_937;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_936;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_938;
      std::complex<double> tmp_939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_939 += Conj(Ye(j1,gO2))*ZER(gI2,j1);
      }
      tmp_938 += tmp_939;
      result += (UM(gI1,1)) * tmp_938;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_942;
      std::complex<double> tmp_943;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_943 += KroneckerDelta(gO1,3 + j1)*Yv(j1,j2);
      }
      tmp_942 += tmp_943;
      tmp_941 += (Conj(ZEL(gI2,j2))) * tmp_942;
   }
   tmp_940 += tmp_941;
   result += (Conj(UP(gI1,1))) * tmp_940;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZEL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_944;
   std::complex<double> tmp_945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_946;
      std::complex<double> tmp_947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_947 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_946 += tmp_947;
      tmp_945 += (Conj(ZE(gI2,j2))) * tmp_946;
   }
   tmp_944 += tmp_945;
   result += (Conj(Mu)*ZP(gI1,0)) * tmp_944;
   std::complex<double> tmp_948;
   std::complex<double> tmp_949;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_950;
      std::complex<double> tmp_951;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_952;
         std::complex<double> tmp_953;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_953 += Conj(Ye(j3,j1))*Yv(j2,j1);
         }
         tmp_952 += tmp_953;
         tmp_951 += (KroneckerDelta(gO2,3 + j2)) * tmp_952;
      }
      tmp_950 += tmp_951;
      tmp_949 += (Conj(ZE(gI2,3 + j3))) * tmp_950;
   }
   tmp_948 += tmp_949;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_948;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_954;
      std::complex<double> tmp_955;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_955 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_954 += tmp_955;
      result += (ZP(gI1,0)) * tmp_954;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_956;
      std::complex<double> tmp_957;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_958;
         std::complex<double> tmp_959;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_959 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_958 += tmp_959;
         tmp_957 += (Conj(ZE(gI2,j2))) * tmp_958;
      }
      tmp_956 += tmp_957;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_956;
   }
   std::complex<double> tmp_960;
   std::complex<double> tmp_961;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_962;
      std::complex<double> tmp_963;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_963 += KroneckerDelta(gO2,3 + j1)*TYv(j1,j2);
      }
      tmp_962 += tmp_963;
      tmp_961 += (Conj(ZE(gI2,j2))) * tmp_962;
   }
   tmp_960 += tmp_961;
   result += (ZP(gI1,1)) * tmp_960;
   std::complex<double> tmp_964;
   std::complex<double> tmp_965;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_966;
      std::complex<double> tmp_967;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_968;
         std::complex<double> tmp_969;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_969 += Conj(Ye(j3,j1))*Yv(j2,j1);
         }
         tmp_968 += tmp_969;
         tmp_967 += (KroneckerDelta(gO2,3 + j2)) * tmp_968;
      }
      tmp_966 += tmp_967;
      tmp_965 += (Conj(ZE(gI2,3 + j3))) * tmp_966;
   }
   tmp_964 += tmp_965;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_964;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_970;
      std::complex<double> tmp_971;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_971 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_970 += tmp_971;
      result += (Mu*ZP(gI1,1)) * tmp_970;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_972;
      std::complex<double> tmp_973;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_974;
         std::complex<double> tmp_975;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_975 += Conj(Yv(j1,gO2))*Yv(j1,j2);
         }
         tmp_974 += tmp_975;
         tmp_973 += (Conj(ZE(gI2,j2))) * tmp_974;
      }
      tmp_972 += tmp_973;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_972;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_976;
      std::complex<double> tmp_977;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_977 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_976 += tmp_977;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_976;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_978;
      std::complex<double> tmp_979;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_979 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_978 += tmp_979;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_978;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_980;
      std::complex<double> tmp_981;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_981 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_980 += tmp_981;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_980;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_982;
      std::complex<double> tmp_983;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_983 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_982 += tmp_983;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_982;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_984;
      std::complex<double> tmp_985;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_985 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_984 += tmp_985;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_984;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_986;
      std::complex<double> tmp_987;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_987 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_986 += tmp_987;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_986;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_988;
   std::complex<double> tmp_990;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_991;
      std::complex<double> tmp_992;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_992 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_991 += tmp_992;
      tmp_990 += (Conj(ZE(gI2,j2))) * tmp_991;
   }
   tmp_988 += tmp_990;
   std::complex<double> tmp_989;
   std::complex<double> tmp_993;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_994;
      std::complex<double> tmp_995;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_995 += Conj(Yv(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_994 += tmp_995;
      tmp_993 += (ZE(gI1,j4)) * tmp_994;
   }
   tmp_989 += tmp_993;
   result += (-1) * tmp_988 * tmp_989;
   if (gO1 < 3) {
      std::complex<double> tmp_996;
      std::complex<double> tmp_997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_997 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_996 += tmp_997;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_996;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_998;
      std::complex<double> tmp_999;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_999 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_998 += tmp_999;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_998;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1000;
      std::complex<double> tmp_1001;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1001 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1000 += tmp_1001;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1000;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1002;
      std::complex<double> tmp_1003;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1003 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1002 += tmp_1003;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1002;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1004;
      std::complex<double> tmp_1005;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1005 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1004 += tmp_1005;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1004;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1006;
      std::complex<double> tmp_1007;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1007 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1006 += tmp_1007;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1006;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1008;
      std::complex<double> tmp_1010;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1010 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_1008 += tmp_1010;
      std::complex<double> tmp_1009;
      std::complex<double> tmp_1011;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1011 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_1009 += tmp_1011;
      result += (-1) * tmp_1008 * tmp_1009;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1012;
      std::complex<double> tmp_1013;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1013 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1012 += tmp_1013;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1012;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1014;
      std::complex<double> tmp_1015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1015 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1014 += tmp_1015;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1014;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1016;
      std::complex<double> tmp_1017;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1017 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1016 += tmp_1017;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1016;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1018;
      std::complex<double> tmp_1019;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1019 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1018 += tmp_1019;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1018;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1020;
      std::complex<double> tmp_1021;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1021 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1020 += tmp_1021;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1020;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1022;
      std::complex<double> tmp_1023;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1023 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1022 += tmp_1023;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1022;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1024;
      std::complex<double> tmp_1026;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1026 += KroneckerDelta(gO2,3 + j1)*Yv(j1,gO1);
      }
      tmp_1024 += tmp_1026;
      std::complex<double> tmp_1025;
      std::complex<double> tmp_1027;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1028;
         std::complex<double> tmp_1029;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1029 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_1028 += tmp_1029;
         tmp_1027 += (ZU(gI1,j4)) * tmp_1028;
      }
      tmp_1025 += tmp_1027;
      result += (-1) * tmp_1024 * tmp_1025;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1030;
      std::complex<double> tmp_1032;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1033;
         std::complex<double> tmp_1034;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1034 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_1033 += tmp_1034;
         tmp_1032 += (Conj(ZU(gI2,j2))) * tmp_1033;
      }
      tmp_1030 += tmp_1032;
      std::complex<double> tmp_1031;
      std::complex<double> tmp_1035;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1035 += Conj(Yv(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1031 += tmp_1035;
      result += (-1) * tmp_1030 * tmp_1031;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1036;
   std::complex<double> tmp_1038;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1039;
      std::complex<double> tmp_1040;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1040 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_1039 += tmp_1040;
      tmp_1038 += (Conj(ZV(gI2,j2))) * tmp_1039;
   }
   tmp_1036 += tmp_1038;
   std::complex<double> tmp_1037;
   std::complex<double> tmp_1041;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1042;
      std::complex<double> tmp_1043;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1043 += Conj(Yv(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1042 += tmp_1043;
      tmp_1041 += (ZV(gI1,j4)) * tmp_1042;
   }
   tmp_1037 += tmp_1041;
   result += (-1) * tmp_1036 * tmp_1037;
   if (gO1 < 3) {
      std::complex<double> tmp_1044;
      std::complex<double> tmp_1045;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1045 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1044 += tmp_1045;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1044;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1046;
      std::complex<double> tmp_1047;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1047 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1046 += tmp_1047;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1046;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1048;
      std::complex<double> tmp_1049;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1049 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1048 += tmp_1049;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1048;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1050;
      std::complex<double> tmp_1051;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1051 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1050 += tmp_1051;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1050;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1052;
      std::complex<double> tmp_1054;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1054 += KroneckerDelta(gO2,3 + j1)*Yv(j1,gO1);
      }
      tmp_1052 += tmp_1054;
      std::complex<double> tmp_1053;
      std::complex<double> tmp_1055;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1056;
         std::complex<double> tmp_1057;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1057 += Conj(Yv(j3,j4))*Conj(ZV(gI2,3 + j3));
         }
         tmp_1056 += tmp_1057;
         tmp_1055 += (ZV(gI1,j4)) * tmp_1056;
      }
      tmp_1053 += tmp_1055;
      result += (-1) * tmp_1052 * tmp_1053;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1058;
      std::complex<double> tmp_1060;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1061;
         std::complex<double> tmp_1062;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1062 += Yv(j1,j2)*ZV(gI1,3 + j1);
         }
         tmp_1061 += tmp_1062;
         tmp_1060 += (Conj(ZV(gI2,j2))) * tmp_1061;
      }
      tmp_1058 += tmp_1060;
      std::complex<double> tmp_1059;
      std::complex<double> tmp_1063;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1063 += Conj(Yv(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1059 += tmp_1063;
      result += (-1) * tmp_1058 * tmp_1059;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1064;
      std::complex<double> tmp_1066;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1066 += Yv(j1,gO1)*ZV(gI1,3 + j1);
      }
      tmp_1064 += tmp_1066;
      std::complex<double> tmp_1065;
      std::complex<double> tmp_1067;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1067 += Conj(Yv(j3,gO2))*Conj(ZV(gI2,3 + j3));
      }
      tmp_1065 += tmp_1067;
      result += (-1) * tmp_1064 * tmp_1065;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.15*Conj(ZV(gI2,gO2))*Sqr(g1)*ZV(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjSeHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1068;
   std::complex<double> tmp_1069;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1070;
      std::complex<double> tmp_1071;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1072;
         std::complex<double> tmp_1073;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1073 += Conj(Yv(j1,j3))*Mv(j1,j2);
         }
         tmp_1072 += tmp_1073;
         tmp_1071 += (KroneckerDelta(gO2,3 + j2)) * tmp_1072;
      }
      tmp_1070 += tmp_1071;
      tmp_1069 += (ZE(gI1,j3)) * tmp_1070;
   }
   tmp_1068 += tmp_1069;
   std::complex<double> tmp_1074;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1075;
      std::complex<double> tmp_1076;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1077;
         std::complex<double> tmp_1078;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1078 += Conj(Yv(j1,j3))*Mv(j2,j1);
         }
         tmp_1077 += tmp_1078;
         tmp_1076 += (KroneckerDelta(gO2,3 + j2)) * tmp_1077;
      }
      tmp_1075 += tmp_1076;
      tmp_1074 += (ZE(gI1,j3)) * tmp_1075;
   }
   tmp_1068 += tmp_1074;
   result += (0.5*ZP(gI2,1)) * tmp_1068;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjSvAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1079;
   std::complex<double> tmp_1080;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1081;
      std::complex<double> tmp_1082;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1083;
         std::complex<double> tmp_1084;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1084 += Conj(Yv(j1,j3))*Mv(j1,j2);
         }
         tmp_1083 += tmp_1084;
         tmp_1082 += (KroneckerDelta(gO2,3 + j2)) * tmp_1083;
      }
      tmp_1081 += tmp_1082;
      tmp_1080 += (ZV(gI1,j3)) * tmp_1081;
   }
   tmp_1079 += tmp_1080;
   result += (std::complex<double>(0.,0.35355339059327373)*ZA(gI2,1)) *
      tmp_1079;
   std::complex<double> tmp_1085;
   std::complex<double> tmp_1086;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1087;
      std::complex<double> tmp_1088;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1089;
         std::complex<double> tmp_1090;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1090 += Conj(Yv(j1,j3))*Mv(j2,j1);
         }
         tmp_1089 += tmp_1090;
         tmp_1088 += (KroneckerDelta(gO2,3 + j2)) * tmp_1089;
      }
      tmp_1087 += tmp_1088;
      tmp_1086 += (ZV(gI1,j3)) * tmp_1087;
   }
   tmp_1085 += tmp_1086;
   result += (std::complex<double>(0.,0.35355339059327373)*ZA(gI2,1)) *
      tmp_1085;
   if (gO2 < 3) {
      std::complex<double> tmp_1091;
      std::complex<double> tmp_1092;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1093;
         std::complex<double> tmp_1094;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1094 += Conj(Yv(j1,gO2))*Mv(j1,j2);
         }
         tmp_1093 += tmp_1094;
         tmp_1092 += (ZV(gI1,3 + j2)) * tmp_1093;
      }
      tmp_1091 += tmp_1092;
      result += (std::complex<double>(0.,0.35355339059327373)*ZA(gI2,1)) *
         tmp_1091;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1095;
      std::complex<double> tmp_1096;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1097;
         std::complex<double> tmp_1098;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1098 += Conj(Yv(j1,gO2))*Mv(j2,j1);
         }
         tmp_1097 += tmp_1098;
         tmp_1096 += (ZV(gI1,3 + j2)) * tmp_1097;
      }
      tmp_1095 += tmp_1096;
      result += (std::complex<double>(0.,0.35355339059327373)*ZA(gI2,1)) *
         tmp_1095;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1099;
   std::complex<double> tmp_1100;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1101;
      std::complex<double> tmp_1102;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1103;
         std::complex<double> tmp_1104;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1104 += Conj(Yv(j1,j3))*Mv(j1,j2);
         }
         tmp_1103 += tmp_1104;
         tmp_1102 += (KroneckerDelta(gO2,3 + j2)) * tmp_1103;
      }
      tmp_1101 += tmp_1102;
      tmp_1100 += (ZV(gI1,j3)) * tmp_1101;
   }
   tmp_1099 += tmp_1100;
   result += (-0.35355339059327373*ZH(gI2,1)) * tmp_1099;
   std::complex<double> tmp_1105;
   std::complex<double> tmp_1106;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1107;
      std::complex<double> tmp_1108;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1109;
         std::complex<double> tmp_1110;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1110 += Conj(Yv(j1,j3))*Mv(j2,j1);
         }
         tmp_1109 += tmp_1110;
         tmp_1108 += (KroneckerDelta(gO2,3 + j2)) * tmp_1109;
      }
      tmp_1107 += tmp_1108;
      tmp_1106 += (ZV(gI1,j3)) * tmp_1107;
   }
   tmp_1105 += tmp_1106;
   result += (-0.35355339059327373*ZH(gI2,1)) * tmp_1105;
   if (gO2 < 3) {
      std::complex<double> tmp_1111;
      std::complex<double> tmp_1112;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1113;
         std::complex<double> tmp_1114;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1114 += Conj(Yv(j1,gO2))*Mv(j1,j2);
         }
         tmp_1113 += tmp_1114;
         tmp_1112 += (ZV(gI1,3 + j2)) * tmp_1113;
      }
      tmp_1111 += tmp_1112;
      result += (-0.35355339059327373*ZH(gI2,1)) * tmp_1111;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1115;
      std::complex<double> tmp_1116;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1117;
         std::complex<double> tmp_1118;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1118 += Conj(Yv(j1,gO2))*Mv(j2,j1);
         }
         tmp_1117 += tmp_1118;
         tmp_1116 += (ZV(gI1,3 + j2)) * tmp_1117;
      }
      tmp_1115 += tmp_1116;
      result += (-0.35355339059327373*ZH(gI2,1)) * tmp_1115;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1119;
   std::complex<double> tmp_1120;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1121;
      std::complex<double> tmp_1122;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1122 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_1121 += tmp_1122;
      tmp_1120 += (Conj(ZV(gI1,j2))) * tmp_1121;
   }
   tmp_1119 += tmp_1120;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(Mu)*ZA(gI2,0))
      * tmp_1119;
   if (gO2 < 3) {
      std::complex<double> tmp_1123;
      std::complex<double> tmp_1124;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1124 += Conj(Yv(j1,gO2))*Conj(ZV(gI1,3 + j1));
      }
      tmp_1123 += tmp_1124;
      result += (std::complex<double>(0.,0.7071067811865475)*Mu*ZA(gI2,0)) *
         tmp_1123;
   }
   std::complex<double> tmp_1125;
   std::complex<double> tmp_1126;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1127;
      std::complex<double> tmp_1128;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1128 += KroneckerDelta(gO2,3 + j1)*TYv(j1,j2);
      }
      tmp_1127 += tmp_1128;
      tmp_1126 += (Conj(ZV(gI1,j2))) * tmp_1127;
   }
   tmp_1125 += tmp_1126;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
      tmp_1125;
   if (gO2 < 3) {
      std::complex<double> tmp_1129;
      std::complex<double> tmp_1130;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1130 += Conj(ZV(gI1,3 + j1))*Conj(TYv(j1,gO2));
      }
      tmp_1129 += tmp_1130;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_1129;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1131;
   std::complex<double> tmp_1132;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1133;
      std::complex<double> tmp_1134;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1134 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_1133 += tmp_1134;
      tmp_1132 += (Conj(ZV(gI1,j2))) * tmp_1133;
   }
   tmp_1131 += tmp_1132;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,0)) * tmp_1131;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1135;
      std::complex<double> tmp_1136;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1136 += Conj(Yv(j1,gO2))*Conj(ZV(gI1,3 + j1));
      }
      tmp_1135 += tmp_1136;
      result += (0.7071067811865475*Mu*ZH(gI2,0)) * tmp_1135;
   }
   std::complex<double> tmp_1137;
   std::complex<double> tmp_1138;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1139;
      std::complex<double> tmp_1140;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1140 += KroneckerDelta(gO2,3 + j1)*TYv(j1,j2);
      }
      tmp_1139 += tmp_1140;
      tmp_1138 += (Conj(ZV(gI1,j2))) * tmp_1139;
   }
   tmp_1137 += tmp_1138;
   result += (-0.7071067811865475*ZH(gI2,1)) * tmp_1137;
   std::complex<double> tmp_1141;
   std::complex<double> tmp_1142;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1143;
      std::complex<double> tmp_1144;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1145;
         std::complex<double> tmp_1146;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1146 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_1145 += tmp_1146;
         tmp_1144 += (KroneckerDelta(gO2,3 + j2)) * tmp_1145;
      }
      tmp_1143 += tmp_1144;
      tmp_1142 += (Conj(ZV(gI1,3 + j3))) * tmp_1143;
   }
   tmp_1141 += tmp_1142;
   result += (-(vu*ZH(gI2,1))) * tmp_1141;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1147;
      std::complex<double> tmp_1148;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1148 += Conj(ZV(gI1,3 + j1))*Conj(TYv(j1,gO2));
      }
      tmp_1147 += tmp_1148;
      result += (-0.7071067811865475*ZH(gI2,1)) * tmp_1147;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1149;
      std::complex<double> tmp_1150;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1151;
         std::complex<double> tmp_1152;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1152 += Conj(Yv(j1,gO2))*Yv(j1,j2);
         }
         tmp_1151 += tmp_1152;
         tmp_1150 += (Conj(ZV(gI1,j2))) * tmp_1151;
      }
      tmp_1149 += tmp_1150;
      result += (-(vu*ZH(gI2,1))) * tmp_1149;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1153;
      std::complex<double> tmp_1154;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1154 += Conj(Yv(j1,gO2))*UV(gI1,3 + j1);
      }
      tmp_1153 += tmp_1154;
      result += (-ZN(gI2,3)) * tmp_1153;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1155;
   std::complex<double> tmp_1156;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1157;
      std::complex<double> tmp_1158;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1158 += KroneckerDelta(gO1,3 + j1)*Yv(j1,j2);
      }
      tmp_1157 += tmp_1158;
      tmp_1156 += (Conj(UV(gI1,j2))) * tmp_1157;
   }
   tmp_1155 += tmp_1156;
   result += (-Conj(ZN(gI2,3))) * tmp_1155;
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(UV(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(UV(gI1,gO1))*Conj(ZN(gI2,1));
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

std::complex<double> CLASSNAME::CpconjUSvVZSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW());
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
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,0)*ZA(gI2,0)
      - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(
      gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1)
      + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2
      ,1)))) + KroneckerDelta(1,gO1)*(-5*KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*
      ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr
      (g2))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(1,gO1)*(KroneckerDelta
      (0,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*
      (ZH(gI1,0)*ZH(gI2,0) - 3*ZH(gI1,1)*ZH(gI2,1))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(1,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) +
      KroneckerDelta(0,gO2)*(-3*ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*(3*Sqr(
      g1) + 5*Sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(ZP(gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2)
      )*ZP(gI2,0) + 5*vu*Sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*Sqr(g2)*ZP(gI2,0) +
      vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1)))) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*
      (vu*(3*Sqr(g1) - 5*Sqr(g2))*ZP(gI2,0) - 5*vd*Sqr(g2)*ZP(gI2,1)) - ZP(gI1,1)*
      (5*vd*Sqr(g2)*ZP(gI2,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(1,gO2)*(ZH(gI1,0)*(vu*
      ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(vd*ZH(gI2,0) - 3*vu*ZH(gI2,1))) +
      KroneckerDelta(0,gO2)*(ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,0)*(
      -3*vd*ZH(gI2,0) + vu*ZH(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*(KroneckerDelta(0,gO2)*UM(gI1,1)*UP(gI2,0) +
      KroneckerDelta(1,gO2)*UM(gI1,0)*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*(Conj(UM(gI2,1))*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) + Conj(UM(gI2,0))*Conj(UP(gI1,1))*KroneckerDelta(1,gO1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1159;
   std::complex<double> tmp_1160;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1161;
      std::complex<double> tmp_1162;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1162 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1161 += tmp_1162;
      tmp_1160 += (ZDL(gI1,j2)) * tmp_1161;
   }
   tmp_1159 += tmp_1160;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1159;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1163;
   std::complex<double> tmp_1164;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1165;
      std::complex<double> tmp_1166;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1166 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1165 += tmp_1166;
      tmp_1164 += (Conj(ZDL(gI2,j2))) * tmp_1165;
   }
   tmp_1163 += tmp_1164;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1163;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1167;
   std::complex<double> tmp_1168;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1169;
      std::complex<double> tmp_1170;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1170 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1169 += tmp_1170;
      tmp_1168 += (ZEL(gI1,j2)) * tmp_1169;
   }
   tmp_1167 += tmp_1168;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1167;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1171;
   std::complex<double> tmp_1172;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1173;
      std::complex<double> tmp_1174;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1174 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1173 += tmp_1174;
      tmp_1172 += (Conj(ZEL(gI2,j2))) * tmp_1173;
   }
   tmp_1171 += tmp_1172;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1171;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1175;
   std::complex<double> tmp_1176;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1177;
      std::complex<double> tmp_1178;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1178 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1177 += tmp_1178;
      tmp_1176 += (ZUL(gI1,j2)) * tmp_1177;
   }
   tmp_1175 += tmp_1176;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1175;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1179;
   std::complex<double> tmp_1180;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1181;
      std::complex<double> tmp_1182;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1182 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1181 += tmp_1182;
      tmp_1180 += (Conj(ZUL(gI2,j2))) * tmp_1181;
   }
   tmp_1179 += tmp_1180;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1179;

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO2)*(ZN(gI1,2)*(3.872983346207417*g1*ZN(gI2,
      0) - 5*g2*ZN(gI2,1)) + (3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(
      gI2,2)) - KroneckerDelta(1,gO2)*(ZN(gI1,3)*(3.872983346207417*g1*ZN(gI2,0) -
      5*g2*ZN(gI2,1)) + (3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,
      3)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(Conj(ZN(gI1,2))*(3.872983346207417*g1*Conj(ZN(gI2,0)) - 5*g2*
      Conj(ZN(gI2,1)))*KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2)
      )*KroneckerDelta(0,gO1) - 3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0
      ))*KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*
      KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(
      1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(Conj(ZN(gI2,2))*
      KroneckerDelta(0,gO1) - Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1183;
   std::complex<double> tmp_1184;
   std::complex<double> tmp_1185;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1185 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1184 += tmp_1185;
   tmp_1183 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1184;
   std::complex<double> tmp_1186;
   std::complex<double> tmp_1187;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1187 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1186 += tmp_1187;
   tmp_1183 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1186;
   std::complex<double> tmp_1188;
   std::complex<double> tmp_1189;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1189 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1188 += tmp_1189;
   tmp_1183 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1188;
   std::complex<double> tmp_1190;
   std::complex<double> tmp_1191;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1191 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1190 += tmp_1191;
   tmp_1183 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1190;
   std::complex<double> tmp_1192;
   std::complex<double> tmp_1193;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1193 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1192 += tmp_1193;
   tmp_1183 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1192;
   std::complex<double> tmp_1194;
   std::complex<double> tmp_1195;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1195 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1194 += tmp_1195;
   tmp_1183 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1194;
   std::complex<double> tmp_1196;
   std::complex<double> tmp_1197;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1198;
      std::complex<double> tmp_1199;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1200;
         std::complex<double> tmp_1201;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1201 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1200 += tmp_1201;
         tmp_1199 += (ZD(gI1,3 + j2)) * tmp_1200;
      }
      tmp_1198 += tmp_1199;
      tmp_1197 += (Conj(ZD(gI2,3 + j3))) * tmp_1198;
   }
   tmp_1196 += tmp_1197;
   tmp_1183 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1196;
   std::complex<double> tmp_1202;
   std::complex<double> tmp_1203;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1204;
      std::complex<double> tmp_1205;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1206;
         std::complex<double> tmp_1207;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1207 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1206 += tmp_1207;
         tmp_1205 += (Conj(ZD(gI2,j2))) * tmp_1206;
      }
      tmp_1204 += tmp_1205;
      tmp_1203 += (ZD(gI1,j3)) * tmp_1204;
   }
   tmp_1202 += tmp_1203;
   tmp_1183 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1202;
   result += (std::complex<double>(0,-1)) * tmp_1183;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1208;
   std::complex<double> tmp_1209;
   std::complex<double> tmp_1210;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1210 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1209 += tmp_1210;
   tmp_1208 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1209;
   std::complex<double> tmp_1211;
   std::complex<double> tmp_1212;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1212 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1211 += tmp_1212;
   tmp_1208 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1211;
   std::complex<double> tmp_1213;
   std::complex<double> tmp_1214;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1214 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1213 += tmp_1214;
   tmp_1208 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1213;
   std::complex<double> tmp_1215;
   std::complex<double> tmp_1216;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1216 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1215 += tmp_1216;
   tmp_1208 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1215;
   std::complex<double> tmp_1217;
   std::complex<double> tmp_1218;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1218 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1217 += tmp_1218;
   tmp_1208 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1217;
   std::complex<double> tmp_1219;
   std::complex<double> tmp_1220;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1220 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1219 += tmp_1220;
   tmp_1208 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1219;
   std::complex<double> tmp_1221;
   std::complex<double> tmp_1222;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1223;
      std::complex<double> tmp_1224;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1225;
         std::complex<double> tmp_1226;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1226 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1225 += tmp_1226;
         tmp_1224 += (ZE(gI1,3 + j2)) * tmp_1225;
      }
      tmp_1223 += tmp_1224;
      tmp_1222 += (Conj(ZE(gI2,3 + j3))) * tmp_1223;
   }
   tmp_1221 += tmp_1222;
   tmp_1208 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1221;
   std::complex<double> tmp_1227;
   std::complex<double> tmp_1228;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1229;
      std::complex<double> tmp_1230;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1231;
         std::complex<double> tmp_1232;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1232 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1231 += tmp_1232;
         tmp_1230 += (Conj(ZE(gI2,j2))) * tmp_1231;
      }
      tmp_1229 += tmp_1230;
      tmp_1228 += (ZE(gI1,j3)) * tmp_1229;
   }
   tmp_1227 += tmp_1228;
   tmp_1208 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1227;
   result += (std::complex<double>(0,-1)) * tmp_1208;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1233;
   std::complex<double> tmp_1234;
   std::complex<double> tmp_1235;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1235 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1234 += tmp_1235;
   tmp_1233 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1234;
   std::complex<double> tmp_1236;
   std::complex<double> tmp_1237;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1237 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1236 += tmp_1237;
   tmp_1233 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1236;
   std::complex<double> tmp_1238;
   std::complex<double> tmp_1239;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1239 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1238 += tmp_1239;
   tmp_1233 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1238;
   std::complex<double> tmp_1240;
   std::complex<double> tmp_1241;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1241 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1240 += tmp_1241;
   tmp_1233 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1240;
   std::complex<double> tmp_1242;
   std::complex<double> tmp_1243;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1243 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1242 += tmp_1243;
   tmp_1233 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1242;
   std::complex<double> tmp_1244;
   std::complex<double> tmp_1245;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1245 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1244 += tmp_1245;
   tmp_1233 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1244;
   std::complex<double> tmp_1246;
   std::complex<double> tmp_1247;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1248;
      std::complex<double> tmp_1249;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1250;
         std::complex<double> tmp_1251;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1251 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1250 += tmp_1251;
         tmp_1249 += (ZU(gI1,3 + j2)) * tmp_1250;
      }
      tmp_1248 += tmp_1249;
      tmp_1247 += (Conj(ZU(gI2,3 + j3))) * tmp_1248;
   }
   tmp_1246 += tmp_1247;
   tmp_1233 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1246;
   std::complex<double> tmp_1252;
   std::complex<double> tmp_1253;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1254;
      std::complex<double> tmp_1255;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1256;
         std::complex<double> tmp_1257;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1257 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1256 += tmp_1257;
         tmp_1255 += (Conj(ZU(gI2,j2))) * tmp_1256;
      }
      tmp_1254 += tmp_1255;
      tmp_1253 += (ZU(gI1,j3)) * tmp_1254;
   }
   tmp_1252 += tmp_1253;
   tmp_1233 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1252;
   result += (std::complex<double>(0,-1)) * tmp_1233;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1258;
   std::complex<double> tmp_1259;
   std::complex<double> tmp_1260;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1260 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1259 += tmp_1260;
   tmp_1258 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1259;
   std::complex<double> tmp_1261;
   std::complex<double> tmp_1262;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1262 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1261 += tmp_1262;
   tmp_1258 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1261;
   std::complex<double> tmp_1263;
   std::complex<double> tmp_1264;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1264 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1263 += tmp_1264;
   tmp_1258 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1263;
   std::complex<double> tmp_1265;
   std::complex<double> tmp_1266;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1266 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1265 += tmp_1266;
   tmp_1258 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1265;
   std::complex<double> tmp_1267;
   std::complex<double> tmp_1268;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1269;
      std::complex<double> tmp_1270;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1271;
         std::complex<double> tmp_1272;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1272 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_1271 += tmp_1272;
         tmp_1270 += (ZV(gI1,3 + j2)) * tmp_1271;
      }
      tmp_1269 += tmp_1270;
      tmp_1268 += (Conj(ZV(gI2,3 + j3))) * tmp_1269;
   }
   tmp_1267 += tmp_1268;
   tmp_1258 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1267;
   std::complex<double> tmp_1273;
   std::complex<double> tmp_1274;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1275;
      std::complex<double> tmp_1276;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1277;
         std::complex<double> tmp_1278;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1278 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_1277 += tmp_1278;
         tmp_1276 += (Conj(ZV(gI2,j2))) * tmp_1277;
      }
      tmp_1275 += tmp_1276;
      tmp_1274 += (ZV(gI1,j3)) * tmp_1275;
   }
   tmp_1273 += tmp_1274;
   tmp_1258 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1273;
   result += (std::complex<double>(0,-1)) * tmp_1258;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1279;
   std::complex<double> tmp_1280;
   std::complex<double> tmp_1281;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1281 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1280 += tmp_1281;
   tmp_1279 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1280;
   std::complex<double> tmp_1282;
   std::complex<double> tmp_1283;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1283 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1282 += tmp_1283;
   tmp_1279 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1282;
   std::complex<double> tmp_1284;
   std::complex<double> tmp_1285;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1285 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1284 += tmp_1285;
   tmp_1279 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1284;
   std::complex<double> tmp_1286;
   std::complex<double> tmp_1287;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1287 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1286 += tmp_1287;
   tmp_1279 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1286;
   std::complex<double> tmp_1288;
   std::complex<double> tmp_1289;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1289 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1288 += tmp_1289;
   tmp_1279 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1288;
   std::complex<double> tmp_1290;
   std::complex<double> tmp_1291;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1291 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1290 += tmp_1291;
   tmp_1279 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1290;
   std::complex<double> tmp_1292;
   std::complex<double> tmp_1293;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1294;
      std::complex<double> tmp_1295;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1295 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1294 += tmp_1295;
      tmp_1293 += (Conj(ZD(gI2,j2))) * tmp_1294;
   }
   tmp_1292 += tmp_1293;
   tmp_1279 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1292;
   std::complex<double> tmp_1296;
   std::complex<double> tmp_1297;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1298;
      std::complex<double> tmp_1299;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1299 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1298 += tmp_1299;
      tmp_1297 += (Conj(ZD(gI2,j2))) * tmp_1298;
   }
   tmp_1296 += tmp_1297;
   tmp_1279 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1296;
   std::complex<double> tmp_1300;
   std::complex<double> tmp_1301;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1302;
      std::complex<double> tmp_1303;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1303 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1302 += tmp_1303;
      tmp_1301 += (ZD(gI1,j2)) * tmp_1302;
   }
   tmp_1300 += tmp_1301;
   tmp_1279 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,
      gO2)*Mu) * tmp_1300;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1306;
      std::complex<double> tmp_1307;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1307 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1306 += tmp_1307;
      tmp_1305 += (ZD(gI1,j2)) * tmp_1306;
   }
   tmp_1304 += tmp_1305;
   tmp_1279 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1304;
   std::complex<double> tmp_1308;
   std::complex<double> tmp_1309;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1310;
      std::complex<double> tmp_1311;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1312;
         std::complex<double> tmp_1313;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1313 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1312 += tmp_1313;
         tmp_1311 += (ZD(gI1,3 + j2)) * tmp_1312;
      }
      tmp_1310 += tmp_1311;
      tmp_1309 += (Conj(ZD(gI2,3 + j3))) * tmp_1310;
   }
   tmp_1308 += tmp_1309;
   tmp_1279 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1308
      ;
   std::complex<double> tmp_1314;
   std::complex<double> tmp_1315;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1316;
      std::complex<double> tmp_1317;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1318;
         std::complex<double> tmp_1319;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1319 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1318 += tmp_1319;
         tmp_1317 += (Conj(ZD(gI2,j2))) * tmp_1318;
      }
      tmp_1316 += tmp_1317;
      tmp_1315 += (ZD(gI1,j3)) * tmp_1316;
   }
   tmp_1314 += tmp_1315;
   tmp_1279 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1314
      ;
   result += (std::complex<double>(0,-1)) * tmp_1279;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1320;
   std::complex<double> tmp_1321;
   std::complex<double> tmp_1322;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1322 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1321 += tmp_1322;
   tmp_1320 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1321;
   std::complex<double> tmp_1323;
   std::complex<double> tmp_1324;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1324 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1323 += tmp_1324;
   tmp_1320 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1323;
   std::complex<double> tmp_1325;
   std::complex<double> tmp_1326;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1326 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1325 += tmp_1326;
   tmp_1320 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1325;
   std::complex<double> tmp_1327;
   std::complex<double> tmp_1328;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1328 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1327 += tmp_1328;
   tmp_1320 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1327;
   std::complex<double> tmp_1329;
   std::complex<double> tmp_1330;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1330 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1329 += tmp_1330;
   tmp_1320 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1329;
   std::complex<double> tmp_1331;
   std::complex<double> tmp_1332;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1332 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1331 += tmp_1332;
   tmp_1320 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1331;
   std::complex<double> tmp_1333;
   std::complex<double> tmp_1334;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1335;
      std::complex<double> tmp_1336;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1336 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1335 += tmp_1336;
      tmp_1334 += (Conj(ZE(gI2,j2))) * tmp_1335;
   }
   tmp_1333 += tmp_1334;
   tmp_1320 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1333;
   std::complex<double> tmp_1337;
   std::complex<double> tmp_1338;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1339;
      std::complex<double> tmp_1340;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1340 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1339 += tmp_1340;
      tmp_1338 += (Conj(ZE(gI2,j2))) * tmp_1339;
   }
   tmp_1337 += tmp_1338;
   tmp_1320 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1337;
   std::complex<double> tmp_1341;
   std::complex<double> tmp_1342;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1343;
      std::complex<double> tmp_1344;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1344 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1343 += tmp_1344;
      tmp_1342 += (ZE(gI1,j2)) * tmp_1343;
   }
   tmp_1341 += tmp_1342;
   tmp_1320 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,
      gO2)*Mu) * tmp_1341;
   std::complex<double> tmp_1345;
   std::complex<double> tmp_1346;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1347;
      std::complex<double> tmp_1348;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1348 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1347 += tmp_1348;
      tmp_1346 += (ZE(gI1,j2)) * tmp_1347;
   }
   tmp_1345 += tmp_1346;
   tmp_1320 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1345;
   std::complex<double> tmp_1349;
   std::complex<double> tmp_1350;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1351;
      std::complex<double> tmp_1352;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1353;
         std::complex<double> tmp_1354;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1354 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1353 += tmp_1354;
         tmp_1352 += (ZE(gI1,3 + j2)) * tmp_1353;
      }
      tmp_1351 += tmp_1352;
      tmp_1350 += (Conj(ZE(gI2,3 + j3))) * tmp_1351;
   }
   tmp_1349 += tmp_1350;
   tmp_1320 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1349
      ;
   std::complex<double> tmp_1355;
   std::complex<double> tmp_1356;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1357;
      std::complex<double> tmp_1358;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1359;
         std::complex<double> tmp_1360;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1360 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1359 += tmp_1360;
         tmp_1358 += (Conj(ZE(gI2,j2))) * tmp_1359;
      }
      tmp_1357 += tmp_1358;
      tmp_1356 += (ZE(gI1,j3)) * tmp_1357;
   }
   tmp_1355 += tmp_1356;
   tmp_1320 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1355
      ;
   result += (std::complex<double>(0,-1)) * tmp_1320;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1361;
   std::complex<double> tmp_1362;
   std::complex<double> tmp_1363;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1363 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1362 += tmp_1363;
   tmp_1361 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1362;
   std::complex<double> tmp_1364;
   std::complex<double> tmp_1365;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1365 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1364 += tmp_1365;
   tmp_1361 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1364;
   std::complex<double> tmp_1366;
   std::complex<double> tmp_1367;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1367 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1366 += tmp_1367;
   tmp_1361 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1366;
   std::complex<double> tmp_1368;
   std::complex<double> tmp_1369;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1369 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1368 += tmp_1369;
   tmp_1361 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1368;
   std::complex<double> tmp_1370;
   std::complex<double> tmp_1371;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1371 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1370 += tmp_1371;
   tmp_1361 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1370;
   std::complex<double> tmp_1372;
   std::complex<double> tmp_1373;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1373 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1372 += tmp_1373;
   tmp_1361 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1372;
   std::complex<double> tmp_1374;
   std::complex<double> tmp_1375;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1376;
      std::complex<double> tmp_1377;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1377 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1376 += tmp_1377;
      tmp_1375 += (Conj(ZU(gI2,j2))) * tmp_1376;
   }
   tmp_1374 += tmp_1375;
   tmp_1361 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(0,gO2)) * tmp_1374;
   std::complex<double> tmp_1378;
   std::complex<double> tmp_1379;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1380;
      std::complex<double> tmp_1381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1381 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1380 += tmp_1381;
      tmp_1379 += (Conj(ZU(gI2,j2))) * tmp_1380;
   }
   tmp_1378 += tmp_1379;
   tmp_1361 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1378;
   std::complex<double> tmp_1382;
   std::complex<double> tmp_1383;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1384;
      std::complex<double> tmp_1385;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1385 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1384 += tmp_1385;
      tmp_1383 += (ZU(gI1,j2)) * tmp_1384;
   }
   tmp_1382 += tmp_1383;
   tmp_1361 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,
      gO2)*Mu) * tmp_1382;
   std::complex<double> tmp_1386;
   std::complex<double> tmp_1387;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1388;
      std::complex<double> tmp_1389;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1389 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1388 += tmp_1389;
      tmp_1387 += (ZU(gI1,j2)) * tmp_1388;
   }
   tmp_1386 += tmp_1387;
   tmp_1361 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1386;
   std::complex<double> tmp_1390;
   std::complex<double> tmp_1391;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1392;
      std::complex<double> tmp_1393;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1394;
         std::complex<double> tmp_1395;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1395 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1394 += tmp_1395;
         tmp_1393 += (ZU(gI1,3 + j2)) * tmp_1394;
      }
      tmp_1392 += tmp_1393;
      tmp_1391 += (Conj(ZU(gI2,3 + j3))) * tmp_1392;
   }
   tmp_1390 += tmp_1391;
   tmp_1361 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1390
      ;
   std::complex<double> tmp_1396;
   std::complex<double> tmp_1397;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1398;
      std::complex<double> tmp_1399;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1400;
         std::complex<double> tmp_1401;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1401 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1400 += tmp_1401;
         tmp_1399 += (Conj(ZU(gI2,j2))) * tmp_1400;
      }
      tmp_1398 += tmp_1399;
      tmp_1397 += (ZU(gI1,j3)) * tmp_1398;
   }
   tmp_1396 += tmp_1397;
   tmp_1361 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1396
      ;
   result += (std::complex<double>(0,-1)) * tmp_1361;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1402;
   std::complex<double> tmp_1403;
   std::complex<double> tmp_1404;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1404 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1403 += tmp_1404;
   tmp_1402 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1403;
   std::complex<double> tmp_1405;
   std::complex<double> tmp_1406;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1406 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1405 += tmp_1406;
   tmp_1402 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1405;
   std::complex<double> tmp_1407;
   std::complex<double> tmp_1408;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1408 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1407 += tmp_1408;
   tmp_1402 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1407;
   std::complex<double> tmp_1409;
   std::complex<double> tmp_1410;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1410 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1409 += tmp_1410;
   tmp_1402 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1409;
   std::complex<double> tmp_1411;
   std::complex<double> tmp_1412;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1413;
      std::complex<double> tmp_1414;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1414 += Yv(j1,j2)*ZV(gI1,3 + j1);
      }
      tmp_1413 += tmp_1414;
      tmp_1412 += (Conj(ZV(gI2,j2))) * tmp_1413;
   }
   tmp_1411 += tmp_1412;
   tmp_1402 += (std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(0,gO2)) * tmp_1411;
   std::complex<double> tmp_1415;
   std::complex<double> tmp_1416;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1417;
      std::complex<double> tmp_1418;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1418 += ZV(gI1,3 + j1)*TYv(j1,j2);
      }
      tmp_1417 += tmp_1418;
      tmp_1416 += (Conj(ZV(gI2,j2))) * tmp_1417;
   }
   tmp_1415 += tmp_1416;
   tmp_1402 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1415;
   std::complex<double> tmp_1419;
   std::complex<double> tmp_1420;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1421;
      std::complex<double> tmp_1422;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1422 += Conj(Yv(j1,j2))*Conj(ZV(gI2,3 + j1));
      }
      tmp_1421 += tmp_1422;
      tmp_1420 += (ZV(gI1,j2)) * tmp_1421;
   }
   tmp_1419 += tmp_1420;
   tmp_1402 += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,
      gO2)*Mu) * tmp_1419;
   std::complex<double> tmp_1423;
   std::complex<double> tmp_1424;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1425;
      std::complex<double> tmp_1426;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1426 += Conj(ZV(gI2,3 + j1))*Conj(TYv(j1,j2));
      }
      tmp_1425 += tmp_1426;
      tmp_1424 += (ZV(gI1,j2)) * tmp_1425;
   }
   tmp_1423 += tmp_1424;
   tmp_1402 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1423;
   std::complex<double> tmp_1427;
   std::complex<double> tmp_1428;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1429;
      std::complex<double> tmp_1430;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1431;
         std::complex<double> tmp_1432;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1432 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_1431 += tmp_1432;
         tmp_1430 += (ZV(gI1,3 + j2)) * tmp_1431;
      }
      tmp_1429 += tmp_1430;
      tmp_1428 += (Conj(ZV(gI2,3 + j3))) * tmp_1429;
   }
   tmp_1427 += tmp_1428;
   tmp_1402 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1427
      ;
   std::complex<double> tmp_1433;
   std::complex<double> tmp_1434;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1435;
      std::complex<double> tmp_1436;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1437;
         std::complex<double> tmp_1438;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1438 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_1437 += tmp_1438;
         tmp_1436 += (Conj(ZV(gI2,j2))) * tmp_1437;
      }
      tmp_1435 += tmp_1436;
      tmp_1434 += (ZV(gI1,j3)) * tmp_1435;
   }
   tmp_1433 += tmp_1434;
   tmp_1402 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1433
      ;
   result += (std::complex<double>(0,-1)) * tmp_1402;

   return result;
}

std::complex<double> CLASSNAME::CpUhhSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1439;
   std::complex<double> tmp_1440;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1441;
      std::complex<double> tmp_1442;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1443;
         std::complex<double> tmp_1444;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1444 += Conj(Mv(j1,j3))*Yv(j1,j2);
         }
         tmp_1443 += tmp_1444;
         tmp_1442 += (Conj(ZV(gI1,j2))) * tmp_1443;
      }
      tmp_1441 += tmp_1442;
      tmp_1440 += (Conj(ZV(gI2,3 + j3))) * tmp_1441;
   }
   tmp_1439 += tmp_1440;
   std::complex<double> tmp_1445;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1446;
      std::complex<double> tmp_1447;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1448;
         std::complex<double> tmp_1449;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1449 += Conj(Mv(j1,j3))*Yv(j1,j2);
         }
         tmp_1448 += tmp_1449;
         tmp_1447 += (Conj(ZV(gI2,j2))) * tmp_1448;
      }
      tmp_1446 += tmp_1447;
      tmp_1445 += (Conj(ZV(gI1,3 + j3))) * tmp_1446;
   }
   tmp_1439 += tmp_1445;
   std::complex<double> tmp_1450;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1451;
      std::complex<double> tmp_1452;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1453;
         std::complex<double> tmp_1454;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1454 += Conj(Mv(j3,j1))*Yv(j1,j2);
         }
         tmp_1453 += tmp_1454;
         tmp_1452 += (Conj(ZV(gI1,j2))) * tmp_1453;
      }
      tmp_1451 += tmp_1452;
      tmp_1450 += (Conj(ZV(gI2,3 + j3))) * tmp_1451;
   }
   tmp_1439 += tmp_1450;
   std::complex<double> tmp_1455;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1456;
      std::complex<double> tmp_1457;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1458;
         std::complex<double> tmp_1459;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1459 += Conj(Mv(j3,j1))*Yv(j1,j2);
         }
         tmp_1458 += tmp_1459;
         tmp_1457 += (Conj(ZV(gI2,j2))) * tmp_1458;
      }
      tmp_1456 += tmp_1457;
      tmp_1455 += (Conj(ZV(gI1,3 + j3))) * tmp_1456;
   }
   tmp_1439 += tmp_1455;
   result += (-0.35355339059327373*KroneckerDelta(1,gO2)) * tmp_1439;

   return result;
}

std::complex<double> CLASSNAME::CpUhhFvFvPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1460;
   std::complex<double> tmp_1461;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1462;
      std::complex<double> tmp_1463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1463 += Conj(Yv(j1,j2))*UV(gI2,3 + j1);
      }
      tmp_1462 += tmp_1463;
      tmp_1461 += (UV(gI1,j2)) * tmp_1462;
   }
   tmp_1460 += tmp_1461;
   std::complex<double> tmp_1464;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1465;
      std::complex<double> tmp_1466;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1466 += Conj(Yv(j1,j2))*UV(gI1,3 + j1);
      }
      tmp_1465 += tmp_1466;
      tmp_1464 += (UV(gI2,j2)) * tmp_1465;
   }
   tmp_1460 += tmp_1464;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1460;

   return result;
}

std::complex<double> CLASSNAME::CpUhhFvFvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1467;
   std::complex<double> tmp_1468;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1469;
      std::complex<double> tmp_1470;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1470 += Conj(UV(gI1,3 + j1))*Yv(j1,j2);
      }
      tmp_1469 += tmp_1470;
      tmp_1468 += (Conj(UV(gI2,j2))) * tmp_1469;
   }
   tmp_1467 += tmp_1468;
   std::complex<double> tmp_1471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1472;
      std::complex<double> tmp_1473;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1473 += Conj(UV(gI2,3 + j1))*Yv(j1,j2);
      }
      tmp_1472 += tmp_1473;
      tmp_1471 += (Conj(UV(gI1,j2))) * tmp_1472;
   }
   tmp_1467 += tmp_1471;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1467;

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
      gI2,1));

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
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(1,gO1)*(KroneckerDelta
      (0,gO2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*
      (ZA(gI1,0)*ZA(gI2,0) - 3*ZA(gI1,1)*ZA(gI2,1))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(1,gO2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) +
      KroneckerDelta(0,gO2)*(-3*ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(
      gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) - KroneckerDelta(0,gO2)*((3*Sqr(g1)
      + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2
      ,1))) + KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP
      (gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(
      g2))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,0)*ZH(gI2,0)
      - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vu*KroneckerDelta(0,gO2) + vd*
      KroneckerDelta(1,gO2))*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(0,gO2)*ZA(gI2,0) -
      KroneckerDelta(1,gO2)*ZA(gI2,1))*(vd*ZH(gI1,0) - vu*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*g2*(KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) + KroneckerDelta(1,gO2)*UM(gI1,0)*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gI2,1))*
      Conj(UP(gI1,0))*KroneckerDelta(0,gO1) + Conj(UM(gI2,0))*Conj(UP(gI1,1))*
      KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1474;
   std::complex<double> tmp_1475;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1476;
      std::complex<double> tmp_1477;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1477 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1476 += tmp_1477;
      tmp_1475 += (ZDL(gI1,j2)) * tmp_1476;
   }
   tmp_1474 += tmp_1475;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1474;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1478;
   std::complex<double> tmp_1479;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1480;
      std::complex<double> tmp_1481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1481 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1480 += tmp_1481;
      tmp_1479 += (Conj(ZDL(gI2,j2))) * tmp_1480;
   }
   tmp_1478 += tmp_1479;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1478;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1482;
   std::complex<double> tmp_1483;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1484;
      std::complex<double> tmp_1485;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1485 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1484 += tmp_1485;
      tmp_1483 += (ZEL(gI1,j2)) * tmp_1484;
   }
   tmp_1482 += tmp_1483;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1482;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1486;
   std::complex<double> tmp_1487;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1488;
      std::complex<double> tmp_1489;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1489 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1488 += tmp_1489;
      tmp_1487 += (Conj(ZEL(gI2,j2))) * tmp_1488;
   }
   tmp_1486 += tmp_1487;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1486;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1490;
   std::complex<double> tmp_1491;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1492;
      std::complex<double> tmp_1493;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1493 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1492 += tmp_1493;
      tmp_1491 += (ZUL(gI1,j2)) * tmp_1492;
   }
   tmp_1490 += tmp_1491;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_1490;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1494;
   std::complex<double> tmp_1495;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1496;
      std::complex<double> tmp_1497;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1497 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1496 += tmp_1497;
      tmp_1495 += (Conj(ZUL(gI2,j2))) * tmp_1496;
   }
   tmp_1494 += tmp_1495;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_1494;

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(KroneckerDelta(0,gO2)*(ZN(gI1,2)*(
      3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (3.872983346207417*g1*ZN(
      gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2)) - KroneckerDelta(1,gO2)*(ZN(gI1,3)*(
      3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (3.872983346207417*g1*ZN(
      gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(Conj(ZN(gI1,2))*(3.872983346207417*g1
      *Conj(ZN(gI2,0)) - 5*g2*Conj(ZN(gI2,1)))*KroneckerDelta(0,gO1) - 5*g2*Conj(
      ZN(gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) - 3.872983346207417*g1*Conj
      (ZN(gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,3))*
      Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))
      *KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(Conj(ZN(gI2,2
      ))*KroneckerDelta(0,gO1) - Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1498;
   std::complex<double> tmp_1499;
   std::complex<double> tmp_1500;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1500 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1499 += tmp_1500;
   tmp_1498 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1499;
   std::complex<double> tmp_1501;
   std::complex<double> tmp_1502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1502 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1501 += tmp_1502;
   tmp_1498 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1501;
   std::complex<double> tmp_1503;
   std::complex<double> tmp_1504;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1504 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1503 += tmp_1504;
   tmp_1498 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1503;
   std::complex<double> tmp_1505;
   std::complex<double> tmp_1506;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1506 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1505 += tmp_1506;
   tmp_1498 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1505;
   std::complex<double> tmp_1507;
   std::complex<double> tmp_1508;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1508 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1507 += tmp_1508;
   tmp_1498 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1507;
   std::complex<double> tmp_1509;
   std::complex<double> tmp_1510;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1510 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1509 += tmp_1510;
   tmp_1498 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1509;
   std::complex<double> tmp_1511;
   std::complex<double> tmp_1512;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1513;
      std::complex<double> tmp_1514;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1515;
         std::complex<double> tmp_1516;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1516 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1515 += tmp_1516;
         tmp_1514 += (ZD(gI1,3 + j2)) * tmp_1515;
      }
      tmp_1513 += tmp_1514;
      tmp_1512 += (Conj(ZD(gI2,3 + j3))) * tmp_1513;
   }
   tmp_1511 += tmp_1512;
   tmp_1498 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1511;
   std::complex<double> tmp_1517;
   std::complex<double> tmp_1518;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1519;
      std::complex<double> tmp_1520;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1521;
         std::complex<double> tmp_1522;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1522 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1521 += tmp_1522;
         tmp_1520 += (Conj(ZD(gI2,j2))) * tmp_1521;
      }
      tmp_1519 += tmp_1520;
      tmp_1518 += (ZD(gI1,j3)) * tmp_1519;
   }
   tmp_1517 += tmp_1518;
   tmp_1498 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1517;
   result += (std::complex<double>(0,-1)) * tmp_1498;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1523;
   std::complex<double> tmp_1524;
   std::complex<double> tmp_1525;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1525 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1524 += tmp_1525;
   tmp_1523 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1524;
   std::complex<double> tmp_1526;
   std::complex<double> tmp_1527;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1527 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1526 += tmp_1527;
   tmp_1523 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1526;
   std::complex<double> tmp_1528;
   std::complex<double> tmp_1529;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1529 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1528 += tmp_1529;
   tmp_1523 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1528;
   std::complex<double> tmp_1530;
   std::complex<double> tmp_1531;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1531 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1530 += tmp_1531;
   tmp_1523 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1530;
   std::complex<double> tmp_1532;
   std::complex<double> tmp_1533;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1533 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1532 += tmp_1533;
   tmp_1523 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1532;
   std::complex<double> tmp_1534;
   std::complex<double> tmp_1535;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1535 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1534 += tmp_1535;
   tmp_1523 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1534;
   std::complex<double> tmp_1536;
   std::complex<double> tmp_1537;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1538;
      std::complex<double> tmp_1539;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1540;
         std::complex<double> tmp_1541;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1541 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1540 += tmp_1541;
         tmp_1539 += (ZE(gI1,3 + j2)) * tmp_1540;
      }
      tmp_1538 += tmp_1539;
      tmp_1537 += (Conj(ZE(gI2,3 + j3))) * tmp_1538;
   }
   tmp_1536 += tmp_1537;
   tmp_1523 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1536;
   std::complex<double> tmp_1542;
   std::complex<double> tmp_1543;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1544;
      std::complex<double> tmp_1545;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1546;
         std::complex<double> tmp_1547;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1547 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1546 += tmp_1547;
         tmp_1545 += (Conj(ZE(gI2,j2))) * tmp_1546;
      }
      tmp_1544 += tmp_1545;
      tmp_1543 += (ZE(gI1,j3)) * tmp_1544;
   }
   tmp_1542 += tmp_1543;
   tmp_1523 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1542;
   result += (std::complex<double>(0,-1)) * tmp_1523;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1548;
   std::complex<double> tmp_1549;
   std::complex<double> tmp_1550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1550 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1549 += tmp_1550;
   tmp_1548 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1549;
   std::complex<double> tmp_1551;
   std::complex<double> tmp_1552;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1552 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1551 += tmp_1552;
   tmp_1548 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1551;
   std::complex<double> tmp_1553;
   std::complex<double> tmp_1554;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1554 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1553 += tmp_1554;
   tmp_1548 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1553;
   std::complex<double> tmp_1555;
   std::complex<double> tmp_1556;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1556 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1555 += tmp_1556;
   tmp_1548 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1555;
   std::complex<double> tmp_1557;
   std::complex<double> tmp_1558;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1558 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1557 += tmp_1558;
   tmp_1548 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1557;
   std::complex<double> tmp_1559;
   std::complex<double> tmp_1560;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1560 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1559 += tmp_1560;
   tmp_1548 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1559;
   std::complex<double> tmp_1561;
   std::complex<double> tmp_1562;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1563;
      std::complex<double> tmp_1564;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1565;
         std::complex<double> tmp_1566;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1566 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1565 += tmp_1566;
         tmp_1564 += (ZU(gI1,3 + j2)) * tmp_1565;
      }
      tmp_1563 += tmp_1564;
      tmp_1562 += (Conj(ZU(gI2,3 + j3))) * tmp_1563;
   }
   tmp_1561 += tmp_1562;
   tmp_1548 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1561;
   std::complex<double> tmp_1567;
   std::complex<double> tmp_1568;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1569;
      std::complex<double> tmp_1570;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1571;
         std::complex<double> tmp_1572;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1572 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1571 += tmp_1572;
         tmp_1570 += (Conj(ZU(gI2,j2))) * tmp_1571;
      }
      tmp_1569 += tmp_1570;
      tmp_1568 += (ZU(gI1,j3)) * tmp_1569;
   }
   tmp_1567 += tmp_1568;
   tmp_1548 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1567;
   result += (std::complex<double>(0,-1)) * tmp_1548;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1573;
   std::complex<double> tmp_1574;
   std::complex<double> tmp_1575;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1575 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1574 += tmp_1575;
   tmp_1573 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1574;
   std::complex<double> tmp_1576;
   std::complex<double> tmp_1577;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1577 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1576 += tmp_1577;
   tmp_1573 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1576;
   std::complex<double> tmp_1578;
   std::complex<double> tmp_1579;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1579 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1578 += tmp_1579;
   tmp_1573 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1578;
   std::complex<double> tmp_1580;
   std::complex<double> tmp_1581;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1581 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1580 += tmp_1581;
   tmp_1573 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1580;
   std::complex<double> tmp_1582;
   std::complex<double> tmp_1583;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1584;
      std::complex<double> tmp_1585;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1586;
         std::complex<double> tmp_1587;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1587 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_1586 += tmp_1587;
         tmp_1585 += (ZV(gI1,3 + j2)) * tmp_1586;
      }
      tmp_1584 += tmp_1585;
      tmp_1583 += (Conj(ZV(gI2,3 + j3))) * tmp_1584;
   }
   tmp_1582 += tmp_1583;
   tmp_1573 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1582;
   std::complex<double> tmp_1588;
   std::complex<double> tmp_1589;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1590;
      std::complex<double> tmp_1591;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1592;
         std::complex<double> tmp_1593;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1593 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_1592 += tmp_1593;
         tmp_1591 += (Conj(ZV(gI2,j2))) * tmp_1592;
      }
      tmp_1590 += tmp_1591;
      tmp_1589 += (ZV(gI1,j3)) * tmp_1590;
   }
   tmp_1588 += tmp_1589;
   tmp_1573 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1588;
   result += (std::complex<double>(0,-1)) * tmp_1573;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1594;
   std::complex<double> tmp_1595;
   std::complex<double> tmp_1596;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1597;
      std::complex<double> tmp_1598;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1598 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1597 += tmp_1598;
      tmp_1596 += (Conj(ZD(gI2,j2))) * tmp_1597;
   }
   tmp_1595 += tmp_1596;
   tmp_1594 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1595;
   std::complex<double> tmp_1599;
   std::complex<double> tmp_1600;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1601;
      std::complex<double> tmp_1602;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1602 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1601 += tmp_1602;
      tmp_1600 += (Conj(ZD(gI2,j2))) * tmp_1601;
   }
   tmp_1599 += tmp_1600;
   tmp_1594 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1599;
   std::complex<double> tmp_1603;
   std::complex<double> tmp_1604;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1605;
      std::complex<double> tmp_1606;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1606 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1605 += tmp_1606;
      tmp_1604 += (ZD(gI1,j2)) * tmp_1605;
   }
   tmp_1603 += tmp_1604;
   tmp_1594 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1603;
   std::complex<double> tmp_1607;
   std::complex<double> tmp_1608;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1609;
      std::complex<double> tmp_1610;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1610 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1609 += tmp_1610;
      tmp_1608 += (ZD(gI1,j2)) * tmp_1609;
   }
   tmp_1607 += tmp_1608;
   tmp_1594 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1607;
   result += (std::complex<double>(0,-1)) * tmp_1594;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1611;
   std::complex<double> tmp_1612;
   std::complex<double> tmp_1613;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1614;
      std::complex<double> tmp_1615;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1615 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1614 += tmp_1615;
      tmp_1613 += (Conj(ZE(gI2,j2))) * tmp_1614;
   }
   tmp_1612 += tmp_1613;
   tmp_1611 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1612;
   std::complex<double> tmp_1616;
   std::complex<double> tmp_1617;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1618;
      std::complex<double> tmp_1619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1619 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1618 += tmp_1619;
      tmp_1617 += (Conj(ZE(gI2,j2))) * tmp_1618;
   }
   tmp_1616 += tmp_1617;
   tmp_1611 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1616;
   std::complex<double> tmp_1620;
   std::complex<double> tmp_1621;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1622;
      std::complex<double> tmp_1623;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1623 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1622 += tmp_1623;
      tmp_1621 += (ZE(gI1,j2)) * tmp_1622;
   }
   tmp_1620 += tmp_1621;
   tmp_1611 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1620;
   std::complex<double> tmp_1624;
   std::complex<double> tmp_1625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1626;
      std::complex<double> tmp_1627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1627 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1626 += tmp_1627;
      tmp_1625 += (ZE(gI1,j2)) * tmp_1626;
   }
   tmp_1624 += tmp_1625;
   tmp_1611 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1624;
   result += (std::complex<double>(0,-1)) * tmp_1611;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1628;
   std::complex<double> tmp_1629;
   std::complex<double> tmp_1630;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1631;
      std::complex<double> tmp_1632;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1632 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1631 += tmp_1632;
      tmp_1630 += (Conj(ZU(gI2,j2))) * tmp_1631;
   }
   tmp_1629 += tmp_1630;
   tmp_1628 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(0,gO2)) * tmp_1629;
   std::complex<double> tmp_1633;
   std::complex<double> tmp_1634;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1635;
      std::complex<double> tmp_1636;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1636 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1635 += tmp_1636;
      tmp_1634 += (Conj(ZU(gI2,j2))) * tmp_1635;
   }
   tmp_1633 += tmp_1634;
   tmp_1628 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1633;
   std::complex<double> tmp_1637;
   std::complex<double> tmp_1638;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1639;
      std::complex<double> tmp_1640;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1640 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1639 += tmp_1640;
      tmp_1638 += (ZU(gI1,j2)) * tmp_1639;
   }
   tmp_1637 += tmp_1638;
   tmp_1628 += (-0.7071067811865475*KroneckerDelta(0,gO2)*Mu) * tmp_1637;
   std::complex<double> tmp_1641;
   std::complex<double> tmp_1642;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1643;
      std::complex<double> tmp_1644;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1644 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1643 += tmp_1644;
      tmp_1642 += (ZU(gI1,j2)) * tmp_1643;
   }
   tmp_1641 += tmp_1642;
   tmp_1628 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1641;
   result += (std::complex<double>(0,-1)) * tmp_1628;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1645;
   std::complex<double> tmp_1646;
   std::complex<double> tmp_1647;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1648;
      std::complex<double> tmp_1649;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1649 += Yv(j1,j2)*ZV(gI1,3 + j1);
      }
      tmp_1648 += tmp_1649;
      tmp_1647 += (Conj(ZV(gI2,j2))) * tmp_1648;
   }
   tmp_1646 += tmp_1647;
   tmp_1645 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(0,gO2)) * tmp_1646;
   std::complex<double> tmp_1650;
   std::complex<double> tmp_1651;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1652;
      std::complex<double> tmp_1653;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1653 += ZV(gI1,3 + j1)*TYv(j1,j2);
      }
      tmp_1652 += tmp_1653;
      tmp_1651 += (Conj(ZV(gI2,j2))) * tmp_1652;
   }
   tmp_1650 += tmp_1651;
   tmp_1645 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1650;
   std::complex<double> tmp_1654;
   std::complex<double> tmp_1655;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1656;
      std::complex<double> tmp_1657;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1657 += Conj(Yv(j1,j2))*Conj(ZV(gI2,3 + j1));
      }
      tmp_1656 += tmp_1657;
      tmp_1655 += (ZV(gI1,j2)) * tmp_1656;
   }
   tmp_1654 += tmp_1655;
   tmp_1645 += (-0.7071067811865475*KroneckerDelta(0,gO2)*Mu) * tmp_1654;
   std::complex<double> tmp_1658;
   std::complex<double> tmp_1659;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1660;
      std::complex<double> tmp_1661;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1661 += Conj(ZV(gI2,3 + j1))*Conj(TYv(j1,j2));
      }
      tmp_1660 += tmp_1661;
      tmp_1659 += (ZV(gI1,j2)) * tmp_1660;
   }
   tmp_1658 += tmp_1659;
   tmp_1645 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1658;
   result += (std::complex<double>(0,-1)) * tmp_1645;

   return result;
}

std::complex<double> CLASSNAME::CpUAhSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1662;
   std::complex<double> tmp_1663;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1664;
      std::complex<double> tmp_1665;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1666;
         std::complex<double> tmp_1667;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1667 += Conj(Mv(j1,j3))*Yv(j1,j2);
         }
         tmp_1666 += tmp_1667;
         tmp_1665 += (Conj(ZV(gI1,j2))) * tmp_1666;
      }
      tmp_1664 += tmp_1665;
      tmp_1663 += (Conj(ZV(gI2,3 + j3))) * tmp_1664;
   }
   tmp_1662 += tmp_1663;
   std::complex<double> tmp_1668;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1669;
      std::complex<double> tmp_1670;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1671;
         std::complex<double> tmp_1672;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1672 += Conj(Mv(j1,j3))*Yv(j1,j2);
         }
         tmp_1671 += tmp_1672;
         tmp_1670 += (Conj(ZV(gI2,j2))) * tmp_1671;
      }
      tmp_1669 += tmp_1670;
      tmp_1668 += (Conj(ZV(gI1,3 + j3))) * tmp_1669;
   }
   tmp_1662 += tmp_1668;
   std::complex<double> tmp_1673;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1674;
      std::complex<double> tmp_1675;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1676;
         std::complex<double> tmp_1677;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1677 += Conj(Mv(j3,j1))*Yv(j1,j2);
         }
         tmp_1676 += tmp_1677;
         tmp_1675 += (Conj(ZV(gI1,j2))) * tmp_1676;
      }
      tmp_1674 += tmp_1675;
      tmp_1673 += (Conj(ZV(gI2,3 + j3))) * tmp_1674;
   }
   tmp_1662 += tmp_1673;
   std::complex<double> tmp_1678;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1679;
      std::complex<double> tmp_1680;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1681;
         std::complex<double> tmp_1682;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1682 += Conj(Mv(j3,j1))*Yv(j1,j2);
         }
         tmp_1681 += tmp_1682;
         tmp_1680 += (Conj(ZV(gI2,j2))) * tmp_1681;
      }
      tmp_1679 += tmp_1680;
      tmp_1678 += (Conj(ZV(gI1,3 + j3))) * tmp_1679;
   }
   tmp_1662 += tmp_1678;
   result += (std::complex<double>(0.,-0.35355339059327373)*KroneckerDelta(1,
      gO2)) * tmp_1662;

   return result;
}

std::complex<double> CLASSNAME::CpUAhFvFvPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1683;
   std::complex<double> tmp_1684;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1685;
      std::complex<double> tmp_1686;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1686 += Conj(Yv(j1,j2))*UV(gI2,3 + j1);
      }
      tmp_1685 += tmp_1686;
      tmp_1684 += (UV(gI1,j2)) * tmp_1685;
   }
   tmp_1683 += tmp_1684;
   std::complex<double> tmp_1687;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1688;
      std::complex<double> tmp_1689;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1689 += Conj(Yv(j1,j2))*UV(gI1,3 + j1);
      }
      tmp_1688 += tmp_1689;
      tmp_1687 += (UV(gI2,j2)) * tmp_1688;
   }
   tmp_1683 += tmp_1687;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_1683;

   return result;
}

std::complex<double> CLASSNAME::CpUAhFvFvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1692;
      std::complex<double> tmp_1693;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1693 += Conj(UV(gI1,3 + j1))*Yv(j1,j2);
      }
      tmp_1692 += tmp_1693;
      tmp_1691 += (Conj(UV(gI2,j2))) * tmp_1692;
   }
   tmp_1690 += tmp_1691;
   std::complex<double> tmp_1694;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1695;
      std::complex<double> tmp_1696;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1696 += Conj(UV(gI2,3 + j1))*Yv(j1,j2);
      }
      tmp_1695 += tmp_1696;
      tmp_1694 += (Conj(UV(gI1,j2))) * tmp_1695;
   }
   tmp_1690 += tmp_1694;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_1690;

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
      KroneckerDelta(1,gO2)*ZP(gI2,1));

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
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*(vd*KroneckerDelta(0,gO2) - vu*
      KroneckerDelta(1,gO2))*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargWmCgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*(5*g2
      *Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

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

   result = 0.05*g2*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*(5*g2
      *Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*Sqr(g2)*(ZA(
      gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) - KroneckerDelta(0,gO2)*((3*Sqr(g1)
      + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2
      ,1))) + KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*Sqr(g2)*(ZA(gI1,1)*ZA
      (gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(
      g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(1,gO1)*(KroneckerDelta
      (0,gO2)*ZP(gI1,1)*ZP(gI2,0) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*ZP(gI2,0) - 2
      *ZP(gI1,1)*ZP(gI2,1))) + KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*ZP(gI1
      ,0)*ZP(gI2,1) + KroneckerDelta(0,gO2)*(-2*ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP
      (gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*Sqr(g2)*(ZH(
      gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1)
      + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2
      ,1)))) + KroneckerDelta(1,gO1)*(-5*KroneckerDelta(0,gO2)*Sqr(g2)*(ZH(gI1,1)*
      ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr
      (g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*Sqr(g2)*(vu*ZA(gI2,0) + vd*ZA(gI2,1))*
      (KroneckerDelta(1,gO2)*ZP(gI1,0) - KroneckerDelta(0,gO2)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(ZH(gI2,1)*((-3*vu*Sqr(g1) + 5*vu*Sqr
      (g2))*ZP(gI1,0) + 5*vd*Sqr(g2)*ZP(gI1,1)) + ZH(gI2,0)*(vd*(3*Sqr(g1) + 5*Sqr
      (g2))*ZP(gI1,0) + 5*vu*Sqr(g2)*ZP(gI1,1)))) - KroneckerDelta(1,gO2)*(ZH(gI2,
      0)*(5*vu*Sqr(g2)*ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)) + ZH(gI2
      ,1)*(5*vd*Sqr(g2)*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1697;
   std::complex<double> tmp_1698;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1699;
      std::complex<double> tmp_1700;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1700 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1699 += tmp_1700;
      tmp_1698 += (ZUL(gI1,j2)) * tmp_1699;
   }
   tmp_1697 += tmp_1698;
   result += (KroneckerDelta(0,gO2)) * tmp_1697;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
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
      tmp_1702 += (Conj(ZDL(gI2,j2))) * tmp_1703;
   }
   tmp_1701 += tmp_1702;
   result += (KroneckerDelta(1,gO1)) * tmp_1701;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*KroneckerDelta(1,gO2)*(UP(gI2,1)*(1.0954451150103321*g1*ZN(gI1
      ,0) + 1.4142135623730951*g2*ZN(gI1,1)) + 2*g2*UP(gI2,0)*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(UM(gI2,1))*(1.0954451150103321*g1*Conj(ZN(gI1,0)) +
      1.4142135623730951*g2*Conj(ZN(gI1,1))) - 2*g2*Conj(UM(gI2,0))*Conj(ZN(gI1,2)
      ))*KroneckerDelta(0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1705;
   std::complex<double> tmp_1706;
   std::complex<double> tmp_1707;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1707 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1706 += tmp_1707;
   tmp_1705 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1706;
   std::complex<double> tmp_1708;
   std::complex<double> tmp_1709;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1709 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1708 += tmp_1709;
   tmp_1705 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1708;
   std::complex<double> tmp_1710;
   std::complex<double> tmp_1711;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1711 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1710 += tmp_1711;
   tmp_1705 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1710;
   std::complex<double> tmp_1712;
   std::complex<double> tmp_1713;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1713 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1712 += tmp_1713;
   tmp_1705 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1712;
   std::complex<double> tmp_1714;
   std::complex<double> tmp_1715;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1715 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1714 += tmp_1715;
   tmp_1705 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1714;
   std::complex<double> tmp_1716;
   std::complex<double> tmp_1717;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1717 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1716 += tmp_1717;
   tmp_1705 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1716;
   std::complex<double> tmp_1718;
   std::complex<double> tmp_1719;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1720;
      std::complex<double> tmp_1721;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1722;
         std::complex<double> tmp_1723;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1723 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1722 += tmp_1723;
         tmp_1721 += (ZD(gI1,3 + j2)) * tmp_1722;
      }
      tmp_1720 += tmp_1721;
      tmp_1719 += (Conj(ZD(gI2,3 + j3))) * tmp_1720;
   }
   tmp_1718 += tmp_1719;
   tmp_1705 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1718;
   std::complex<double> tmp_1724;
   std::complex<double> tmp_1725;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1726;
      std::complex<double> tmp_1727;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1728;
         std::complex<double> tmp_1729;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1729 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1728 += tmp_1729;
         tmp_1727 += (Conj(ZD(gI2,j2))) * tmp_1728;
      }
      tmp_1726 += tmp_1727;
      tmp_1725 += (ZD(gI1,j3)) * tmp_1726;
   }
   tmp_1724 += tmp_1725;
   tmp_1705 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1724;
   result += (std::complex<double>(0,-1)) * tmp_1705;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1730;
   std::complex<double> tmp_1731;
   std::complex<double> tmp_1732;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1732 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1731 += tmp_1732;
   tmp_1730 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1731;
   std::complex<double> tmp_1733;
   std::complex<double> tmp_1734;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1734 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1733 += tmp_1734;
   tmp_1730 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1733;
   std::complex<double> tmp_1735;
   std::complex<double> tmp_1736;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1736 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1735 += tmp_1736;
   tmp_1730 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1735;
   std::complex<double> tmp_1737;
   std::complex<double> tmp_1738;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1738 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1737 += tmp_1738;
   tmp_1730 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1737;
   std::complex<double> tmp_1739;
   std::complex<double> tmp_1740;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1740 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1739 += tmp_1740;
   tmp_1730 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1739;
   std::complex<double> tmp_1741;
   std::complex<double> tmp_1742;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1742 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1741 += tmp_1742;
   tmp_1730 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1741;
   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1745;
      std::complex<double> tmp_1746;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1747;
         std::complex<double> tmp_1748;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1748 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1747 += tmp_1748;
         tmp_1746 += (ZE(gI1,3 + j2)) * tmp_1747;
      }
      tmp_1745 += tmp_1746;
      tmp_1744 += (Conj(ZE(gI2,3 + j3))) * tmp_1745;
   }
   tmp_1743 += tmp_1744;
   tmp_1730 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1743;
   std::complex<double> tmp_1749;
   std::complex<double> tmp_1750;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1751;
      std::complex<double> tmp_1752;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1753;
         std::complex<double> tmp_1754;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1754 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_1753 += tmp_1754;
         tmp_1752 += (Conj(ZE(gI2,j2))) * tmp_1753;
      }
      tmp_1751 += tmp_1752;
      tmp_1750 += (ZE(gI1,j3)) * tmp_1751;
   }
   tmp_1749 += tmp_1750;
   tmp_1730 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1749;
   result += (std::complex<double>(0,-1)) * tmp_1730;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1755;
   std::complex<double> tmp_1756;
   std::complex<double> tmp_1757;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1757 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1756 += tmp_1757;
   tmp_1755 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1756;
   std::complex<double> tmp_1758;
   std::complex<double> tmp_1759;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1759 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1758 += tmp_1759;
   tmp_1755 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1758;
   std::complex<double> tmp_1760;
   std::complex<double> tmp_1761;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1761 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1760 += tmp_1761;
   tmp_1755 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1760;
   std::complex<double> tmp_1762;
   std::complex<double> tmp_1763;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1763 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1762 += tmp_1763;
   tmp_1755 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1762;
   std::complex<double> tmp_1764;
   std::complex<double> tmp_1765;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1765 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1764 += tmp_1765;
   tmp_1755 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1764;
   std::complex<double> tmp_1766;
   std::complex<double> tmp_1767;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1767 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1766 += tmp_1767;
   tmp_1755 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1766;
   std::complex<double> tmp_1768;
   std::complex<double> tmp_1769;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1770;
      std::complex<double> tmp_1771;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1772;
         std::complex<double> tmp_1773;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1773 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1772 += tmp_1773;
         tmp_1771 += (ZU(gI1,3 + j2)) * tmp_1772;
      }
      tmp_1770 += tmp_1771;
      tmp_1769 += (Conj(ZU(gI2,3 + j3))) * tmp_1770;
   }
   tmp_1768 += tmp_1769;
   tmp_1755 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1768;
   std::complex<double> tmp_1774;
   std::complex<double> tmp_1775;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1776;
      std::complex<double> tmp_1777;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1778;
         std::complex<double> tmp_1779;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1779 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1778 += tmp_1779;
         tmp_1777 += (Conj(ZU(gI2,j2))) * tmp_1778;
      }
      tmp_1776 += tmp_1777;
      tmp_1775 += (ZU(gI1,j3)) * tmp_1776;
   }
   tmp_1774 += tmp_1775;
   tmp_1755 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1774;
   result += (std::complex<double>(0,-1)) * tmp_1755;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1780;
   std::complex<double> tmp_1781;
   std::complex<double> tmp_1782;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1782 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1781 += tmp_1782;
   tmp_1780 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1781;
   std::complex<double> tmp_1783;
   std::complex<double> tmp_1784;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1784 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1783 += tmp_1784;
   tmp_1780 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1783;
   std::complex<double> tmp_1785;
   std::complex<double> tmp_1786;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1786 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1785 += tmp_1786;
   tmp_1780 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1785;
   std::complex<double> tmp_1787;
   std::complex<double> tmp_1788;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1788 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1787 += tmp_1788;
   tmp_1780 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1787;
   std::complex<double> tmp_1789;
   std::complex<double> tmp_1790;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1791;
      std::complex<double> tmp_1792;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1793;
         std::complex<double> tmp_1794;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1794 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_1793 += tmp_1794;
         tmp_1792 += (ZV(gI1,3 + j2)) * tmp_1793;
      }
      tmp_1791 += tmp_1792;
      tmp_1790 += (Conj(ZV(gI2,3 + j3))) * tmp_1791;
   }
   tmp_1789 += tmp_1790;
   tmp_1780 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1789;
   std::complex<double> tmp_1795;
   std::complex<double> tmp_1796;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1797;
      std::complex<double> tmp_1798;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1799;
         std::complex<double> tmp_1800;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1800 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1799 += tmp_1800;
         tmp_1798 += (Conj(ZV(gI2,j2))) * tmp_1799;
      }
      tmp_1797 += tmp_1798;
      tmp_1796 += (ZV(gI1,j3)) * tmp_1797;
   }
   tmp_1795 += tmp_1796;
   tmp_1780 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1795;
   result += (std::complex<double>(0,-1)) * tmp_1780;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1801;
   std::complex<double> tmp_1802;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1803;
      std::complex<double> tmp_1804;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1804 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1803 += tmp_1804;
      tmp_1802 += (UV(gI1,j2)) * tmp_1803;
   }
   tmp_1801 += tmp_1802;
   result += (KroneckerDelta(0,gO2)) * tmp_1801;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmFvFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1805;
   std::complex<double> tmp_1806;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1807;
      std::complex<double> tmp_1808;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1808 += Conj(UV(gI1,3 + j1))*Yv(j1,j2);
      }
      tmp_1807 += tmp_1808;
      tmp_1806 += (Conj(ZEL(gI2,j2))) * tmp_1807;
   }
   tmp_1805 += tmp_1806;
   result += (KroneckerDelta(1,gO1)) * tmp_1805;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1809;
   std::complex<double> tmp_1810;
   std::complex<double> tmp_1811;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1811 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1810 += tmp_1811;
   tmp_1809 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1810;
   std::complex<double> tmp_1812;
   std::complex<double> tmp_1813;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1813 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1812 += tmp_1813;
   tmp_1809 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1812;
   std::complex<double> tmp_1814;
   std::complex<double> tmp_1815;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1816;
      std::complex<double> tmp_1817;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1817 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1816 += tmp_1817;
      tmp_1815 += (Conj(ZD(gI2,j2))) * tmp_1816;
   }
   tmp_1814 += tmp_1815;
   tmp_1809 += (std::complex<double>(0,1)*Conj(Mu)*KroneckerDelta(0,gO2)) *
      tmp_1814;
   std::complex<double> tmp_1818;
   std::complex<double> tmp_1819;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1820;
      std::complex<double> tmp_1821;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1821 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1820 += tmp_1821;
      tmp_1819 += (Conj(ZD(gI2,j2))) * tmp_1820;
   }
   tmp_1818 += tmp_1819;
   tmp_1809 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1818;
   std::complex<double> tmp_1822;
   std::complex<double> tmp_1823;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1824;
      std::complex<double> tmp_1825;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1825 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1824 += tmp_1825;
      tmp_1823 += (ZU(gI1,j2)) * tmp_1824;
   }
   tmp_1822 += tmp_1823;
   tmp_1809 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1822;
   std::complex<double> tmp_1826;
   std::complex<double> tmp_1827;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1828;
      std::complex<double> tmp_1829;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1829 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1828 += tmp_1829;
      tmp_1827 += (ZU(gI1,j2)) * tmp_1828;
   }
   tmp_1826 += tmp_1827;
   tmp_1809 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1826;
   std::complex<double> tmp_1830;
   std::complex<double> tmp_1831;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1832;
      std::complex<double> tmp_1833;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1834;
         std::complex<double> tmp_1835;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1835 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1834 += tmp_1835;
         tmp_1833 += (ZU(gI1,3 + j2)) * tmp_1834;
      }
      tmp_1832 += tmp_1833;
      tmp_1831 += (Conj(ZD(gI2,3 + j3))) * tmp_1832;
   }
   tmp_1830 += tmp_1831;
   tmp_1809 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_1830;
   std::complex<double> tmp_1836;
   std::complex<double> tmp_1837;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1838;
      std::complex<double> tmp_1839;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1840;
         std::complex<double> tmp_1841;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1841 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1840 += tmp_1841;
         tmp_1839 += (ZU(gI1,3 + j2)) * tmp_1840;
      }
      tmp_1838 += tmp_1839;
      tmp_1837 += (Conj(ZD(gI2,3 + j3))) * tmp_1838;
   }
   tmp_1836 += tmp_1837;
   tmp_1809 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_1836;
   std::complex<double> tmp_1842;
   std::complex<double> tmp_1843;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1844;
      std::complex<double> tmp_1845;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1846;
         std::complex<double> tmp_1847;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1847 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1846 += tmp_1847;
         tmp_1845 += (Conj(ZD(gI2,j2))) * tmp_1846;
      }
      tmp_1844 += tmp_1845;
      tmp_1843 += (ZU(gI1,j3)) * tmp_1844;
   }
   tmp_1842 += tmp_1843;
   tmp_1809 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1842;
   std::complex<double> tmp_1848;
   std::complex<double> tmp_1849;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1850;
      std::complex<double> tmp_1851;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1852;
         std::complex<double> tmp_1853;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1853 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1852 += tmp_1853;
         tmp_1851 += (Conj(ZD(gI2,j2))) * tmp_1852;
      }
      tmp_1850 += tmp_1851;
      tmp_1849 += (ZU(gI1,j3)) * tmp_1850;
   }
   tmp_1848 += tmp_1849;
   tmp_1809 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_1848;
   result += (std::complex<double>(0,-1)) * tmp_1809;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1854;
   std::complex<double> tmp_1855;
   std::complex<double> tmp_1856;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1856 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1855 += tmp_1856;
   tmp_1854 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1855;
   std::complex<double> tmp_1857;
   std::complex<double> tmp_1858;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1858 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1857 += tmp_1858;
   tmp_1854 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1857;
   std::complex<double> tmp_1859;
   std::complex<double> tmp_1860;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1861;
      std::complex<double> tmp_1862;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1862 += Yv(j1,j2)*ZV(gI1,3 + j1);
      }
      tmp_1861 += tmp_1862;
      tmp_1860 += (Conj(ZE(gI2,j2))) * tmp_1861;
   }
   tmp_1859 += tmp_1860;
   tmp_1854 += (std::complex<double>(0,1)*Conj(Mu)*KroneckerDelta(0,gO2)) *
      tmp_1859;
   std::complex<double> tmp_1863;
   std::complex<double> tmp_1864;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1865;
      std::complex<double> tmp_1866;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1866 += ZV(gI1,3 + j1)*TYv(j1,j2);
      }
      tmp_1865 += tmp_1866;
      tmp_1864 += (Conj(ZE(gI2,j2))) * tmp_1865;
   }
   tmp_1863 += tmp_1864;
   tmp_1854 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1863;
   std::complex<double> tmp_1867;
   std::complex<double> tmp_1868;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1869;
      std::complex<double> tmp_1870;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1870 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1869 += tmp_1870;
      tmp_1868 += (ZV(gI1,j2)) * tmp_1869;
   }
   tmp_1867 += tmp_1868;
   tmp_1854 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1867;
   std::complex<double> tmp_1871;
   std::complex<double> tmp_1872;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1873;
      std::complex<double> tmp_1874;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1874 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1873 += tmp_1874;
      tmp_1872 += (ZV(gI1,j2)) * tmp_1873;
   }
   tmp_1871 += tmp_1872;
   tmp_1854 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1871;
   std::complex<double> tmp_1875;
   std::complex<double> tmp_1876;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1877;
      std::complex<double> tmp_1878;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1879;
         std::complex<double> tmp_1880;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1880 += Conj(Ye(j3,j1))*Yv(j2,j1);
         }
         tmp_1879 += tmp_1880;
         tmp_1878 += (ZV(gI1,3 + j2)) * tmp_1879;
      }
      tmp_1877 += tmp_1878;
      tmp_1876 += (Conj(ZE(gI2,3 + j3))) * tmp_1877;
   }
   tmp_1875 += tmp_1876;
   tmp_1854 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_1875;
   std::complex<double> tmp_1881;
   std::complex<double> tmp_1882;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1883;
      std::complex<double> tmp_1884;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1885;
         std::complex<double> tmp_1886;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1886 += Conj(Ye(j3,j1))*Yv(j2,j1);
         }
         tmp_1885 += tmp_1886;
         tmp_1884 += (ZV(gI1,3 + j2)) * tmp_1885;
      }
      tmp_1883 += tmp_1884;
      tmp_1882 += (Conj(ZE(gI2,3 + j3))) * tmp_1883;
   }
   tmp_1881 += tmp_1882;
   tmp_1854 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_1881;
   std::complex<double> tmp_1887;
   std::complex<double> tmp_1888;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1889;
      std::complex<double> tmp_1890;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1891;
         std::complex<double> tmp_1892;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1892 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1891 += tmp_1892;
         tmp_1890 += (Conj(ZE(gI2,j2))) * tmp_1891;
      }
      tmp_1889 += tmp_1890;
      tmp_1888 += (ZV(gI1,j3)) * tmp_1889;
   }
   tmp_1887 += tmp_1888;
   tmp_1854 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1887;
   std::complex<double> tmp_1893;
   std::complex<double> tmp_1894;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1895;
      std::complex<double> tmp_1896;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1897;
         std::complex<double> tmp_1898;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1898 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_1897 += tmp_1898;
         tmp_1896 += (Conj(ZE(gI2,j2))) * tmp_1897;
      }
      tmp_1895 += tmp_1896;
      tmp_1894 += (ZV(gI1,j3)) * tmp_1895;
   }
   tmp_1893 += tmp_1894;
   tmp_1854 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_1893;
   result += (std::complex<double>(0,-1)) * tmp_1854;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1899;
   std::complex<double> tmp_1900;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1901;
      std::complex<double> tmp_1902;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1903;
         std::complex<double> tmp_1904;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1904 += Conj(Mv(j1,j3))*Yv(j1,j2);
         }
         tmp_1903 += tmp_1904;
         tmp_1902 += (Conj(ZE(gI2,j2))) * tmp_1903;
      }
      tmp_1901 += tmp_1902;
      tmp_1900 += (Conj(ZV(gI1,3 + j3))) * tmp_1901;
   }
   tmp_1899 += tmp_1900;
   std::complex<double> tmp_1905;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1906;
      std::complex<double> tmp_1907;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1908;
         std::complex<double> tmp_1909;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1909 += Conj(Mv(j3,j1))*Yv(j1,j2);
         }
         tmp_1908 += tmp_1909;
         tmp_1907 += (Conj(ZE(gI2,j2))) * tmp_1908;
      }
      tmp_1906 += tmp_1907;
      tmp_1905 += (Conj(ZV(gI1,3 + j3))) * tmp_1906;
   }
   tmp_1899 += tmp_1905;
   result += (0.5*KroneckerDelta(1,gO2)) * tmp_1899;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZA(gI2,0) +
      KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) - KroneckerDelta(1,gO2)*ZH(
      gI2,1));

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

std::complex<double> CLASSNAME::CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const
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

double CLASSNAME::CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.1*KroneckerDelta(gI1,gI2)*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

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

   std::complex<double> tmp_1910;
   std::complex<double> tmp_1911;
   std::complex<double> tmp_1912;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1912 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1911 += tmp_1912;
   tmp_1910 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1911;
   std::complex<double> tmp_1913;
   std::complex<double> tmp_1914;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1914 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1913 += tmp_1914;
   tmp_1910 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1913;
   std::complex<double> tmp_1915;
   std::complex<double> tmp_1916;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1916 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1915 += tmp_1916;
   tmp_1910 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1915;
   std::complex<double> tmp_1917;
   std::complex<double> tmp_1918;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1918 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1917 += tmp_1918;
   tmp_1910 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1917;
   result += (std::complex<double>(0,-1)) * tmp_1910;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1919;
   std::complex<double> tmp_1920;
   std::complex<double> tmp_1921;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1921 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1920 += tmp_1921;
   tmp_1919 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1920;
   std::complex<double> tmp_1922;
   std::complex<double> tmp_1923;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1923 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1922 += tmp_1923;
   tmp_1919 += (std::complex<double>(0.,-0.7745966692414834)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1922;
   std::complex<double> tmp_1924;
   std::complex<double> tmp_1925;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1925 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1924 += tmp_1925;
   tmp_1919 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1924;
   std::complex<double> tmp_1926;
   std::complex<double> tmp_1927;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1927 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1926 += tmp_1927;
   tmp_1919 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1926;
   result += (std::complex<double>(0,-1)) * tmp_1919;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1928;
   std::complex<double> tmp_1929;
   std::complex<double> tmp_1930;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1930 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1929 += tmp_1930;
   tmp_1928 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1929;
   std::complex<double> tmp_1931;
   std::complex<double> tmp_1932;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1932 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1931 += tmp_1932;
   tmp_1928 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1931;
   std::complex<double> tmp_1933;
   std::complex<double> tmp_1934;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1934 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1933 += tmp_1934;
   tmp_1928 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1933;
   std::complex<double> tmp_1935;
   std::complex<double> tmp_1936;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1936 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1935 += tmp_1936;
   tmp_1928 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1935;
   result += (std::complex<double>(0,-1)) * tmp_1928;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1937;
   std::complex<double> tmp_1938;
   std::complex<double> tmp_1939;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1939 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1938 += tmp_1939;
   tmp_1937 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1938;
   std::complex<double> tmp_1940;
   std::complex<double> tmp_1941;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1941 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1940 += tmp_1941;
   tmp_1937 += (std::complex<double>(0.,0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1940;
   std::complex<double> tmp_1942;
   std::complex<double> tmp_1943;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1943 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1942 += tmp_1943;
   tmp_1937 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1942;
   result += (std::complex<double>(0,-1)) * tmp_1937;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1944;
   std::complex<double> tmp_1945;
   std::complex<double> tmp_1946;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1946 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1945 += tmp_1946;
   tmp_1944 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1945;
   std::complex<double> tmp_1947;
   std::complex<double> tmp_1948;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1948 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1947 += tmp_1948;
   tmp_1944 += (-3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1947;
   result += (0.16666666666666666) * tmp_1944;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1949;
   std::complex<double> tmp_1950;
   std::complex<double> tmp_1951;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1951 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1950 += tmp_1951;
   tmp_1949 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1950;
   std::complex<double> tmp_1952;
   std::complex<double> tmp_1953;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1953 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1952 += tmp_1953;
   tmp_1949 += (-(g2*Cos(ThetaW())) + 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1952;
   result += (0.5) * tmp_1949;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1954;
   std::complex<double> tmp_1955;
   std::complex<double> tmp_1956;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1956 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1955 += tmp_1956;
   tmp_1954 += (-3.0983866769659336*g1*Sin(ThetaW())) * tmp_1955;
   std::complex<double> tmp_1957;
   std::complex<double> tmp_1958;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1958 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1957 += tmp_1958;
   tmp_1954 += (3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1957;
   result += (0.16666666666666666) * tmp_1954;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1959;
   std::complex<double> tmp_1960;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1960 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1959 += tmp_1960;
   result += (0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))) *
      tmp_1959;

   return result;
}

std::complex<double> CLASSNAME::CpVZFvFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1961;
   std::complex<double> tmp_1962;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1962 += Conj(UV(gI2,j1))*UV(gI1,j1);
   }
   tmp_1961 += tmp_1962;
   result += (0.5*(-(g2*Cos(ThetaW())) - 0.7745966692414834*g1*Sin(ThetaW())))
      * tmp_1961;

   return result;
}

std::complex<double> CLASSNAME::CpVZFvFvPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1963;
   std::complex<double> tmp_1964;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1964 += Conj(UV(gI1,j1))*UV(gI2,j1);
   }
   tmp_1963 += tmp_1964;
   result += (0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))) *
      tmp_1963;

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

   result = 0.3872983346207417*g1*g2*Sin(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1)
      );

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

std::complex<double> CLASSNAME::CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const
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

std::complex<double> CLASSNAME::CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(ZA(gI2,0)*ZP(gI1,0) + ZA(gI2,1)*ZP
      (gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1965;
   std::complex<double> tmp_1966;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1966 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1965 += tmp_1966;
   result += (-0.7071067811865475*g2) * tmp_1965;

   return result;
}

double CLASSNAME::CpconjVWmbarFuFdPR(unsigned , unsigned ) const
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

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1967;
   std::complex<double> tmp_1968;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1968 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1967 += tmp_1968;
   result += (0.5*Sqr(g2)) * tmp_1967;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1969;
   std::complex<double> tmp_1970;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1970 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1969 += tmp_1970;
   result += (0.5*Sqr(g2)) * tmp_1969;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1971;
   std::complex<double> tmp_1972;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1972 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1971 += tmp_1972;
   result += (0.5*Sqr(g2)) * tmp_1971;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1973;
   std::complex<double> tmp_1974;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1974 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1973 += tmp_1974;
   result += (0.5*Sqr(g2)) * tmp_1973;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmFvFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1975;
   std::complex<double> tmp_1976;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1976 += Conj(ZEL(gI2,j1))*UV(gI1,j1);
   }
   tmp_1975 += tmp_1976;
   result += (-0.7071067811865475*g2) * tmp_1975;

   return result;
}

double CLASSNAME::CpconjVWmFvFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1977;
   std::complex<double> tmp_1978;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1978 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1977 += tmp_1978;
   result += (0.7071067811865475*g2) * tmp_1977;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1979;
   std::complex<double> tmp_1980;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1980 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1979 += tmp_1980;
   result += (0.7071067811865475*g2) * tmp_1979;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVPHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Sin(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1)
      );

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

std::complex<double> CLASSNAME::CpUChiconjHpmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(Conj(UM(gI2,1))*(5.477225575051661*g1*KroneckerDelta(0,gO2) +
      7.0710678118654755*g2*KroneckerDelta(1,gO2)) - 10*g2*Conj(UM(gI2,0))*
      KroneckerDelta(2,gO2))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(10*g2*KroneckerDelta(3,gO1)*UP(gI2,0) + 1.4142135623730951*(
      3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*KroneckerDelta(1,gO1))*UP(
      gI2,1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(Conj(ZN(gI2,2))*(3.872983346207417*g1*KroneckerDelta(0,gO2) -
      5*g2*KroneckerDelta(1,gO2))*ZH(gI1,0) - 5*g2*Conj(ZN(gI2,1))*KroneckerDelta(
      2,gO2)*ZH(gI1,0) - 3.872983346207417*g1*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2
      )*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,3))*KroneckerDelta(1,gO2)*ZH(gI1,1) + 5*g2*
      Conj(ZN(gI2,1))*KroneckerDelta(3,gO2)*ZH(gI1,1) + 3.872983346207417*g1*Conj(
      ZN(gI2,0))*(KroneckerDelta(2,gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1
      )));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(KroneckerDelta(2,gO1)*ZH(gI1,0)*(3.872983346207417*g1*ZN(gI2,0
      ) - 5*g2*ZN(gI2,1)) + KroneckerDelta(3,gO1)*ZH(gI1,1)*(-3.872983346207417*g1
      *ZN(gI2,0) + 5*g2*ZN(gI2,1)) + (3.872983346207417*g1*KroneckerDelta(0,gO1) -
      5*g2*KroneckerDelta(1,gO1))*(ZH(gI1,0)*ZN(gI2,2) - ZH(gI1,1)*ZN(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(Conj(ZN(gI1,2))*(3.872983346207417*g1
      *KroneckerDelta(0,gO2) - 5*g2*KroneckerDelta(1,gO2))*ZA(gI2,0) - 5*g2*Conj(
      ZN(gI1,1))*KroneckerDelta(2,gO2)*ZA(gI2,0) - 3.872983346207417*g1*Conj(ZN(
      gI1,3))*KroneckerDelta(0,gO2)*ZA(gI2,1) + 5*g2*Conj(ZN(gI1,3))*
      KroneckerDelta(1,gO2)*ZA(gI2,1) + 5*g2*Conj(ZN(gI1,1))*KroneckerDelta(3,gO2)
      *ZA(gI2,1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(KroneckerDelta(2,gO2)*ZA(
      gI2,0) - KroneckerDelta(3,gO2)*ZA(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(KroneckerDelta(2,gO1)*ZA(gI2,0)*(
      3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1)) + KroneckerDelta(3,gO1)*ZA(
      gI2,1)*(-3.872983346207417*g1*ZN(gI1,0) + 5*g2*ZN(gI1,1)) + (
      3.872983346207417*g1*KroneckerDelta(0,gO1) - 5*g2*KroneckerDelta(1,gO1))*(ZA
      (gI2,0)*ZN(gI1,2) - ZA(gI2,1)*ZN(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1981;
   std::complex<double> tmp_1982;
   std::complex<double> tmp_1983;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1983 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1982 += tmp_1983;
   tmp_1981 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO2)) * tmp_1982;
   std::complex<double> tmp_1984;
   std::complex<double> tmp_1985;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1985 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1984 += tmp_1985;
   tmp_1981 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO2)) * tmp_1984;
   std::complex<double> tmp_1986;
   std::complex<double> tmp_1987;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1988;
      std::complex<double> tmp_1989;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1989 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1988 += tmp_1989;
      tmp_1987 += (Conj(ZDL(gI2,j2))) * tmp_1988;
   }
   tmp_1986 += tmp_1987;
   tmp_1981 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1986;
   result += (std::complex<double>(0,-1)) * tmp_1981;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1990;
   std::complex<double> tmp_1991;
   std::complex<double> tmp_1992;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1992 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1991 += tmp_1992;
   tmp_1990 += (std::complex<double>(0.,-0.3651483716701107)*g1*KroneckerDelta(
      0,gO1)) * tmp_1991;
   std::complex<double> tmp_1993;
   std::complex<double> tmp_1994;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1995;
      std::complex<double> tmp_1996;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1996 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1995 += tmp_1996;
      tmp_1994 += (ZD(gI1,j2)) * tmp_1995;
   }
   tmp_1993 += tmp_1994;
   tmp_1990 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1993;
   result += (std::complex<double>(0,-1)) * tmp_1990;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1997;
   std::complex<double> tmp_1998;
   std::complex<double> tmp_1999;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1999 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1998 += tmp_1999;
   tmp_1997 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(0
      ,gO2)) * tmp_1998;
   std::complex<double> tmp_2000;
   std::complex<double> tmp_2001;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2001 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2000 += tmp_2001;
   tmp_1997 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO2)) * tmp_2000;
   std::complex<double> tmp_2002;
   std::complex<double> tmp_2003;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2004;
      std::complex<double> tmp_2005;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2005 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2004 += tmp_2005;
      tmp_2003 += (Conj(ZEL(gI2,j2))) * tmp_2004;
   }
   tmp_2002 += tmp_2003;
   tmp_1997 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_2002;
   result += (std::complex<double>(0,-1)) * tmp_1997;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2006;
   std::complex<double> tmp_2007;
   std::complex<double> tmp_2008;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2008 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_2007 += tmp_2008;
   tmp_2006 += (std::complex<double>(0.,-1.0954451150103321)*g1*KroneckerDelta(
      0,gO1)) * tmp_2007;
   std::complex<double> tmp_2009;
   std::complex<double> tmp_2010;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2011;
      std::complex<double> tmp_2012;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2012 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2011 += tmp_2012;
      tmp_2010 += (ZE(gI1,j2)) * tmp_2011;
   }
   tmp_2009 += tmp_2010;
   tmp_2006 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_2009;
   result += (std::complex<double>(0,-1)) * tmp_2006;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2013;
   std::complex<double> tmp_2014;
   std::complex<double> tmp_2015;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2015 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2014 += tmp_2015;
   tmp_2013 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO2)) * tmp_2014;
   std::complex<double> tmp_2016;
   std::complex<double> tmp_2017;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2017 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2016 += tmp_2017;
   tmp_2013 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      1,gO2)) * tmp_2016;
   std::complex<double> tmp_2018;
   std::complex<double> tmp_2019;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2020;
      std::complex<double> tmp_2021;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2021 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2020 += tmp_2021;
      tmp_2019 += (Conj(ZUL(gI2,j2))) * tmp_2020;
   }
   tmp_2018 += tmp_2019;
   tmp_2013 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_2018;
   result += (std::complex<double>(0,-1)) * tmp_2013;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2022;
   std::complex<double> tmp_2023;
   std::complex<double> tmp_2024;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2024 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2023 += tmp_2024;
   tmp_2022 += (std::complex<double>(0.,0.7302967433402214)*g1*KroneckerDelta(0
      ,gO1)) * tmp_2023;
   std::complex<double> tmp_2025;
   std::complex<double> tmp_2026;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2027;
      std::complex<double> tmp_2028;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2028 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2027 += tmp_2028;
      tmp_2026 += (ZU(gI1,j2)) * tmp_2027;
   }
   tmp_2025 += tmp_2026;
   tmp_2022 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_2025;
   result += (std::complex<double>(0,-1)) * tmp_2022;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2029;
   std::complex<double> tmp_2030;
   std::complex<double> tmp_2031;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2031 += Conj(UV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2030 += tmp_2031;
   tmp_2029 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(0
      ,gO2)) * tmp_2030;
   std::complex<double> tmp_2032;
   std::complex<double> tmp_2033;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2033 += Conj(UV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2032 += tmp_2033;
   tmp_2029 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      1,gO2)) * tmp_2032;
   std::complex<double> tmp_2034;
   std::complex<double> tmp_2035;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2036;
      std::complex<double> tmp_2037;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2037 += Yv(j1,j2)*ZV(gI1,3 + j1);
      }
      tmp_2036 += tmp_2037;
      tmp_2035 += (Conj(UV(gI2,j2))) * tmp_2036;
   }
   tmp_2034 += tmp_2035;
   tmp_2029 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_2034;
   result += (std::complex<double>(0,-1)) * tmp_2029;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2038;
   std::complex<double> tmp_2039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2040;
      std::complex<double> tmp_2041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2041 += Conj(Yv(j1,j2))*UV(gI2,3 + j1);
      }
      tmp_2040 += tmp_2041;
      tmp_2039 += (ZV(gI1,j2)) * tmp_2040;
   }
   tmp_2038 += tmp_2039;
   result += (-KroneckerDelta(3,gO1)) * tmp_2038;

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

std::complex<double> CLASSNAME::CpUFvconjHpmFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2042;
   std::complex<double> tmp_2043;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2044;
      std::complex<double> tmp_2045;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2045 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_2044 += tmp_2045;
      tmp_2043 += (Conj(ZEL(gI2,j2))) * tmp_2044;
   }
   tmp_2042 += tmp_2043;
   result += (ZP(gI1,1)) * tmp_2042;

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjHpmFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2046;
      std::complex<double> tmp_2047;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2047 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_2046 += tmp_2047;
      result += (ZP(gI1,0)) * tmp_2046;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvhhFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2048;
   std::complex<double> tmp_2049;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2050;
      std::complex<double> tmp_2051;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2051 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_2050 += tmp_2051;
      tmp_2049 += (Conj(UV(gI2,j2))) * tmp_2050;
   }
   tmp_2048 += tmp_2049;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2048;
   if (gO2 < 3) {
      std::complex<double> tmp_2052;
      std::complex<double> tmp_2053;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2053 += Conj(UV(gI2,3 + j1))*Yv(j1,gO2);
      }
      tmp_2052 += tmp_2053;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2052;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvhhFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2056;
      std::complex<double> tmp_2057;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2057 += Conj(Yv(j1,j2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_2056 += tmp_2057;
      tmp_2055 += (UV(gI2,j2)) * tmp_2056;
   }
   tmp_2054 += tmp_2055;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2054;
   if (gO1 < 3) {
      std::complex<double> tmp_2058;
      std::complex<double> tmp_2059;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2059 += Conj(Yv(j1,gO1))*UV(gI2,3 + j1);
      }
      tmp_2058 += tmp_2059;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2058;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjSeChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2060;
      std::complex<double> tmp_2061;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2061 += Ye(j1,gO2)*ZE(gI1,3 + j1);
      }
      tmp_2060 += tmp_2061;
      result += (Conj(UM(gI2,1))) * tmp_2060;
   }
   if (gO2 < 3) {
      result += -(g2*Conj(UM(gI2,0))*ZE(gI1,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjSeChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2062;
   std::complex<double> tmp_2063;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2064;
      std::complex<double> tmp_2065;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2065 += Conj(Yv(j1,j2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_2064 += tmp_2065;
      tmp_2063 += (ZE(gI1,j2)) * tmp_2064;
   }
   tmp_2062 += tmp_2063;
   result += (UP(gI2,1)) * tmp_2062;

   return result;
}

std::complex<double> CLASSNAME::CpUFvFvAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2066;
   std::complex<double> tmp_2067;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2068;
      std::complex<double> tmp_2069;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2069 += KroneckerDelta(gO2,3 + j1)*Yv(j1,j2);
      }
      tmp_2068 += tmp_2069;
      tmp_2067 += (Conj(UV(gI1,j2))) * tmp_2068;
   }
   tmp_2066 += tmp_2067;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
      tmp_2066;
   if (gO2 < 3) {
      std::complex<double> tmp_2070;
      std::complex<double> tmp_2071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2071 += Conj(UV(gI1,3 + j1))*Yv(j1,gO2);
      }
      tmp_2070 += tmp_2071;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_2070;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvFvAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2074;
      std::complex<double> tmp_2075;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2075 += Conj(Yv(j1,j2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_2074 += tmp_2075;
      tmp_2073 += (UV(gI1,j2)) * tmp_2074;
   }
   tmp_2072 += tmp_2073;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_2072
      ;
   if (gO1 < 3) {
      std::complex<double> tmp_2076;
      std::complex<double> tmp_2077;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2077 += Conj(Yv(j1,gO1))*UV(gI1,3 + j1);
      }
      tmp_2076 += tmp_2077;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_2076;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjSvChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2078;
      std::complex<double> tmp_2079;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2079 += Yv(j1,gO2)*ZV(gI1,3 + j1);
      }
      tmp_2078 += tmp_2079;
      result += (-Conj(ZN(gI2,3))) * tmp_2078;
   }
   if (gO2 < 3) {
      result += 0.5477225575051661*g1*Conj(ZN(gI2,0))*ZV(gI1,gO2);
   }
   if (gO2 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*ZV(gI1,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjSvChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2080;
   std::complex<double> tmp_2081;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2082;
      std::complex<double> tmp_2083;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2083 += Conj(Yv(j1,j2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_2082 += tmp_2083;
      tmp_2081 += (ZV(gI1,j2)) * tmp_2082;
   }
   tmp_2080 += tmp_2081;
   result += (-ZN(gI2,3)) * tmp_2080;

   return result;
}

double CLASSNAME::CpUFvconjVWmFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUFvconjVWmFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvVZFvPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.5*g2*Cos(ThetaW())*UV(gI2,gO2);
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Sin(ThetaW())*UV(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUFvVZFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.5*g2*Conj(UV(gI2,gO1))*Cos(ThetaW());
   }
   if (gO1 < 3) {
      result += -0.3872983346207417*g1*Conj(UV(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gI1,1))*
      KroneckerDelta(0,gO2)*ZA(gI2,0) + Conj(UM(gI1,0))*KroneckerDelta(1,gO2)*ZA(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*g2*(KroneckerDelta(1,
      gO1)*UP(gI1,0)*ZA(gI2,0) + KroneckerDelta(0,gO1)*UP(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*(Conj(UM(gI2,1))*KroneckerDelta(0,gO2)*ZH(
      gI1,0) + Conj(UM(gI2,0))*KroneckerDelta(1,gO2)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*(KroneckerDelta(1,gO1)*UP(gI2,0)*ZH(gI1,0) +
      KroneckerDelta(0,gO1)*UP(gI2,1)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(10*g2*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2) +
      1.4142135623730951*(3.872983346207417*g1*Conj(ZN(gI2,0)) + 5*g2*Conj(ZN(gI2,
      1)))*KroneckerDelta(1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(KroneckerDelta(1,gO1)*(5.477225575051661*g1*ZN(gI2,0) +
      7.0710678118654755*g2*ZN(gI2,1)) - 10*g2*KroneckerDelta(0,gO1)*ZN(gI2,2))*ZP
      (gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2084;
   std::complex<double> tmp_2085;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2086;
      std::complex<double> tmp_2087;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2087 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_2086 += tmp_2087;
      tmp_2085 += (Conj(ZD(gI2,j2))) * tmp_2086;
   }
   tmp_2084 += tmp_2085;
   result += (KroneckerDelta(1,gO2)) * tmp_2084;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2088;
   std::complex<double> tmp_2089;
   std::complex<double> tmp_2090;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2090 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_2089 += tmp_2090;
   tmp_2088 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_2089
      ;
   std::complex<double> tmp_2091;
   std::complex<double> tmp_2092;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2093;
      std::complex<double> tmp_2094;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2094 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2093 += tmp_2094;
      tmp_2092 += (ZUL(gI1,j2)) * tmp_2093;
   }
   tmp_2091 += tmp_2092;
   tmp_2088 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_2091;
   result += (std::complex<double>(0,-1)) * tmp_2088;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2095;
   std::complex<double> tmp_2096;
   std::complex<double> tmp_2097;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2097 += Conj(ZDL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2096 += tmp_2097;
   tmp_2095 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_2096
      ;
   std::complex<double> tmp_2098;
   std::complex<double> tmp_2099;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2100;
      std::complex<double> tmp_2101;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2101 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2100 += tmp_2101;
      tmp_2099 += (Conj(ZDL(gI2,j2))) * tmp_2100;
   }
   tmp_2098 += tmp_2099;
   tmp_2095 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_2098;
   result += (std::complex<double>(0,-1)) * tmp_2095;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2102;
   std::complex<double> tmp_2103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2104;
      std::complex<double> tmp_2105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2105 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2104 += tmp_2105;
      tmp_2103 += (ZU(gI1,j2)) * tmp_2104;
   }
   tmp_2102 += tmp_2103;
   result += (KroneckerDelta(1,gO1)) * tmp_2102;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2106;
   std::complex<double> tmp_2107;
   std::complex<double> tmp_2108;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2108 += Conj(ZEL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2107 += tmp_2108;
   tmp_2106 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_2107
      ;
   std::complex<double> tmp_2109;
   std::complex<double> tmp_2110;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2111;
      std::complex<double> tmp_2112;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2112 += Yv(j1,j2)*ZV(gI1,3 + j1);
      }
      tmp_2111 += tmp_2112;
      tmp_2110 += (Conj(ZEL(gI2,j2))) * tmp_2111;
   }
   tmp_2109 += tmp_2110;
   tmp_2106 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_2109;
   result += (std::complex<double>(0,-1)) * tmp_2106;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2113;
   std::complex<double> tmp_2114;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2115;
      std::complex<double> tmp_2116;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2116 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2115 += tmp_2116;
      tmp_2114 += (ZV(gI1,j2)) * tmp_2115;
   }
   tmp_2113 += tmp_2114;
   result += (KroneckerDelta(1,gO1)) * tmp_2113;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaSeFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2117;
   std::complex<double> tmp_2118;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2119;
      std::complex<double> tmp_2120;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2120 += Conj(UV(gI2,3 + j1))*Yv(j1,j2);
      }
      tmp_2119 += tmp_2120;
      tmp_2118 += (Conj(ZE(gI1,j2))) * tmp_2119;
   }
   tmp_2117 += tmp_2118;
   result += (KroneckerDelta(1,gO2)) * tmp_2117;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaSeFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2121;
   std::complex<double> tmp_2122;
   std::complex<double> tmp_2123;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2123 += Conj(ZE(gI1,j1))*UV(gI2,j1);
   }
   tmp_2122 += tmp_2123;
   tmp_2121 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_2122
      ;
   std::complex<double> tmp_2124;
   std::complex<double> tmp_2125;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2126;
      std::complex<double> tmp_2127;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2127 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2126 += tmp_2127;
      tmp_2125 += (UV(gI2,j2)) * tmp_2126;
   }
   tmp_2124 += tmp_2125;
   tmp_2121 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_2124;
   result += (std::complex<double>(0,-1)) * tmp_2121;

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

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2128;
      std::complex<double> tmp_2129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2129 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_2128 += tmp_2129;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2128;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2130;
      std::complex<double> tmp_2131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2131 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_2130 += tmp_2131;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2130;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2132;
      std::complex<double> tmp_2133;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2133 += Conj(UV(gI2,j2))*Ye(gO2,j2);
      }
      tmp_2132 += tmp_2133;
      result += (ZP(gI1,0)) * tmp_2132;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHpmFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2134;
      std::complex<double> tmp_2135;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2135 += Conj(Yv(j1,gO1))*UV(gI2,3 + j1);
      }
      tmp_2134 += tmp_2135;
      result += (ZP(gI1,1)) * tmp_2134;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2136;
      std::complex<double> tmp_2137;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2137 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2136 += tmp_2137;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_2136;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2138;
      std::complex<double> tmp_2139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2139 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_2138 += tmp_2139;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_2138;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2140;
      std::complex<double> tmp_2141;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2141 += Conj(ZV(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2140 += tmp_2141;
      result += (Conj(UM(gI2,1))) * tmp_2140;
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
      std::complex<double> tmp_2142;
      std::complex<double> tmp_2143;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2143 += Conj(Yv(j1,gO1))*Conj(ZV(gI1,3 + j1));
      }
      tmp_2142 += tmp_2143;
      result += (UP(gI2,1)) * tmp_2142;
   }

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
      std::complex<double> tmp_2144;
      std::complex<double> tmp_2145;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2145 += Conj(ZE(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2144 += tmp_2145;
      result += (-Conj(ZN(gI2,2))) * tmp_2144;
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
      std::complex<double> tmp_2146;
      std::complex<double> tmp_2147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2147 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2146 += tmp_2147;
      result += (-ZN(gI2,2)) * tmp_2146;
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
      result += -0.7071067811865475*g2*Conj(UV(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2148;
      std::complex<double> tmp_2149;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2149 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_2148 += tmp_2149;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2148;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2150;
      std::complex<double> tmp_2151;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2151 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_2150 += tmp_2151;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2150;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2152;
      std::complex<double> tmp_2153;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2153 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_2152 += tmp_2153;
      result += (ZP(gI1,0)) * tmp_2152;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2154;
      std::complex<double> tmp_2155;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2155 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_2154 += tmp_2155;
      result += (ZP(gI1,1)) * tmp_2154;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2156;
      std::complex<double> tmp_2157;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2157 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2156 += tmp_2157;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_2156;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2158;
      std::complex<double> tmp_2159;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2159 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_2158 += tmp_2159;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_2158;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2160;
      std::complex<double> tmp_2161;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2161 += Conj(ZU(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2160 += tmp_2161;
      result += (Conj(UM(gI2,1))) * tmp_2160;
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
      std::complex<double> tmp_2162;
      std::complex<double> tmp_2163;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2163 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2162 += tmp_2163;
      result += (UP(gI2,1)) * tmp_2162;
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
      std::complex<double> tmp_2164;
      std::complex<double> tmp_2165;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2165 += Conj(ZD(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2164 += tmp_2165;
      result += (-Conj(ZN(gI2,2))) * tmp_2164;
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
      std::complex<double> tmp_2166;
      std::complex<double> tmp_2167;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2167 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2166 += tmp_2167;
      result += (-ZN(gI2,2)) * tmp_2166;
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

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2168;
      std::complex<double> tmp_2169;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2169 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2168 += tmp_2169;
      result += (ZP(gI1,1)) * tmp_2168;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2170;
      std::complex<double> tmp_2171;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2171 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_2170 += tmp_2171;
      result += (ZP(gI1,0)) * tmp_2170;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2172;
      std::complex<double> tmp_2173;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2173 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2172 += tmp_2173;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2172;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2174;
      std::complex<double> tmp_2175;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2175 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_2174 += tmp_2175;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2174;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2176;
      std::complex<double> tmp_2177;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2177 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2176 += tmp_2177;
      result += (Conj(UP(gI1,1))) * tmp_2176;
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
      std::complex<double> tmp_2178;
      std::complex<double> tmp_2179;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2179 += Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2178 += tmp_2179;
      result += (UM(gI1,1)) * tmp_2178;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2180;
      std::complex<double> tmp_2181;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2181 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_2180 += tmp_2181;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_2180;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2182;
      std::complex<double> tmp_2183;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2183 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_2182 += tmp_2183;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_2182;
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
      std::complex<double> tmp_2184;
      std::complex<double> tmp_2185;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2185 += Conj(ZU(gI1,j2))*Yu(gO2,j2);
      }
      tmp_2184 += tmp_2185;
      result += (-Conj(ZN(gI2,3))) * tmp_2184;
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
      std::complex<double> tmp_2186;
      std::complex<double> tmp_2187;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2187 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2186 += tmp_2187;
      result += (-ZN(gI2,3)) * tmp_2186;
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

   std::complex<double> tmp_2188;
   std::complex<double> tmp_2189;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2189 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2188 += tmp_2189;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_2188;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2190;
   std::complex<double> tmp_2191;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2191 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_2190 += tmp_2191;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2190;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2192;
   std::complex<double> tmp_2193;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2193 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2192 += tmp_2193;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_2192;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2194;
   std::complex<double> tmp_2195;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2195 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2194 += tmp_2195;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2194;

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

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2196;
   std::complex<double> tmp_2197;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2198;
      std::complex<double> tmp_2199;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2199 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2198 += tmp_2199;
      tmp_2197 += (Conj(ZEL(gI2,j2))) * tmp_2198;
   }
   tmp_2196 += tmp_2197;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2196;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2200;
   std::complex<double> tmp_2201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2202;
      std::complex<double> tmp_2203;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2203 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2202 += tmp_2203;
      tmp_2201 += (ZEL(gO1,j2)) * tmp_2202;
   }
   tmp_2200 += tmp_2201;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2200;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2204;
   std::complex<double> tmp_2205;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2206;
      std::complex<double> tmp_2207;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2207 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2206 += tmp_2207;
      tmp_2205 += (Conj(UV(gI2,j2))) * tmp_2206;
   }
   tmp_2204 += tmp_2205;
   result += (ZP(gI1,0)) * tmp_2204;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2208;
   std::complex<double> tmp_2209;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2210;
      std::complex<double> tmp_2211;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2211 += Conj(Yv(j1,j2))*UV(gI2,3 + j1);
      }
      tmp_2210 += tmp_2211;
      tmp_2209 += (ZEL(gO1,j2)) * tmp_2210;
   }
   tmp_2208 += tmp_2209;
   result += (ZP(gI1,1)) * tmp_2208;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2212;
   std::complex<double> tmp_2213;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2214;
      std::complex<double> tmp_2215;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2215 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2214 += tmp_2215;
      tmp_2213 += (Conj(ZEL(gI1,j2))) * tmp_2214;
   }
   tmp_2212 += tmp_2213;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_2212;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2216;
   std::complex<double> tmp_2217;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2218;
      std::complex<double> tmp_2219;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2219 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_2218 += tmp_2219;
      tmp_2217 += (ZEL(gO1,j2)) * tmp_2218;
   }
   tmp_2216 += tmp_2217;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_2216
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2220;
   std::complex<double> tmp_2221;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2222;
      std::complex<double> tmp_2223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2223 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2222 += tmp_2223;
      tmp_2221 += (Conj(ZV(gI1,j2))) * tmp_2222;
   }
   tmp_2220 += tmp_2221;
   result += (Conj(UM(gI2,1))) * tmp_2220;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2224;
   std::complex<double> tmp_2225;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2225 += Conj(ZV(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2224 += tmp_2225;
   result += (-(g2*UP(gI2,0))) * tmp_2224;
   std::complex<double> tmp_2226;
   std::complex<double> tmp_2227;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2228;
      std::complex<double> tmp_2229;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2229 += Conj(Yv(j1,j2))*Conj(ZV(gI1,3 + j1));
      }
      tmp_2228 += tmp_2229;
      tmp_2227 += (ZEL(gO1,j2)) * tmp_2228;
   }
   tmp_2226 += tmp_2227;
   result += (UP(gI2,1)) * tmp_2226;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2230;
   std::complex<double> tmp_2231;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2231 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_2230 += tmp_2231;
   result += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_2230;
   std::complex<double> tmp_2232;
   std::complex<double> tmp_2233;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2234;
      std::complex<double> tmp_2235;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2235 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2234 += tmp_2235;
      tmp_2233 += (Conj(ZE(gI1,j2))) * tmp_2234;
   }
   tmp_2232 += tmp_2233;
   result += (-Conj(ZN(gI2,2))) * tmp_2232;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2236;
   std::complex<double> tmp_2237;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2237 += Conj(ZE(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2236 += tmp_2237;
   result += (0.7071067811865475*(0.7745966692414834*g1*ZN(gI2,0) + g2*ZN(gI2,1
      ))) * tmp_2236;
   std::complex<double> tmp_2238;
   std::complex<double> tmp_2239;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2240;
      std::complex<double> tmp_2241;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2241 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2240 += tmp_2241;
      tmp_2239 += (ZEL(gO1,j2)) * tmp_2240;
   }
   tmp_2238 += tmp_2239;
   result += (-ZN(gI2,2)) * tmp_2238;

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

double CLASSNAME::CpbarFeVWmFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2242;
   std::complex<double> tmp_2243;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2243 += Conj(UV(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_2242 += tmp_2243;
   result += (-0.7071067811865475*g2) * tmp_2242;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2244;
   std::complex<double> tmp_2245;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2246;
      std::complex<double> tmp_2247;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2247 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2246 += tmp_2247;
      tmp_2245 += (Conj(ZDL(gI2,j2))) * tmp_2246;
   }
   tmp_2244 += tmp_2245;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2244;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2248;
   std::complex<double> tmp_2249;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2250;
      std::complex<double> tmp_2251;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2251 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2250 += tmp_2251;
      tmp_2249 += (ZDL(gO1,j2)) * tmp_2250;
   }
   tmp_2248 += tmp_2249;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2248;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2252;
   std::complex<double> tmp_2253;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2254;
      std::complex<double> tmp_2255;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2255 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2254 += tmp_2255;
      tmp_2253 += (Conj(ZUL(gI2,j2))) * tmp_2254;
   }
   tmp_2252 += tmp_2253;
   result += (ZP(gI1,0)) * tmp_2252;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2256;
   std::complex<double> tmp_2257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2258;
      std::complex<double> tmp_2259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2259 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2258 += tmp_2259;
      tmp_2257 += (ZDL(gO1,j2)) * tmp_2258;
   }
   tmp_2256 += tmp_2257;
   result += (ZP(gI1,1)) * tmp_2256;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2260;
   std::complex<double> tmp_2261;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2262;
      std::complex<double> tmp_2263;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2263 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2262 += tmp_2263;
      tmp_2261 += (Conj(ZDL(gI1,j2))) * tmp_2262;
   }
   tmp_2260 += tmp_2261;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_2260;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2264;
   std::complex<double> tmp_2265;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2266;
      std::complex<double> tmp_2267;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2267 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_2266 += tmp_2267;
      tmp_2265 += (ZDL(gO1,j2)) * tmp_2266;
   }
   tmp_2264 += tmp_2265;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_2264
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2268;
   std::complex<double> tmp_2269;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2270;
      std::complex<double> tmp_2271;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2271 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2270 += tmp_2271;
      tmp_2269 += (Conj(ZU(gI1,j2))) * tmp_2270;
   }
   tmp_2268 += tmp_2269;
   result += (Conj(UM(gI2,1))) * tmp_2268;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2272;
   std::complex<double> tmp_2273;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2273 += Conj(ZU(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2272 += tmp_2273;
   result += (-(g2*UP(gI2,0))) * tmp_2272;
   std::complex<double> tmp_2274;
   std::complex<double> tmp_2275;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2276;
      std::complex<double> tmp_2277;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2277 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2276 += tmp_2277;
      tmp_2275 += (ZDL(gO1,j2)) * tmp_2276;
   }
   tmp_2274 += tmp_2275;
   result += (UP(gI2,1)) * tmp_2274;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2278;
   std::complex<double> tmp_2279;
   std::complex<double> tmp_2280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2280 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2279 += tmp_2280;
   tmp_2278 += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_2279;
   std::complex<double> tmp_2281;
   std::complex<double> tmp_2282;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2283;
      std::complex<double> tmp_2284;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2284 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2283 += tmp_2284;
      tmp_2282 += (Conj(ZD(gI1,j2))) * tmp_2283;
   }
   tmp_2281 += tmp_2282;
   tmp_2278 += (-3*Conj(ZN(gI2,2))) * tmp_2281;
   result += (0.3333333333333333) * tmp_2278;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2285;
   std::complex<double> tmp_2286;
   std::complex<double> tmp_2287;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2287 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2286 += tmp_2287;
   tmp_2285 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) - 3*g2*ZN(
      gI2,1))) * tmp_2286;
   std::complex<double> tmp_2288;
   std::complex<double> tmp_2289;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2290;
      std::complex<double> tmp_2291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2291 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2290 += tmp_2291;
      tmp_2289 += (ZDL(gO1,j2)) * tmp_2290;
   }
   tmp_2288 += tmp_2289;
   tmp_2285 += (-6*ZN(gI2,2)) * tmp_2288;
   result += (0.16666666666666666) * tmp_2285;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2292;
   std::complex<double> tmp_2293;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2293 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2292 += tmp_2293;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2292;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2294;
   std::complex<double> tmp_2295;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2295 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2294 += tmp_2295;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2294;

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

   std::complex<double> tmp_2296;
   std::complex<double> tmp_2297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2297 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2296 += tmp_2297;
   result += (-0.7071067811865475*g2) * tmp_2296;

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

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2298;
   std::complex<double> tmp_2299;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2300;
      std::complex<double> tmp_2301;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2301 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2300 += tmp_2301;
      tmp_2299 += (Conj(ZDL(gI2,j2))) * tmp_2300;
   }
   tmp_2298 += tmp_2299;
   result += (ZP(gI1,1)) * tmp_2298;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2302;
   std::complex<double> tmp_2303;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2304;
      std::complex<double> tmp_2305;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2305 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2304 += tmp_2305;
      tmp_2303 += (ZUL(gO1,j2)) * tmp_2304;
   }
   tmp_2302 += tmp_2303;
   result += (ZP(gI1,0)) * tmp_2302;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2306;
   std::complex<double> tmp_2307;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2308;
      std::complex<double> tmp_2309;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2309 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2308 += tmp_2309;
      tmp_2307 += (Conj(ZUL(gI2,j2))) * tmp_2308;
   }
   tmp_2306 += tmp_2307;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2306;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2310;
   std::complex<double> tmp_2311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2312;
      std::complex<double> tmp_2313;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2313 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2312 += tmp_2313;
      tmp_2311 += (ZUL(gO1,j2)) * tmp_2312;
   }
   tmp_2310 += tmp_2311;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2310;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2314;
   std::complex<double> tmp_2315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2316;
      std::complex<double> tmp_2317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2317 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2316 += tmp_2317;
      tmp_2315 += (Conj(ZD(gI2,j2))) * tmp_2316;
   }
   tmp_2314 += tmp_2315;
   result += (Conj(UP(gI1,1))) * tmp_2314;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2318;
   std::complex<double> tmp_2319;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2319 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2318 += tmp_2319;
   result += (-(g2*UM(gI1,0))) * tmp_2318;
   std::complex<double> tmp_2320;
   std::complex<double> tmp_2321;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2322;
      std::complex<double> tmp_2323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2323 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2322 += tmp_2323;
      tmp_2321 += (ZUL(gO1,j2)) * tmp_2322;
   }
   tmp_2320 += tmp_2321;
   result += (UM(gI1,1)) * tmp_2320;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2324;
   std::complex<double> tmp_2325;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2326;
      std::complex<double> tmp_2327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2327 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2326 += tmp_2327;
      tmp_2325 += (Conj(ZUL(gI1,j2))) * tmp_2326;
   }
   tmp_2324 += tmp_2325;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
      tmp_2324;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2328;
   std::complex<double> tmp_2329;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2330;
      std::complex<double> tmp_2331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2331 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2330 += tmp_2331;
      tmp_2329 += (ZUL(gO1,j2)) * tmp_2330;
   }
   tmp_2328 += tmp_2329;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_2328
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2332;
   std::complex<double> tmp_2333;
   std::complex<double> tmp_2334;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2334 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2333 += tmp_2334;
   tmp_2332 += (2.1908902300206643*g1*Conj(ZN(gI2,0))) * tmp_2333;
   std::complex<double> tmp_2335;
   std::complex<double> tmp_2336;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2337;
      std::complex<double> tmp_2338;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2338 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2337 += tmp_2338;
      tmp_2336 += (Conj(ZU(gI1,j2))) * tmp_2337;
   }
   tmp_2335 += tmp_2336;
   tmp_2332 += (-3*Conj(ZN(gI2,3))) * tmp_2335;
   result += (0.3333333333333333) * tmp_2332;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2339;
   std::complex<double> tmp_2340;
   std::complex<double> tmp_2341;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2341 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2340 += tmp_2341;
   tmp_2339 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) + 3*g2*ZN(
      gI2,1))) * tmp_2340;
   std::complex<double> tmp_2342;
   std::complex<double> tmp_2343;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2344;
      std::complex<double> tmp_2345;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2345 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2344 += tmp_2345;
      tmp_2343 += (ZUL(gO1,j2)) * tmp_2344;
   }
   tmp_2342 += tmp_2343;
   tmp_2339 += (-6*ZN(gI2,3)) * tmp_2342;
   result += (0.16666666666666666) * tmp_2339;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2346;
   std::complex<double> tmp_2347;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2347 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2346 += tmp_2347;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2346;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2348;
   std::complex<double> tmp_2349;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2349 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2348 += tmp_2349;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2348;

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

   std::complex<double> tmp_2350;
   std::complex<double> tmp_2351;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2351 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2350 += tmp_2351;
   result += (-0.7071067811865475*g2) * tmp_2350;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2352;
   std::complex<double> tmp_2353;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2353 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2352 += tmp_2353;
   result += (-0.5) * tmp_2352;
   std::complex<double> tmp_2354;
   std::complex<double> tmp_2355;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2355 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2354 += tmp_2355;
   result += (-1) * tmp_2354;
   std::complex<double> tmp_2356;
   std::complex<double> tmp_2357;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2357 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2356 += tmp_2357;
   result += (-0.5) * tmp_2356;
   std::complex<double> tmp_2358;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2359;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2359 += (Conj(CpconjUSdFuChaPL(gO2,gI1,gI2))*
            CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPR(gO2,gI1,gI2))*
            CpconjUSdFuChaPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha(gI2));
      }
      tmp_2358 += tmp_2359;
   }
   result += tmp_2358;
   std::complex<double> tmp_2360;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2361;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2361 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2360 += tmp_2361;
   }
   result += tmp_2360;
   std::complex<double> tmp_2362;
   std::complex<double> tmp_2363;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2364;
      std::complex<double> tmp_2365;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2365 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2364 += tmp_2365;
      tmp_2363 += (MFd(gI1)) * tmp_2364;
   }
   tmp_2362 += tmp_2363;
   result += (-2) * tmp_2362;
   std::complex<double> tmp_2366;
   std::complex<double> tmp_2367;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2368;
      std::complex<double> tmp_2369;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2369 += B0(p,MFu(gI1),MCha(gI2))*(Conj(CpconjUSdFuChaPR(gO2,
            gI1,gI2))*CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPL(gO2,
            gI1,gI2))*CpconjUSdFuChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2368 += tmp_2369;
      tmp_2367 += (MFu(gI1)) * tmp_2368;
   }
   tmp_2366 += tmp_2367;
   result += (-2) * tmp_2366;
   std::complex<double> tmp_2370;
   std::complex<double> tmp_2371;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2371 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2370 += tmp_2371;
   result += (-1) * tmp_2370;
   std::complex<double> tmp_2372;
   std::complex<double> tmp_2373;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2373 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2372 += tmp_2373;
   result += (-1) * tmp_2372;
   std::complex<double> tmp_2374;
   std::complex<double> tmp_2375;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2375 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2374 += tmp_2375;
   result += (-1) * tmp_2374;
   std::complex<double> tmp_2376;
   std::complex<double> tmp_2377;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2377 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2376 += tmp_2377;
   result += (-1) * tmp_2376;
   std::complex<double> tmp_2378;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2379;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2379 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2378 += tmp_2379;
   }
   result += tmp_2378;
   std::complex<double> tmp_2380;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2381;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2381 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2380 += tmp_2381;
   }
   result += tmp_2380;
   std::complex<double> tmp_2382;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2383;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2383 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2382 += tmp_2383;
   }
   result += tmp_2382;
   std::complex<double> tmp_2384;
   std::complex<double> tmp_2385;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2385 += (Conj(CpconjUSdGluFdPL(gO2,1,gI2))*CpconjUSdGluFdPL(gO1,1,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,1,gI2))*CpconjUSdGluFdPR(gO1,1,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2384 += tmp_2385;
   result += (1.3333333333333333) * tmp_2384;
   std::complex<double> tmp_2386;
   std::complex<double> tmp_2387;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2387 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2386 += tmp_2387;
   result += (1.3333333333333333) * tmp_2386;
   std::complex<double> tmp_2388;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2388 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2388;
   std::complex<double> tmp_2389;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2389 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2389;
   std::complex<double> tmp_2390;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2390 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2390;
   std::complex<double> tmp_2391;
   std::complex<double> tmp_2392;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2392 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,1,gI2))*
         CpconjUSdGluFdPL(gO1,1,gI2) + Conj(CpconjUSdGluFdPL(gO2,1,gI2))*
         CpconjUSdGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_2391 += tmp_2392;
   result += (-2.6666666666666665*MGlu) * tmp_2391;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2393;
   std::complex<double> tmp_2394;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2394 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2393 += tmp_2394;
   result += (-0.5) * tmp_2393;
   std::complex<double> tmp_2395;
   std::complex<double> tmp_2396;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2396 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2395 += tmp_2396;
   result += (-1) * tmp_2395;
   std::complex<double> tmp_2397;
   std::complex<double> tmp_2398;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2398 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2397 += tmp_2398;
   result += (-0.5) * tmp_2397;
   std::complex<double> tmp_2399;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2400;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2400 += (Conj(CpconjUSubarChaFdPL(gO2,gI1,gI2))*
            CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(CpconjUSubarChaFdPR(gO2,gI1,
            gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFd(gI2));
      }
      tmp_2399 += tmp_2400;
   }
   result += tmp_2399;
   std::complex<double> tmp_2401;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2402;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2402 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2401 += tmp_2402;
   }
   result += tmp_2401;
   std::complex<double> tmp_2403;
   std::complex<double> tmp_2404;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2405;
      std::complex<double> tmp_2406;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2406 += B0(p,MCha(gI1),MFd(gI2))*(Conj(CpconjUSubarChaFdPR(
            gO2,gI1,gI2))*CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChaFdPL(gO2,gI1,gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2405 += tmp_2406;
      tmp_2404 += (MCha(gI1)) * tmp_2405;
   }
   tmp_2403 += tmp_2404;
   result += (-2) * tmp_2403;
   std::complex<double> tmp_2407;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2408;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2408 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2407 += tmp_2408;
   }
   result += tmp_2407;
   std::complex<double> tmp_2409;
   std::complex<double> tmp_2410;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2411;
      std::complex<double> tmp_2412;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2412 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2411 += tmp_2412;
      tmp_2410 += (MFu(gI1)) * tmp_2411;
   }
   tmp_2409 += tmp_2410;
   result += (-2) * tmp_2409;
   std::complex<double> tmp_2413;
   std::complex<double> tmp_2414;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2414 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2413 += tmp_2414;
   result += (-1) * tmp_2413;
   std::complex<double> tmp_2415;
   std::complex<double> tmp_2416;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2416 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2415 += tmp_2416;
   result += (-1) * tmp_2415;
   std::complex<double> tmp_2417;
   std::complex<double> tmp_2418;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2418 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2417 += tmp_2418;
   result += (-1) * tmp_2417;
   std::complex<double> tmp_2419;
   std::complex<double> tmp_2420;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2420 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2419 += tmp_2420;
   result += (-1) * tmp_2419;
   std::complex<double> tmp_2421;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2422;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2422 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2421 += tmp_2422;
   }
   result += tmp_2421;
   std::complex<double> tmp_2423;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2424;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2424 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2423 += tmp_2424;
   }
   result += tmp_2423;
   std::complex<double> tmp_2425;
   std::complex<double> tmp_2426;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2426 += (Conj(CpconjUSuGluFuPL(gO2,1,gI2))*CpconjUSuGluFuPL(gO1,1,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,1,gI2))*CpconjUSuGluFuPR(gO1,1,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2425 += tmp_2426;
   result += (1.3333333333333333) * tmp_2425;
   std::complex<double> tmp_2427;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2427 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2427;
   std::complex<double> tmp_2428;
   std::complex<double> tmp_2429;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2429 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2428 += tmp_2429;
   result += (1.3333333333333333) * tmp_2428;
   std::complex<double> tmp_2430;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2430 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2430;
   std::complex<double> tmp_2431;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2431 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2431;
   std::complex<double> tmp_2432;
   std::complex<double> tmp_2433;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2433 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,1,gI2))*
         CpconjUSuGluFuPL(gO1,1,gI2) + Conj(CpconjUSuGluFuPL(gO2,1,gI2))*
         CpconjUSuGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_2432 += tmp_2433;
   result += (-2.6666666666666665*MGlu) * tmp_2432;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2434;
   std::complex<double> tmp_2435;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2435 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2434 += tmp_2435;
   result += (-0.5) * tmp_2434;
   std::complex<double> tmp_2436;
   std::complex<double> tmp_2437;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2437 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2436 += tmp_2437;
   result += (-1) * tmp_2436;
   std::complex<double> tmp_2438;
   std::complex<double> tmp_2439;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2439 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2438 += tmp_2439;
   result += (-0.5) * tmp_2438;
   std::complex<double> tmp_2440;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2441;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2441 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_2440 += tmp_2441;
   }
   result += tmp_2440;
   std::complex<double> tmp_2442;
   std::complex<double> tmp_2443;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2444;
      std::complex<double> tmp_2445;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2445 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2444 += tmp_2445;
      tmp_2443 += (MFe(gI1)) * tmp_2444;
   }
   tmp_2442 += tmp_2443;
   result += (-2) * tmp_2442;
   std::complex<double> tmp_2446;
   std::complex<double> tmp_2447;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2447 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2446 += tmp_2447;
   result += (-3) * tmp_2446;
   std::complex<double> tmp_2448;
   std::complex<double> tmp_2449;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2449 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2448 += tmp_2449;
   result += (-1) * tmp_2448;
   std::complex<double> tmp_2450;
   std::complex<double> tmp_2451;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2451 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2450 += tmp_2451;
   result += (-3) * tmp_2450;
   std::complex<double> tmp_2452;
   std::complex<double> tmp_2453;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2453 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2452 += tmp_2453;
   result += (-1) * tmp_2452;
   std::complex<double> tmp_2454;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2455;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2455 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeconjSvHpm(gO2
            ,gI1,gI2))*CpconjUSeconjSvHpm(gO1,gI1,gI2);
      }
      tmp_2454 += tmp_2455;
   }
   result += tmp_2454;
   std::complex<double> tmp_2456;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2457;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2457 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2456 += tmp_2457;
   }
   result += tmp_2456;
   std::complex<double> tmp_2458;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2459;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2459 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2458 += tmp_2459;
   }
   result += tmp_2458;
   std::complex<double> tmp_2460;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2461;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2461 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_2460 += tmp_2461;
   }
   result += tmp_2460;
   std::complex<double> tmp_2462;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2463;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2463 += (Conj(CpconjUSeFvChaPL(gO2,gI1,gI2))*
            CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPR(gO2,gI1,gI2))*
            CpconjUSeFvChaPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha(gI2));
      }
      tmp_2462 += tmp_2463;
   }
   result += tmp_2462;
   std::complex<double> tmp_2464;
   std::complex<double> tmp_2465;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2466;
      std::complex<double> tmp_2467;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2467 += B0(p,MFv(gI1),MCha(gI2))*(Conj(CpconjUSeFvChaPR(gO2,
            gI1,gI2))*CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPL(gO2,
            gI1,gI2))*CpconjUSeFvChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2466 += tmp_2467;
      tmp_2465 += (MFv(gI1)) * tmp_2466;
   }
   tmp_2464 += tmp_2465;
   result += (-2) * tmp_2464;
   std::complex<double> tmp_2468;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2468 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2468;
   std::complex<double> tmp_2469;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2469 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2469;
   std::complex<double> tmp_2470;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2470 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2470;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2471;
   std::complex<double> tmp_2472;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2472 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2471 += tmp_2472;
   result += (-0.5) * tmp_2471;
   std::complex<double> tmp_2473;
   std::complex<double> tmp_2474;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2474 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2473 += tmp_2474;
   result += (-1) * tmp_2473;
   std::complex<double> tmp_2475;
   std::complex<double> tmp_2476;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2476 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2475 += tmp_2476;
   result += (-0.5) * tmp_2475;
   std::complex<double> tmp_2477;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2478;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2478 += (Conj(CpconjUSvbarChaFePL(gO2,gI1,gI2))*
            CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(CpconjUSvbarChaFePR(gO2,gI1,
            gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFe(gI2));
      }
      tmp_2477 += tmp_2478;
   }
   result += tmp_2477;
   std::complex<double> tmp_2479;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2480;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2480 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2479 += tmp_2480;
   }
   result += tmp_2479;
   std::complex<double> tmp_2481;
   std::complex<double> tmp_2482;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2483;
      std::complex<double> tmp_2484;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2484 += B0(p,MCha(gI1),MFe(gI2))*(Conj(CpconjUSvbarChaFePR(
            gO2,gI1,gI2))*CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(
            CpconjUSvbarChaFePL(gO2,gI1,gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2483 += tmp_2484;
      tmp_2482 += (MCha(gI1)) * tmp_2483;
   }
   tmp_2481 += tmp_2482;
   result += (-2) * tmp_2481;
   std::complex<double> tmp_2485;
   std::complex<double> tmp_2486;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2486 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2485 += tmp_2486;
   result += (-3) * tmp_2485;
   std::complex<double> tmp_2487;
   std::complex<double> tmp_2488;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2488 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2487 += tmp_2488;
   result += (-1) * tmp_2487;
   std::complex<double> tmp_2489;
   std::complex<double> tmp_2490;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2490 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2489 += tmp_2490;
   result += (-3) * tmp_2489;
   std::complex<double> tmp_2491;
   std::complex<double> tmp_2492;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2492 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2491 += tmp_2492;
   result += (-1) * tmp_2491;
   std::complex<double> tmp_2493;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2494;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2494 += B0(p,MSe(gI1),MHpm(gI2))*Conj(CpconjUSvconjSeHpm(gO2
            ,gI1,gI2))*CpconjUSvconjSeHpm(gO1,gI1,gI2);
      }
      tmp_2493 += tmp_2494;
   }
   result += tmp_2493;
   std::complex<double> tmp_2495;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2496;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2496 += B0(p,MSv(gI1),MAh(gI2))*Conj(CpconjUSvconjSvAh(gO2,
            gI1,gI2))*CpconjUSvconjSvAh(gO1,gI1,gI2);
      }
      tmp_2495 += tmp_2496;
   }
   result += tmp_2495;
   std::complex<double> tmp_2497;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2498;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2498 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvconjSvhh(gO2,
            gI1,gI2))*CpconjUSvconjSvhh(gO1,gI1,gI2);
      }
      tmp_2497 += tmp_2498;
   }
   result += tmp_2497;
   std::complex<double> tmp_2499;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2500;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2500 += B0(p,MSv(gI1),MAh(gI2))*Conj(CpconjUSvSvAh(gO2,gI1,
            gI2))*CpconjUSvSvAh(gO1,gI1,gI2);
      }
      tmp_2499 += tmp_2500;
   }
   result += tmp_2499;
   std::complex<double> tmp_2501;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2502;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2502 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2501 += tmp_2502;
   }
   result += tmp_2501;
   std::complex<double> tmp_2503;
   std::complex<double> tmp_2504;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2505;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2505 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2504 += tmp_2505;
   }
   tmp_2503 += tmp_2504;
   result += (0.5) * tmp_2503;
   std::complex<double> tmp_2506;
   std::complex<double> tmp_2507;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2508;
      std::complex<double> tmp_2509;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2509 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2508 += tmp_2509;
      tmp_2507 += (MFv(gI1)) * tmp_2508;
   }
   tmp_2506 += tmp_2507;
   result += (-1) * tmp_2506;
   std::complex<double> tmp_2510;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2510 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2510;
   std::complex<double> tmp_2511;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2511 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2511;

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
   std::complex<double> tmp_2512;
   std::complex<double> tmp_2513;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2513 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2512 += tmp_2513;
   result += (-0.5) * tmp_2512;
   std::complex<double> tmp_2514;
   std::complex<double> tmp_2515;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2515 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2514 += tmp_2515;
   result += (-1) * tmp_2514;
   std::complex<double> tmp_2516;
   std::complex<double> tmp_2517;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2517 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2516 += tmp_2517;
   result += (-0.5) * tmp_2516;
   std::complex<double> tmp_2518;
   std::complex<double> tmp_2519;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2520;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2520 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2519 += tmp_2520;
   }
   tmp_2518 += tmp_2519;
   result += (0.5) * tmp_2518;
   std::complex<double> tmp_2521;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2522;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2522 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2521 += tmp_2522;
   }
   result += tmp_2521;
   std::complex<double> tmp_2523;
   std::complex<double> tmp_2524;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2525;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2525 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2524 += tmp_2525;
   }
   tmp_2523 += tmp_2524;
   result += (0.5) * tmp_2523;
   std::complex<double> tmp_2526;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2527;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2527 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*
            CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*
            CpUhhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2526 += tmp_2527;
   }
   result += tmp_2526;
   std::complex<double> tmp_2528;
   std::complex<double> tmp_2529;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2530;
      std::complex<double> tmp_2531;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2531 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2
            ,gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2530 += tmp_2531;
      tmp_2529 += (MCha(gI1)) * tmp_2530;
   }
   tmp_2528 += tmp_2529;
   result += (-2) * tmp_2528;
   std::complex<double> tmp_2532;
   std::complex<double> tmp_2533;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2534;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2534 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2533 += tmp_2534;
   }
   tmp_2532 += tmp_2533;
   result += (3) * tmp_2532;
   std::complex<double> tmp_2535;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2536;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2536 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2535 += tmp_2536;
   }
   result += tmp_2535;
   std::complex<double> tmp_2537;
   std::complex<double> tmp_2538;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2539;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2539 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2538 += tmp_2539;
   }
   tmp_2537 += tmp_2538;
   result += (3) * tmp_2537;
   std::complex<double> tmp_2540;
   std::complex<double> tmp_2541;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2542;
      std::complex<double> tmp_2543;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2543 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2542 += tmp_2543;
      tmp_2541 += (MFd(gI1)) * tmp_2542;
   }
   tmp_2540 += tmp_2541;
   result += (-6) * tmp_2540;
   std::complex<double> tmp_2544;
   std::complex<double> tmp_2545;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2546;
      std::complex<double> tmp_2547;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2547 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2546 += tmp_2547;
      tmp_2545 += (MFe(gI1)) * tmp_2546;
   }
   tmp_2544 += tmp_2545;
   result += (-2) * tmp_2544;
   std::complex<double> tmp_2548;
   std::complex<double> tmp_2549;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2550;
      std::complex<double> tmp_2551;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2551 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2550 += tmp_2551;
      tmp_2549 += (MFu(gI1)) * tmp_2550;
   }
   tmp_2548 += tmp_2549;
   result += (-6) * tmp_2548;
   std::complex<double> tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2554;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2554 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2553 += tmp_2554;
   }
   tmp_2552 += tmp_2553;
   result += (0.5) * tmp_2552;
   std::complex<double> tmp_2555;
   std::complex<double> tmp_2556;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2557;
      std::complex<double> tmp_2558;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2558 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,
            gI1,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2)
            )*CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2557 += tmp_2558;
      tmp_2556 += (MChi(gI1)) * tmp_2557;
   }
   tmp_2555 += tmp_2556;
   result += (-1) * tmp_2555;
   std::complex<double> tmp_2559;
   std::complex<double> tmp_2560;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2560 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2559 += tmp_2560;
   result += (-3) * tmp_2559;
   std::complex<double> tmp_2561;
   std::complex<double> tmp_2562;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2562 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2561 += tmp_2562;
   result += (-1) * tmp_2561;
   std::complex<double> tmp_2563;
   std::complex<double> tmp_2564;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2564 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2563 += tmp_2564;
   result += (-3) * tmp_2563;
   std::complex<double> tmp_2565;
   std::complex<double> tmp_2566;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2566 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2565 += tmp_2566;
   result += (-1) * tmp_2565;
   std::complex<double> tmp_2567;
   std::complex<double> tmp_2568;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2569;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2569 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2568 += tmp_2569;
   }
   tmp_2567 += tmp_2568;
   result += (3) * tmp_2567;
   std::complex<double> tmp_2570;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2571;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2571 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2570 += tmp_2571;
   }
   result += tmp_2570;
   std::complex<double> tmp_2572;
   std::complex<double> tmp_2573;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2574;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2574 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2573 += tmp_2574;
   }
   tmp_2572 += tmp_2573;
   result += (3) * tmp_2572;
   std::complex<double> tmp_2575;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2576;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2576 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2575 += tmp_2576;
   }
   result += tmp_2575;
   std::complex<double> tmp_2577;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2578;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2578 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhSvSv(gO2,gI1,gI2))
            *CpUhhSvSv(gO1,gI1,gI2);
      }
      tmp_2577 += tmp_2578;
   }
   result += tmp_2577;
   std::complex<double> tmp_2579;
   std::complex<double> tmp_2580;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2581;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2581 += (Conj(CpUhhFvFvPL(gO2,gI1,gI2))*CpUhhFvFvPL(gO1,gI1,
            gI2) + Conj(CpUhhFvFvPR(gO2,gI1,gI2))*CpUhhFvFvPR(gO1,gI1,gI2))*G0(p,
            MFv(gI1),MFv(gI2));
      }
      tmp_2580 += tmp_2581;
   }
   tmp_2579 += tmp_2580;
   result += (0.5) * tmp_2579;
   std::complex<double> tmp_2582;
   std::complex<double> tmp_2583;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2584;
      std::complex<double> tmp_2585;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2585 += B0(p,MFv(gI1),MFv(gI2))*(Conj(CpUhhFvFvPR(gO2,gI1,
            gI2))*CpUhhFvFvPL(gO1,gI1,gI2) + Conj(CpUhhFvFvPL(gO2,gI1,gI2))*
            CpUhhFvFvPR(gO1,gI1,gI2))*MFv(gI2);
      }
      tmp_2584 += tmp_2585;
      tmp_2583 += (MFv(gI1)) * tmp_2584;
   }
   tmp_2582 += tmp_2583;
   result += (-1) * tmp_2582;
   std::complex<double> tmp_2586;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2586 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2586;
   std::complex<double> tmp_2587;
   std::complex<double> tmp_2588;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2588 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2587 += tmp_2588;
   result += (2) * tmp_2587;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2589;
   std::complex<double> tmp_2590;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2590 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2589 += tmp_2590;
   result += (-0.5) * tmp_2589;
   std::complex<double> tmp_2591;
   std::complex<double> tmp_2592;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2592 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2591 += tmp_2592;
   result += (-1) * tmp_2591;
   std::complex<double> tmp_2593;
   std::complex<double> tmp_2594;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2594 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2593 += tmp_2594;
   result += (-0.5) * tmp_2593;
   std::complex<double> tmp_2595;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2596;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2596 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2595 += tmp_2596;
   }
   result += tmp_2595;
   std::complex<double> tmp_2597;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2598;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2598 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2597 += tmp_2598;
   }
   result += tmp_2597;
   std::complex<double> tmp_2599;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2600;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2600 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*
            CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*
            CpUAhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2599 += tmp_2600;
   }
   result += tmp_2599;
   std::complex<double> tmp_2601;
   std::complex<double> tmp_2602;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2603;
      std::complex<double> tmp_2604;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2604 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2
            ,gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2603 += tmp_2604;
      tmp_2602 += (MCha(gI1)) * tmp_2603;
   }
   tmp_2601 += tmp_2602;
   result += (-2) * tmp_2601;
   std::complex<double> tmp_2605;
   std::complex<double> tmp_2606;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2607;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2607 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2606 += tmp_2607;
   }
   tmp_2605 += tmp_2606;
   result += (3) * tmp_2605;
   std::complex<double> tmp_2608;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2609;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2609 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2608 += tmp_2609;
   }
   result += tmp_2608;
   std::complex<double> tmp_2610;
   std::complex<double> tmp_2611;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2612;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2612 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2611 += tmp_2612;
   }
   tmp_2610 += tmp_2611;
   result += (3) * tmp_2610;
   std::complex<double> tmp_2613;
   std::complex<double> tmp_2614;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2615;
      std::complex<double> tmp_2616;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2616 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2615 += tmp_2616;
      tmp_2614 += (MFd(gI1)) * tmp_2615;
   }
   tmp_2613 += tmp_2614;
   result += (-6) * tmp_2613;
   std::complex<double> tmp_2617;
   std::complex<double> tmp_2618;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2619;
      std::complex<double> tmp_2620;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2620 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2619 += tmp_2620;
      tmp_2618 += (MFe(gI1)) * tmp_2619;
   }
   tmp_2617 += tmp_2618;
   result += (-2) * tmp_2617;
   std::complex<double> tmp_2621;
   std::complex<double> tmp_2622;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2623;
      std::complex<double> tmp_2624;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2624 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2623 += tmp_2624;
      tmp_2622 += (MFu(gI1)) * tmp_2623;
   }
   tmp_2621 += tmp_2622;
   result += (-6) * tmp_2621;
   std::complex<double> tmp_2625;
   std::complex<double> tmp_2626;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2627;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2627 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2626 += tmp_2627;
   }
   tmp_2625 += tmp_2626;
   result += (0.5) * tmp_2625;
   std::complex<double> tmp_2628;
   std::complex<double> tmp_2629;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2630;
      std::complex<double> tmp_2631;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2631 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,
            gI1,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2)
            )*CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2630 += tmp_2631;
      tmp_2629 += (MChi(gI1)) * tmp_2630;
   }
   tmp_2628 += tmp_2629;
   result += (-1) * tmp_2628;
   std::complex<double> tmp_2632;
   std::complex<double> tmp_2633;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2633 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2632 += tmp_2633;
   result += (-3) * tmp_2632;
   std::complex<double> tmp_2634;
   std::complex<double> tmp_2635;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2635 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2634 += tmp_2635;
   result += (-1) * tmp_2634;
   std::complex<double> tmp_2636;
   std::complex<double> tmp_2637;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2637 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2636 += tmp_2637;
   result += (-3) * tmp_2636;
   std::complex<double> tmp_2638;
   std::complex<double> tmp_2639;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2639 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2638 += tmp_2639;
   result += (-1) * tmp_2638;
   std::complex<double> tmp_2640;
   std::complex<double> tmp_2641;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2642;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2642 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2641 += tmp_2642;
   }
   tmp_2640 += tmp_2641;
   result += (3) * tmp_2640;
   std::complex<double> tmp_2643;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2644;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2644 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2643 += tmp_2644;
   }
   result += tmp_2643;
   std::complex<double> tmp_2645;
   std::complex<double> tmp_2646;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2647;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2647 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2646 += tmp_2647;
   }
   tmp_2645 += tmp_2646;
   result += (3) * tmp_2645;
   std::complex<double> tmp_2648;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2649;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2649 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUAhconjSvSv(gO2,gI1,
            gI2))*CpUAhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2648 += tmp_2649;
   }
   result += tmp_2648;
   std::complex<double> tmp_2650;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2651;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2651 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUAhSvSv(gO2,gI1,gI2))
            *CpUAhSvSv(gO1,gI1,gI2);
      }
      tmp_2650 += tmp_2651;
   }
   result += tmp_2650;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2654;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2654 += (Conj(CpUAhFvFvPL(gO2,gI1,gI2))*CpUAhFvFvPL(gO1,gI1,
            gI2) + Conj(CpUAhFvFvPR(gO2,gI1,gI2))*CpUAhFvFvPR(gO1,gI1,gI2))*G0(p,
            MFv(gI1),MFv(gI2));
      }
      tmp_2653 += tmp_2654;
   }
   tmp_2652 += tmp_2653;
   result += (0.5) * tmp_2652;
   std::complex<double> tmp_2655;
   std::complex<double> tmp_2656;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2657;
      std::complex<double> tmp_2658;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2658 += B0(p,MFv(gI1),MFv(gI2))*(Conj(CpUAhFvFvPR(gO2,gI1,
            gI2))*CpUAhFvFvPL(gO1,gI1,gI2) + Conj(CpUAhFvFvPL(gO2,gI1,gI2))*
            CpUAhFvFvPR(gO1,gI1,gI2))*MFv(gI2);
      }
      tmp_2657 += tmp_2658;
      tmp_2656 += (MFv(gI1)) * tmp_2657;
   }
   tmp_2655 += tmp_2656;
   result += (-1) * tmp_2655;
   std::complex<double> tmp_2659;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2659 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2659;
   std::complex<double> tmp_2660;
   std::complex<double> tmp_2661;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2661 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2660 += tmp_2661;
   result += (2) * tmp_2660;

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
   std::complex<double> tmp_2662;
   std::complex<double> tmp_2663;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2663 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2662 += tmp_2663;
   result += (-0.5) * tmp_2662;
   std::complex<double> tmp_2664;
   std::complex<double> tmp_2665;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2665 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2664 += tmp_2665;
   result += (-1) * tmp_2664;
   std::complex<double> tmp_2666;
   std::complex<double> tmp_2667;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2667 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2666 += tmp_2667;
   result += (-0.5) * tmp_2666;
   std::complex<double> tmp_2668;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2669;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2669 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2668 += tmp_2669;
   }
   result += tmp_2668;
   std::complex<double> tmp_2670;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2671;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2671 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2670 += tmp_2671;
   }
   result += tmp_2670;
   std::complex<double> tmp_2672;
   std::complex<double> tmp_2673;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2674;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2674 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2673 += tmp_2674;
   }
   tmp_2672 += tmp_2673;
   result += (3) * tmp_2672;
   std::complex<double> tmp_2675;
   std::complex<double> tmp_2676;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2677;
      std::complex<double> tmp_2678;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2678 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2677 += tmp_2678;
      tmp_2676 += (MFu(gI1)) * tmp_2677;
   }
   tmp_2675 += tmp_2676;
   result += (-6) * tmp_2675;
   std::complex<double> tmp_2679;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2680;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2680 += (Conj(CpconjUHpmChiChaPL(gO2,gI1,gI2))*
            CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHpmChiChaPR(gO2,gI1,gI2)
            )*CpconjUHpmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_2679 += tmp_2680;
   }
   result += tmp_2679;
   std::complex<double> tmp_2681;
   std::complex<double> tmp_2682;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2683;
      std::complex<double> tmp_2684;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2684 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHpmChiChaPR(
            gO2,gI1,gI2))*CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmChiChaPL(gO2,gI1,gI2))*CpconjUHpmChiChaPR(gO1,gI1,gI2))*MCha
            (gI2);
      }
      tmp_2683 += tmp_2684;
      tmp_2682 += (MChi(gI1)) * tmp_2683;
   }
   tmp_2681 += tmp_2682;
   result += (-2) * tmp_2681;
   std::complex<double> tmp_2685;
   std::complex<double> tmp_2686;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2686 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2685 += tmp_2686;
   result += (-3) * tmp_2685;
   std::complex<double> tmp_2687;
   std::complex<double> tmp_2688;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2688 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2687 += tmp_2688;
   result += (-1) * tmp_2687;
   std::complex<double> tmp_2689;
   std::complex<double> tmp_2690;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2690 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2689 += tmp_2690;
   result += (-3) * tmp_2689;
   std::complex<double> tmp_2691;
   std::complex<double> tmp_2692;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2692 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2691 += tmp_2692;
   result += (-1) * tmp_2691;
   std::complex<double> tmp_2693;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2694;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2694 += (Conj(CpconjUHpmFvFePL(gO2,gI1,gI2))*
            CpconjUHpmFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmFvFePR(gO2,gI1,gI2))*
            CpconjUHpmFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2693 += tmp_2694;
   }
   result += tmp_2693;
   std::complex<double> tmp_2695;
   std::complex<double> tmp_2696;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2697;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2697 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2696 += tmp_2697;
   }
   tmp_2695 += tmp_2696;
   result += (3) * tmp_2695;
   std::complex<double> tmp_2698;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2699;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2699 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2698 += tmp_2699;
   }
   result += tmp_2698;
   std::complex<double> tmp_2700;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2701;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2701 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmSvSe(gO2,gI1,
            gI2))*CpconjUHpmSvSe(gO1,gI1,gI2);
      }
      tmp_2700 += tmp_2701;
   }
   result += tmp_2700;
   std::complex<double> tmp_2702;
   std::complex<double> tmp_2703;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2704;
      std::complex<double> tmp_2705;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2705 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmFvFePR(gO2,
            gI1,gI2))*CpconjUHpmFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmFvFePL(gO2,
            gI1,gI2))*CpconjUHpmFvFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2704 += tmp_2705;
      tmp_2703 += (MFv(gI1)) * tmp_2704;
   }
   tmp_2702 += tmp_2703;
   result += (-2) * tmp_2702;
   std::complex<double> tmp_2706;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2706 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2706;
   std::complex<double> tmp_2707;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2707 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2707;
   std::complex<double> tmp_2708;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2708 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2708;
   std::complex<double> tmp_2709;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2709 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2709;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2710;
   std::complex<double> tmp_2711;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2711 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2710 += tmp_2711;
   result += (0.5) * tmp_2710;
   std::complex<double> tmp_2712;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2712 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2712;
   std::complex<double> tmp_2713;
   std::complex<double> tmp_2714;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2714 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2713 += tmp_2714;
   result += (0.5) * tmp_2713;
   std::complex<double> tmp_2715;
   std::complex<double> tmp_2716;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2717;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2717 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2716 += tmp_2717;
   }
   tmp_2715 += tmp_2716;
   result += (-4) * tmp_2715;
   std::complex<double> tmp_2718;
   std::complex<double> tmp_2719;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2720;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2720 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2719 += tmp_2720;
   }
   tmp_2718 += tmp_2719;
   result += (-4) * tmp_2718;
   std::complex<double> tmp_2721;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2722;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2722 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2722 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2721 += tmp_2722;
   }
   result += tmp_2721;
   std::complex<double> tmp_2723;
   std::complex<double> tmp_2724;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2725;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2725 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2725 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2724 += tmp_2725;
   }
   tmp_2723 += tmp_2724;
   result += (3) * tmp_2723;
   std::complex<double> tmp_2726;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2727;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2727 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2727 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2726 += tmp_2727;
   }
   result += tmp_2726;
   std::complex<double> tmp_2728;
   std::complex<double> tmp_2729;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2730;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2730 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2730 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2729 += tmp_2730;
   }
   tmp_2728 += tmp_2729;
   result += (3) * tmp_2728;
   std::complex<double> tmp_2731;
   std::complex<double> tmp_2732;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2733;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2733 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2733 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2732 += tmp_2733;
   }
   tmp_2731 += tmp_2732;
   result += (0.5) * tmp_2731;
   std::complex<double> tmp_2734;
   std::complex<double> tmp_2735;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2735 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2734 += tmp_2735;
   result += (3) * tmp_2734;
   std::complex<double> tmp_2736;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2736 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2736;
   std::complex<double> tmp_2737;
   std::complex<double> tmp_2738;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2738 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2737 += tmp_2738;
   result += (3) * tmp_2737;
   std::complex<double> tmp_2739;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2739 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2739;
   std::complex<double> tmp_2740;
   std::complex<double> tmp_2741;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2742;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2742 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2741 += tmp_2742;
   }
   tmp_2740 += tmp_2741;
   result += (-12) * tmp_2740;
   std::complex<double> tmp_2743;
   std::complex<double> tmp_2744;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2745;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2745 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2744 += tmp_2745;
   }
   tmp_2743 += tmp_2744;
   result += (-4) * tmp_2743;
   std::complex<double> tmp_2746;
   std::complex<double> tmp_2747;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2748;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2748 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2747 += tmp_2748;
   }
   tmp_2746 += tmp_2747;
   result += (-12) * tmp_2746;
   std::complex<double> tmp_2749;
   std::complex<double> tmp_2750;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2751;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2751 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2750 += tmp_2751;
   }
   tmp_2749 += tmp_2750;
   result += (-4) * tmp_2749;
   std::complex<double> tmp_2752;
   std::complex<double> tmp_2753;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2754;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2754 += (AbsSqr(CpVZFvFvPL(gI1,gI2)) + AbsSqr(CpVZFvFvPR(gI1
            ,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2754 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZFvFvPL(gI1,gI2))*CpVZFvFvPR(gI1,gI2));
      }
      tmp_2753 += tmp_2754;
   }
   tmp_2752 += tmp_2753;
   result += (0.5) * tmp_2752;
   std::complex<double> tmp_2755;
   std::complex<double> tmp_2756;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2756 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2755 += tmp_2756;
   result += (2) * tmp_2755;
   std::complex<double> tmp_2757;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2757 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2757;
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
   std::complex<double> tmp_2758;
   std::complex<double> tmp_2759;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2759 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_2758 += tmp_2759;
   result += (0.5) * tmp_2758;
   std::complex<double> tmp_2760;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2760 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_2760;
   std::complex<double> tmp_2761;
   std::complex<double> tmp_2762;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2762 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2761 += tmp_2762;
   result += (0.5) * tmp_2761;
   std::complex<double> tmp_2763;
   std::complex<double> tmp_2764;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2765;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2765 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_2764 += tmp_2765;
   }
   tmp_2763 += tmp_2764;
   result += (-4) * tmp_2763;
   std::complex<double> tmp_2766;
   std::complex<double> tmp_2767;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2768;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2768 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_2767 += tmp_2768;
   }
   tmp_2766 += tmp_2767;
   result += (-4) * tmp_2766;
   std::complex<double> tmp_2769;
   std::complex<double> tmp_2770;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2771;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2771 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_2771 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_2770 += tmp_2771;
   }
   tmp_2769 += tmp_2770;
   result += (3) * tmp_2769;
   std::complex<double> tmp_2772;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2773;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2773 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_2773 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_2772 += tmp_2773;
   }
   result += tmp_2772;
   std::complex<double> tmp_2774;
   std::complex<double> tmp_2775;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2775 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2774 += tmp_2775;
   result += (3) * tmp_2774;
   std::complex<double> tmp_2776;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2776 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2776;
   std::complex<double> tmp_2777;
   std::complex<double> tmp_2778;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2778 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2777 += tmp_2778;
   result += (3) * tmp_2777;
   std::complex<double> tmp_2779;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2779 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2779;
   std::complex<double> tmp_2780;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2781;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2781 += (AbsSqr(CpconjVWmFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_2781 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmFvFePL(gI1,gI2))*CpconjVWmFvFePR(gI1,gI2));
      }
      tmp_2780 += tmp_2781;
   }
   result += tmp_2780;
   std::complex<double> tmp_2782;
   std::complex<double> tmp_2783;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2784;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2784 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2783 += tmp_2784;
   }
   tmp_2782 += tmp_2783;
   result += (-12) * tmp_2782;
   std::complex<double> tmp_2785;
   std::complex<double> tmp_2786;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2787;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2787 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2786 += tmp_2787;
   }
   tmp_2785 += tmp_2786;
   result += (-4) * tmp_2785;
   std::complex<double> tmp_2788;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2788 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_2788;
   std::complex<double> tmp_2789;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2789 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2789;
   std::complex<double> tmp_2790;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2790 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_2790;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2791;
   std::complex<double> tmp_2792;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2793;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2793 += B0(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2792 += tmp_2793;
   }
   tmp_2791 += tmp_2792;
   result += (2) * tmp_2791;
   std::complex<double> tmp_2794;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2795;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2795 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2794 += tmp_2795;
   }
   result += tmp_2794;
   std::complex<double> tmp_2796;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2797;
      std::complex<double> tmp_2798;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2798 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2797 += tmp_2798;
      tmp_2796 += (MChi(gI1)) * tmp_2797;
   }
   result += tmp_2796;
   std::complex<double> tmp_2799;
   std::complex<double> tmp_2800;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2801;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2801 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2800 += tmp_2801;
   }
   tmp_2799 += tmp_2800;
   result += (6) * tmp_2799;
   std::complex<double> tmp_2802;
   std::complex<double> tmp_2803;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2804;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2804 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2803 += tmp_2804;
   }
   tmp_2802 += tmp_2803;
   result += (2) * tmp_2802;
   std::complex<double> tmp_2805;
   std::complex<double> tmp_2806;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2807;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2807 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2806 += tmp_2807;
   }
   tmp_2805 += tmp_2806;
   result += (6) * tmp_2805;
   std::complex<double> tmp_2808;
   std::complex<double> tmp_2809;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2810;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2810 += B0(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2809 += tmp_2810;
   }
   tmp_2808 += tmp_2809;
   result += (2) * tmp_2808;
   std::complex<double> tmp_2811;
   std::complex<double> tmp_2812;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2812 += B0(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2811 += tmp_2812;
   result += (-8) * tmp_2811;
   std::complex<double> tmp_2813;
   std::complex<double> tmp_2814;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2814 += B0(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2813 += tmp_2814;
   result += (-4) * tmp_2813;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2815;
   std::complex<double> tmp_2816;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2817;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2817 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPR(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2);
      }
      tmp_2816 += tmp_2817;
   }
   tmp_2815 += tmp_2816;
   result += (-1) * tmp_2815;
   std::complex<double> tmp_2818;
   std::complex<double> tmp_2819;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2820;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2820 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_2819 += tmp_2820;
   }
   tmp_2818 += tmp_2819;
   result += (-0.5) * tmp_2818;
   std::complex<double> tmp_2821;
   std::complex<double> tmp_2822;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2823;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2823 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2822 += tmp_2823;
   }
   tmp_2821 += tmp_2822;
   result += (-0.5) * tmp_2821;
   std::complex<double> tmp_2824;
   std::complex<double> tmp_2825;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2826;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2826 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPR(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_2825 += tmp_2826;
   }
   tmp_2824 += tmp_2825;
   result += (-3) * tmp_2824;
   std::complex<double> tmp_2827;
   std::complex<double> tmp_2828;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2829;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2829 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePR(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_2828 += tmp_2829;
   }
   tmp_2827 += tmp_2828;
   result += (-1) * tmp_2827;
   std::complex<double> tmp_2830;
   std::complex<double> tmp_2831;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2832;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2832 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPR(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_2831 += tmp_2832;
   }
   tmp_2830 += tmp_2831;
   result += (-3) * tmp_2830;
   std::complex<double> tmp_2833;
   std::complex<double> tmp_2834;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2835;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2835 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPR(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2);
      }
      tmp_2834 += tmp_2835;
   }
   tmp_2833 += tmp_2834;
   result += (-1) * tmp_2833;
   std::complex<double> tmp_2836;
   std::complex<double> tmp_2837;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2837 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPL(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_2836 += tmp_2837;
   result += (-2) * tmp_2836;
   std::complex<double> tmp_2838;
   std::complex<double> tmp_2839;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2839 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_2838 += tmp_2839;
   result += (-1) * tmp_2838;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2840;
   std::complex<double> tmp_2841;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2842;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2842 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPL(gO1,gI1,gI2);
      }
      tmp_2841 += tmp_2842;
   }
   tmp_2840 += tmp_2841;
   result += (-1) * tmp_2840;
   std::complex<double> tmp_2843;
   std::complex<double> tmp_2844;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2845;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2845 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
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
         tmp_2848 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_2847 += tmp_2848;
   }
   tmp_2846 += tmp_2847;
   result += (-0.5) * tmp_2846;
   std::complex<double> tmp_2849;
   std::complex<double> tmp_2850;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2851;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2851 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_2850 += tmp_2851;
   }
   tmp_2849 += tmp_2850;
   result += (-3) * tmp_2849;
   std::complex<double> tmp_2852;
   std::complex<double> tmp_2853;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2854;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2854 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_2853 += tmp_2854;
   }
   tmp_2852 += tmp_2853;
   result += (-1) * tmp_2852;
   std::complex<double> tmp_2855;
   std::complex<double> tmp_2856;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2857;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2857 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_2856 += tmp_2857;
   }
   tmp_2855 += tmp_2856;
   result += (-3) * tmp_2855;
   std::complex<double> tmp_2858;
   std::complex<double> tmp_2859;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2860;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2860 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPL(gO1,gI1,gI2);
      }
      tmp_2859 += tmp_2860;
   }
   tmp_2858 += tmp_2859;
   result += (-1) * tmp_2858;
   std::complex<double> tmp_2861;
   std::complex<double> tmp_2862;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2862 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_2861 += tmp_2862;
   result += (-2) * tmp_2861;
   std::complex<double> tmp_2863;
   std::complex<double> tmp_2864;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2864 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_2863 += tmp_2864;
   result += (-1) * tmp_2863;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2865;
   std::complex<double> tmp_2866;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2867;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2867 += B0(p,MFe(gI2),MHpm(gI1))*Conj(CpUFvconjHpmFePL(gO2,
            gI1,gI2))*CpUFvconjHpmFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2866 += tmp_2867;
   }
   tmp_2865 += tmp_2866;
   result += (2) * tmp_2865;
   std::complex<double> tmp_2868;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2869;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2869 += B0(p,MFv(gI2),Mhh(gI1))*Conj(CpUFvhhFvPL(gO2,gI1,gI2
            ))*CpUFvhhFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2868 += tmp_2869;
   }
   result += tmp_2868;
   std::complex<double> tmp_2870;
   std::complex<double> tmp_2871;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2872;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2872 += B0(p,MCha(gI2),MSe(gI1))*Conj(CpUFvconjSeChaPL(gO2,
            gI1,gI2))*CpUFvconjSeChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2871 += tmp_2872;
   }
   tmp_2870 += tmp_2871;
   result += (2) * tmp_2870;
   std::complex<double> tmp_2873;
   std::complex<double> tmp_2874;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2875;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2875 += B0(p,MChi(gI2),MSv(gI1))*Conj(CpUFvconjSvChiPL(gO2,
            gI1,gI2))*CpUFvconjSvChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2874 += tmp_2875;
   }
   tmp_2873 += tmp_2874;
   result += (2) * tmp_2873;
   std::complex<double> tmp_2876;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2877;
      std::complex<double> tmp_2878;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2878 += B0(p,MFv(gI1),MAh(gI2))*Conj(CpUFvFvAhPL(gO2,gI1,gI2
            ))*CpUFvFvAhPR(gO1,gI1,gI2);
      }
      tmp_2877 += tmp_2878;
      tmp_2876 += (MFv(gI1)) * tmp_2877;
   }
   result += tmp_2876;
   std::complex<double> tmp_2879;
   std::complex<double> tmp_2880;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2880 += B0(p,MFe(gI2),MVWm)*Conj(CpUFvconjVWmFePR(gO2,gI2))*
         CpUFvconjVWmFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2879 += tmp_2880;
   result += (-8) * tmp_2879;
   std::complex<double> tmp_2881;
   std::complex<double> tmp_2882;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2882 += B0(p,MFv(gI2),MVZ)*Conj(CpUFvVZFvPR(gO2,gI2))*CpUFvVZFvPL(
         gO1,gI2)*MFv(gI2);
   }
   tmp_2881 += tmp_2882;
   result += (-4) * tmp_2881;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2883;
   std::complex<double> tmp_2884;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2885;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2885 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpUFvconjHpmFePR(gO2,
            gI1,gI2))*CpUFvconjHpmFePR(gO1,gI1,gI2);
      }
      tmp_2884 += tmp_2885;
   }
   tmp_2883 += tmp_2884;
   result += (-1) * tmp_2883;
   std::complex<double> tmp_2886;
   std::complex<double> tmp_2887;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2888;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2888 += B1(p,MFv(gI2),Mhh(gI1))*Conj(CpUFvhhFvPR(gO2,gI1,gI2
            ))*CpUFvhhFvPR(gO1,gI1,gI2);
      }
      tmp_2887 += tmp_2888;
   }
   tmp_2886 += tmp_2887;
   result += (-0.5) * tmp_2886;
   std::complex<double> tmp_2889;
   std::complex<double> tmp_2890;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2891;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2891 += B1(p,MCha(gI2),MSe(gI1))*Conj(CpUFvconjSeChaPR(gO2,
            gI1,gI2))*CpUFvconjSeChaPR(gO1,gI1,gI2);
      }
      tmp_2890 += tmp_2891;
   }
   tmp_2889 += tmp_2890;
   result += (-1) * tmp_2889;
   std::complex<double> tmp_2892;
   std::complex<double> tmp_2893;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2894;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2894 += B1(p,MFv(gI1),MAh(gI2))*Conj(CpUFvFvAhPR(gO2,gI1,gI2
            ))*CpUFvFvAhPR(gO1,gI1,gI2);
      }
      tmp_2893 += tmp_2894;
   }
   tmp_2892 += tmp_2893;
   result += (-0.5) * tmp_2892;
   std::complex<double> tmp_2895;
   std::complex<double> tmp_2896;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2897;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2897 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpUFvconjSvChiPR(gO2,
            gI1,gI2))*CpUFvconjSvChiPR(gO1,gI1,gI2);
      }
      tmp_2896 += tmp_2897;
   }
   tmp_2895 += tmp_2896;
   result += (-1) * tmp_2895;
   std::complex<double> tmp_2898;
   std::complex<double> tmp_2899;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2899 += B1(p,MFe(gI2),MVWm)*Conj(CpUFvconjVWmFePL(gO2,gI2))*
         CpUFvconjVWmFePL(gO1,gI2);
   }
   tmp_2898 += tmp_2899;
   result += (-2) * tmp_2898;
   std::complex<double> tmp_2900;
   std::complex<double> tmp_2901;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2901 += B1(p,MFv(gI2),MVZ)*Conj(CpUFvVZFvPL(gO2,gI2))*CpUFvVZFvPL(
         gO1,gI2);
   }
   tmp_2900 += tmp_2901;
   result += (-1) * tmp_2900;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2902;
   std::complex<double> tmp_2903;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2904;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2904 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpUFvconjHpmFePL(gO2,
            gI1,gI2))*CpUFvconjHpmFePL(gO1,gI1,gI2);
      }
      tmp_2903 += tmp_2904;
   }
   tmp_2902 += tmp_2903;
   result += (-1) * tmp_2902;
   std::complex<double> tmp_2905;
   std::complex<double> tmp_2906;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2907;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2907 += B1(p,MFv(gI2),Mhh(gI1))*Conj(CpUFvhhFvPL(gO2,gI1,gI2
            ))*CpUFvhhFvPL(gO1,gI1,gI2);
      }
      tmp_2906 += tmp_2907;
   }
   tmp_2905 += tmp_2906;
   result += (-0.5) * tmp_2905;
   std::complex<double> tmp_2908;
   std::complex<double> tmp_2909;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2910;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2910 += B1(p,MCha(gI2),MSe(gI1))*Conj(CpUFvconjSeChaPL(gO2,
            gI1,gI2))*CpUFvconjSeChaPL(gO1,gI1,gI2);
      }
      tmp_2909 += tmp_2910;
   }
   tmp_2908 += tmp_2909;
   result += (-1) * tmp_2908;
   std::complex<double> tmp_2911;
   std::complex<double> tmp_2912;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2913;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2913 += B1(p,MFv(gI1),MAh(gI2))*Conj(CpUFvFvAhPL(gO2,gI1,gI2
            ))*CpUFvFvAhPL(gO1,gI1,gI2);
      }
      tmp_2912 += tmp_2913;
   }
   tmp_2911 += tmp_2912;
   result += (-0.5) * tmp_2911;
   std::complex<double> tmp_2914;
   std::complex<double> tmp_2915;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2916;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2916 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpUFvconjSvChiPL(gO2,
            gI1,gI2))*CpUFvconjSvChiPL(gO1,gI1,gI2);
      }
      tmp_2915 += tmp_2916;
   }
   tmp_2914 += tmp_2915;
   result += (-1) * tmp_2914;
   std::complex<double> tmp_2917;
   std::complex<double> tmp_2918;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2918 += B1(p,MFe(gI2),MVWm)*Conj(CpUFvconjVWmFePR(gO2,gI2))*
         CpUFvconjVWmFePR(gO1,gI2);
   }
   tmp_2917 += tmp_2918;
   result += (-2) * tmp_2917;
   std::complex<double> tmp_2919;
   std::complex<double> tmp_2920;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2920 += B1(p,MFv(gI2),MVZ)*Conj(CpUFvVZFvPR(gO2,gI2))*CpUFvVZFvPR(
         gO1,gI2);
   }
   tmp_2919 += tmp_2920;
   result += (-1) * tmp_2919;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2921;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2922;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2922 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2921 += tmp_2922;
   }
   result += tmp_2921;
   std::complex<double> tmp_2923;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2924;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2924 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2923 += tmp_2924;
   }
   result += tmp_2923;
   std::complex<double> tmp_2925;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2926;
      std::complex<double> tmp_2927;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2927 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2926 += tmp_2927;
      tmp_2925 += (MCha(gI1)) * tmp_2926;
   }
   result += tmp_2925;
   std::complex<double> tmp_2928;
   std::complex<double> tmp_2929;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2930;
      std::complex<double> tmp_2931;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2931 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2930 += tmp_2931;
      tmp_2929 += (MFu(gI1)) * tmp_2930;
   }
   tmp_2928 += tmp_2929;
   result += (3) * tmp_2928;
   std::complex<double> tmp_2932;
   std::complex<double> tmp_2933;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2934;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2934 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2933 += tmp_2934;
   }
   tmp_2932 += tmp_2933;
   result += (3) * tmp_2932;
   std::complex<double> tmp_2935;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2936;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2936 += B0(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2935 += tmp_2936;
   }
   result += tmp_2935;
   std::complex<double> tmp_2937;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2938;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2938 += B0(p,MFv(gI2),MSe(gI1))*Conj(CpbarUChaSeFvPL(gO2,gI1
            ,gI2))*CpbarUChaSeFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2937 += tmp_2938;
   }
   result += tmp_2937;
   std::complex<double> tmp_2939;
   std::complex<double> tmp_2940;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2940 += B0(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2939 += tmp_2940;
   result += (-4) * tmp_2939;
   std::complex<double> tmp_2941;
   std::complex<double> tmp_2942;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2942 += B0(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2941 += tmp_2942;
   result += (-4) * tmp_2941;
   std::complex<double> tmp_2943;
   std::complex<double> tmp_2944;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2944 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2943 += tmp_2944;
   result += (-4) * tmp_2943;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2945;
   std::complex<double> tmp_2946;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2947;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2947 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2946 += tmp_2947;
   }
   tmp_2945 += tmp_2946;
   result += (-0.5) * tmp_2945;
   std::complex<double> tmp_2948;
   std::complex<double> tmp_2949;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2950;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2950 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_2949 += tmp_2950;
   }
   tmp_2948 += tmp_2949;
   result += (-0.5) * tmp_2948;
   std::complex<double> tmp_2951;
   std::complex<double> tmp_2952;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2953;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2953 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPR(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2);
      }
      tmp_2952 += tmp_2953;
   }
   tmp_2951 += tmp_2952;
   result += (-0.5) * tmp_2951;
   std::complex<double> tmp_2954;
   std::complex<double> tmp_2955;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2956;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2956 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPR(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2955 += tmp_2956;
   }
   tmp_2954 += tmp_2955;
   result += (-1.5) * tmp_2954;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2959;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2959 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPR(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2);
      }
      tmp_2958 += tmp_2959;
   }
   tmp_2957 += tmp_2958;
   result += (-1.5) * tmp_2957;
   std::complex<double> tmp_2960;
   std::complex<double> tmp_2961;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2962;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2962 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePR(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2);
      }
      tmp_2961 += tmp_2962;
   }
   tmp_2960 += tmp_2961;
   result += (-0.5) * tmp_2960;
   std::complex<double> tmp_2963;
   std::complex<double> tmp_2964;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2965;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2965 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUChaSeFvPR(gO2,gI1
            ,gI2))*CpbarUChaSeFvPR(gO1,gI1,gI2);
      }
      tmp_2964 += tmp_2965;
   }
   tmp_2963 += tmp_2964;
   result += (-0.5) * tmp_2963;
   std::complex<double> tmp_2966;
   std::complex<double> tmp_2967;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2967 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_2966 += tmp_2967;
   result += (-1) * tmp_2966;
   std::complex<double> tmp_2968;
   std::complex<double> tmp_2969;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2969 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPL(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_2968 += tmp_2969;
   result += (-1) * tmp_2968;
   std::complex<double> tmp_2970;
   std::complex<double> tmp_2971;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2971 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPL(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_2970 += tmp_2971;
   result += (-1) * tmp_2970;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2972;
   std::complex<double> tmp_2973;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2974;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2974 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_2973 += tmp_2974;
   }
   tmp_2972 += tmp_2973;
   result += (-0.5) * tmp_2972;
   std::complex<double> tmp_2975;
   std::complex<double> tmp_2976;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2977;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2977 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_2976 += tmp_2977;
   }
   tmp_2975 += tmp_2976;
   result += (-0.5) * tmp_2975;
   std::complex<double> tmp_2978;
   std::complex<double> tmp_2979;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2980;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2980 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPL(gO1,gI1,gI2);
      }
      tmp_2979 += tmp_2980;
   }
   tmp_2978 += tmp_2979;
   result += (-0.5) * tmp_2978;
   std::complex<double> tmp_2981;
   std::complex<double> tmp_2982;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2983;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2983 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2);
      }
      tmp_2982 += tmp_2983;
   }
   tmp_2981 += tmp_2982;
   result += (-1.5) * tmp_2981;
   std::complex<double> tmp_2984;
   std::complex<double> tmp_2985;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2986;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2986 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPL(gO1,gI1,gI2);
      }
      tmp_2985 += tmp_2986;
   }
   tmp_2984 += tmp_2985;
   result += (-1.5) * tmp_2984;
   std::complex<double> tmp_2987;
   std::complex<double> tmp_2988;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2989;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2989 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePL(gO1,gI1,gI2);
      }
      tmp_2988 += tmp_2989;
   }
   tmp_2987 += tmp_2988;
   result += (-0.5) * tmp_2987;
   std::complex<double> tmp_2990;
   std::complex<double> tmp_2991;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2992;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2992 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUChaSeFvPL(gO2,gI1
            ,gI2))*CpbarUChaSeFvPL(gO1,gI1,gI2);
      }
      tmp_2991 += tmp_2992;
   }
   tmp_2990 += tmp_2991;
   result += (-0.5) * tmp_2990;
   std::complex<double> tmp_2993;
   std::complex<double> tmp_2994;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2994 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_2993 += tmp_2994;
   result += (-1) * tmp_2993;
   std::complex<double> tmp_2995;
   std::complex<double> tmp_2996;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2996 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_2995 += tmp_2996;
   result += (-1) * tmp_2995;
   std::complex<double> tmp_2997;
   std::complex<double> tmp_2998;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2998 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_2997 += tmp_2998;
   result += (-1) * tmp_2997;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2999;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3000;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3000 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2999 += tmp_3000;
   }
   result += tmp_2999;
   std::complex<double> tmp_3001;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3002;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3002 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3001 += tmp_3002;
   }
   result += tmp_3001;
   std::complex<double> tmp_3003;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3004;
      std::complex<double> tmp_3005;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3005 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3004 += tmp_3005;
      tmp_3003 += (MFe(gI1)) * tmp_3004;
   }
   result += tmp_3003;
   std::complex<double> tmp_3006;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3007;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3007 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3006 += tmp_3007;
   }
   result += tmp_3006;
   std::complex<double> tmp_3008;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3009;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3009 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3008 += tmp_3009;
   }
   result += tmp_3008;
   std::complex<double> tmp_3010;
   std::complex<double> tmp_3011;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3011 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3010 += tmp_3011;
   result += (-4) * tmp_3010;
   std::complex<double> tmp_3012;
   std::complex<double> tmp_3013;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3013 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3012 += tmp_3013;
   result += (-4) * tmp_3012;
   std::complex<double> tmp_3014;
   std::complex<double> tmp_3015;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3015 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3014 += tmp_3015;
   result += (-4) * tmp_3014;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3016;
   std::complex<double> tmp_3017;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3018;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3018 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_3017 += tmp_3018;
   }
   tmp_3016 += tmp_3017;
   result += (-0.5) * tmp_3016;
   std::complex<double> tmp_3019;
   std::complex<double> tmp_3020;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3021;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3021 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3020 += tmp_3021;
   }
   tmp_3019 += tmp_3020;
   result += (-0.5) * tmp_3019;
   std::complex<double> tmp_3022;
   std::complex<double> tmp_3023;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3024;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3024 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3023 += tmp_3024;
   }
   tmp_3022 += tmp_3023;
   result += (-0.5) * tmp_3022;
   std::complex<double> tmp_3025;
   std::complex<double> tmp_3026;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3027;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3027 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPR(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3026 += tmp_3027;
   }
   tmp_3025 += tmp_3026;
   result += (-0.5) * tmp_3025;
   std::complex<double> tmp_3028;
   std::complex<double> tmp_3029;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3030;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3030 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3029 += tmp_3030;
   }
   tmp_3028 += tmp_3029;
   result += (-0.5) * tmp_3028;
   std::complex<double> tmp_3031;
   std::complex<double> tmp_3032;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3032 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_3031 += tmp_3032;
   result += (-1) * tmp_3031;
   std::complex<double> tmp_3033;
   std::complex<double> tmp_3034;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3034 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_3033 += tmp_3034;
   result += (-1) * tmp_3033;
   std::complex<double> tmp_3035;
   std::complex<double> tmp_3036;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3036 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_3035 += tmp_3036;
   result += (-1) * tmp_3035;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3037;
   std::complex<double> tmp_3038;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3039;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3039 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_3038 += tmp_3039;
   }
   tmp_3037 += tmp_3038;
   result += (-0.5) * tmp_3037;
   std::complex<double> tmp_3040;
   std::complex<double> tmp_3041;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3042;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3042 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3041 += tmp_3042;
   }
   tmp_3040 += tmp_3041;
   result += (-0.5) * tmp_3040;
   std::complex<double> tmp_3043;
   std::complex<double> tmp_3044;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3045;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3045 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3044 += tmp_3045;
   }
   tmp_3043 += tmp_3044;
   result += (-0.5) * tmp_3043;
   std::complex<double> tmp_3046;
   std::complex<double> tmp_3047;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3048;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3048 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3047 += tmp_3048;
   }
   tmp_3046 += tmp_3047;
   result += (-0.5) * tmp_3046;
   std::complex<double> tmp_3049;
   std::complex<double> tmp_3050;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3051;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3051 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3050 += tmp_3051;
   }
   tmp_3049 += tmp_3050;
   result += (-0.5) * tmp_3049;
   std::complex<double> tmp_3052;
   std::complex<double> tmp_3053;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3053 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_3052 += tmp_3053;
   result += (-1) * tmp_3052;
   std::complex<double> tmp_3054;
   std::complex<double> tmp_3055;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3055 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_3054 += tmp_3055;
   result += (-1) * tmp_3054;
   std::complex<double> tmp_3056;
   std::complex<double> tmp_3057;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3057 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_3056 += tmp_3057;
   result += (-1) * tmp_3056;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3058;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3059;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3059 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3058 += tmp_3059;
   }
   result += tmp_3058;
   std::complex<double> tmp_3060;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3061;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3061 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3060 += tmp_3061;
   }
   result += tmp_3060;
   std::complex<double> tmp_3062;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3063;
      std::complex<double> tmp_3064;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3064 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3063 += tmp_3064;
      tmp_3062 += (MFd(gI1)) * tmp_3063;
   }
   result += tmp_3062;
   std::complex<double> tmp_3065;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3066;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3066 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3065 += tmp_3066;
   }
   result += tmp_3065;
   std::complex<double> tmp_3067;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3068;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3068 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3067 += tmp_3068;
   }
   result += tmp_3067;
   std::complex<double> tmp_3069;
   std::complex<double> tmp_3070;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3070 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3069 += tmp_3070;
   result += (-5.333333333333333) * tmp_3069;
   std::complex<double> tmp_3071;
   std::complex<double> tmp_3072;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3072 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3071 += tmp_3072;
   result += (-4) * tmp_3071;
   std::complex<double> tmp_3073;
   std::complex<double> tmp_3074;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3074 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3073 += tmp_3074;
   result += (-4) * tmp_3073;
   std::complex<double> tmp_3075;
   std::complex<double> tmp_3076;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3076 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3075 += tmp_3076;
   result += (-4) * tmp_3075;
   std::complex<double> tmp_3077;
   std::complex<double> tmp_3078;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3078 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3077 += tmp_3078;
   result += (1.3333333333333333*MGlu) * tmp_3077;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3079;
   std::complex<double> tmp_3080;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3081;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3081 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3080 += tmp_3081;
   }
   tmp_3079 += tmp_3080;
   result += (-0.5) * tmp_3079;
   std::complex<double> tmp_3082;
   std::complex<double> tmp_3083;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3084;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3084 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3083 += tmp_3084;
   }
   tmp_3082 += tmp_3083;
   result += (-0.5) * tmp_3082;
   std::complex<double> tmp_3085;
   std::complex<double> tmp_3086;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3087;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3087 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3086 += tmp_3087;
   }
   tmp_3085 += tmp_3086;
   result += (-0.5) * tmp_3085;
   std::complex<double> tmp_3088;
   std::complex<double> tmp_3089;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3089 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3088 += tmp_3089;
   result += (-0.6666666666666666) * tmp_3088;
   std::complex<double> tmp_3090;
   std::complex<double> tmp_3091;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3092;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3092 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPR(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3091 += tmp_3092;
   }
   tmp_3090 += tmp_3091;
   result += (-0.5) * tmp_3090;
   std::complex<double> tmp_3093;
   std::complex<double> tmp_3094;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3095;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3095 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3094 += tmp_3095;
   }
   tmp_3093 += tmp_3094;
   result += (-0.5) * tmp_3093;
   std::complex<double> tmp_3096;
   std::complex<double> tmp_3097;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3097 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_3096 += tmp_3097;
   result += (-1.3333333333333333) * tmp_3096;
   std::complex<double> tmp_3098;
   std::complex<double> tmp_3099;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3099 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_3098 += tmp_3099;
   result += (-1) * tmp_3098;
   std::complex<double> tmp_3100;
   std::complex<double> tmp_3101;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3101 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_3100 += tmp_3101;
   result += (-1) * tmp_3100;
   std::complex<double> tmp_3102;
   std::complex<double> tmp_3103;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3103 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_3102 += tmp_3103;
   result += (-1) * tmp_3102;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3104;
   std::complex<double> tmp_3105;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3106;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3106 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3105 += tmp_3106;
   }
   tmp_3104 += tmp_3105;
   result += (-0.5) * tmp_3104;
   std::complex<double> tmp_3107;
   std::complex<double> tmp_3108;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3109;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3109 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3108 += tmp_3109;
   }
   tmp_3107 += tmp_3108;
   result += (-0.5) * tmp_3107;
   std::complex<double> tmp_3110;
   std::complex<double> tmp_3111;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3112;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3112 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3111 += tmp_3112;
   }
   tmp_3110 += tmp_3111;
   result += (-0.5) * tmp_3110;
   std::complex<double> tmp_3113;
   std::complex<double> tmp_3114;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3114 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPL(gO1,gI1,1);
   }
   tmp_3113 += tmp_3114;
   result += (-0.6666666666666666) * tmp_3113;
   std::complex<double> tmp_3115;
   std::complex<double> tmp_3116;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3117;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3117 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_3116 += tmp_3117;
   }
   tmp_3115 += tmp_3116;
   result += (-0.5) * tmp_3115;
   std::complex<double> tmp_3118;
   std::complex<double> tmp_3119;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3120;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3120 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3119 += tmp_3120;
   }
   tmp_3118 += tmp_3119;
   result += (-0.5) * tmp_3118;
   std::complex<double> tmp_3121;
   std::complex<double> tmp_3122;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3122 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_3121 += tmp_3122;
   result += (-1.3333333333333333) * tmp_3121;
   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3124 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_3123 += tmp_3124;
   result += (-1) * tmp_3123;
   std::complex<double> tmp_3125;
   std::complex<double> tmp_3126;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3126 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_3125 += tmp_3126;
   result += (-1) * tmp_3125;
   std::complex<double> tmp_3127;
   std::complex<double> tmp_3128;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3128 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_3127 += tmp_3128;
   result += (-1) * tmp_3127;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3129;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3130;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3130 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3129 += tmp_3130;
   }
   result += tmp_3129;
   std::complex<double> tmp_3131;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3132;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3132 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3131 += tmp_3132;
   }
   result += tmp_3131;
   std::complex<double> tmp_3133;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3134;
      std::complex<double> tmp_3135;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3135 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3134 += tmp_3135;
      tmp_3133 += (MCha(gI1)) * tmp_3134;
   }
   result += tmp_3133;
   std::complex<double> tmp_3136;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3137;
      std::complex<double> tmp_3138;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3138 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3137 += tmp_3138;
      tmp_3136 += (MFu(gI1)) * tmp_3137;
   }
   result += tmp_3136;
   std::complex<double> tmp_3139;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3140;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3140 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3139 += tmp_3140;
   }
   result += tmp_3139;
   std::complex<double> tmp_3141;
   std::complex<double> tmp_3142;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3142 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3141 += tmp_3142;
   result += (-4) * tmp_3141;
   std::complex<double> tmp_3143;
   std::complex<double> tmp_3144;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3144 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3143 += tmp_3144;
   result += (-5.333333333333333) * tmp_3143;
   std::complex<double> tmp_3145;
   std::complex<double> tmp_3146;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3146 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3145 += tmp_3146;
   result += (-4) * tmp_3145;
   std::complex<double> tmp_3147;
   std::complex<double> tmp_3148;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3148 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3147 += tmp_3148;
   result += (-4) * tmp_3147;
   std::complex<double> tmp_3149;
   std::complex<double> tmp_3150;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3150 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3149 += tmp_3150;
   result += (1.3333333333333333*MGlu) * tmp_3149;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3151;
   std::complex<double> tmp_3152;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3153;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3153 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3152 += tmp_3153;
   }
   tmp_3151 += tmp_3152;
   result += (-0.5) * tmp_3151;
   std::complex<double> tmp_3154;
   std::complex<double> tmp_3155;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3156;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3156 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3155 += tmp_3156;
   }
   tmp_3154 += tmp_3155;
   result += (-0.5) * tmp_3154;
   std::complex<double> tmp_3157;
   std::complex<double> tmp_3158;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3159;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3159 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3158 += tmp_3159;
   }
   tmp_3157 += tmp_3158;
   result += (-0.5) * tmp_3157;
   std::complex<double> tmp_3160;
   std::complex<double> tmp_3161;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3162;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3162 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3161 += tmp_3162;
   }
   tmp_3160 += tmp_3161;
   result += (-0.5) * tmp_3160;
   std::complex<double> tmp_3163;
   std::complex<double> tmp_3164;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3164 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3163 += tmp_3164;
   result += (-0.6666666666666666) * tmp_3163;
   std::complex<double> tmp_3165;
   std::complex<double> tmp_3166;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3167;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3167 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3166 += tmp_3167;
   }
   tmp_3165 += tmp_3166;
   result += (-0.5) * tmp_3165;
   std::complex<double> tmp_3168;
   std::complex<double> tmp_3169;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3169 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3168 += tmp_3169;
   result += (-1) * tmp_3168;
   std::complex<double> tmp_3170;
   std::complex<double> tmp_3171;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3171 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_3170 += tmp_3171;
   result += (-1.3333333333333333) * tmp_3170;
   std::complex<double> tmp_3172;
   std::complex<double> tmp_3173;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3173 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3172 += tmp_3173;
   result += (-1) * tmp_3172;
   std::complex<double> tmp_3174;
   std::complex<double> tmp_3175;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3175 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3174 += tmp_3175;
   result += (-1) * tmp_3174;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3176;
   std::complex<double> tmp_3177;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3178;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3178 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3177 += tmp_3178;
   }
   tmp_3176 += tmp_3177;
   result += (-0.5) * tmp_3176;
   std::complex<double> tmp_3179;
   std::complex<double> tmp_3180;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3181;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3181 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3180 += tmp_3181;
   }
   tmp_3179 += tmp_3180;
   result += (-0.5) * tmp_3179;
   std::complex<double> tmp_3182;
   std::complex<double> tmp_3183;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3184;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3184 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3183 += tmp_3184;
   }
   tmp_3182 += tmp_3183;
   result += (-0.5) * tmp_3182;
   std::complex<double> tmp_3185;
   std::complex<double> tmp_3186;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3187;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3187 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3186 += tmp_3187;
   }
   tmp_3185 += tmp_3186;
   result += (-0.5) * tmp_3185;
   std::complex<double> tmp_3188;
   std::complex<double> tmp_3189;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3189 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3188 += tmp_3189;
   result += (-0.6666666666666666) * tmp_3188;
   std::complex<double> tmp_3190;
   std::complex<double> tmp_3191;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3192;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3192 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3191 += tmp_3192;
   }
   tmp_3190 += tmp_3191;
   result += (-0.5) * tmp_3190;
   std::complex<double> tmp_3193;
   std::complex<double> tmp_3194;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3194 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3193 += tmp_3194;
   result += (-1) * tmp_3193;
   std::complex<double> tmp_3195;
   std::complex<double> tmp_3196;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3196 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_3195 += tmp_3196;
   result += (-1.3333333333333333) * tmp_3195;
   std::complex<double> tmp_3197;
   std::complex<double> tmp_3198;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3198 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3197 += tmp_3198;
   result += (-1) * tmp_3197;
   std::complex<double> tmp_3199;
   std::complex<double> tmp_3200;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3200 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3199 += tmp_3200;
   result += (-1) * tmp_3199;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3201;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3202;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3202 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpGluconjSdFdPL(gI1,gI2
            ))*CpGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_3201 += tmp_3202;
   }
   result += tmp_3201;
   std::complex<double> tmp_3203;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3204;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3204 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpGluconjSuFuPL(gI1,gI2
            ))*CpGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_3203 += tmp_3204;
   }
   result += tmp_3203;
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpGluVGGluPR())*CpGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_3205;
   std::complex<double> tmp_3206;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3207;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3207 += AbsSqr(CpGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_3206 += tmp_3207;
   }
   tmp_3205 += tmp_3206;
   result += (-0.5) * tmp_3205;
   std::complex<double> tmp_3208;
   std::complex<double> tmp_3209;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3210;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3210 += AbsSqr(CpGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_3209 += tmp_3210;
   }
   tmp_3208 += tmp_3209;
   result += (-0.5) * tmp_3208;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_3211;
   std::complex<double> tmp_3212;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3213;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3213 += AbsSqr(CpGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_3212 += tmp_3213;
   }
   tmp_3211 += tmp_3212;
   result += (-0.5) * tmp_3211;
   std::complex<double> tmp_3214;
   std::complex<double> tmp_3215;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3216;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3216 += AbsSqr(CpGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_3215 += tmp_3216;
   }
   tmp_3214 += tmp_3215;
   result += (-0.5) * tmp_3214;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3217;
   std::complex<double> tmp_3218;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3218 += AbsSqr(CpVZhhAh(gI1,1))*B00(p,MAh(1),Mhh(gI1));
   }
   tmp_3217 += tmp_3218;
   result += (-4) * tmp_3217;
   std::complex<double> tmp_3219;
   std::complex<double> tmp_3220;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3220 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_3219 += tmp_3220;
   result += (0.5) * tmp_3219;
   std::complex<double> tmp_3221;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3222;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3222 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_3222 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_3221 += tmp_3222;
   }
   result += tmp_3221;
   std::complex<double> tmp_3223;
   std::complex<double> tmp_3224;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3225;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3225 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3225 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_3224 += tmp_3225;
   }
   tmp_3223 += tmp_3224;
   result += (0.5) * tmp_3223;
   std::complex<double> tmp_3226;
   std::complex<double> tmp_3227;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3227 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3226 += tmp_3227;
   result += (3) * tmp_3226;
   std::complex<double> tmp_3228;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3228 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3228;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3230 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3229 += tmp_3230;
   result += (3) * tmp_3229;
   std::complex<double> tmp_3231;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3231 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_3231;
   std::complex<double> tmp_3232;
   std::complex<double> tmp_3233;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3234;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3234 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_3233 += tmp_3234;
   }
   tmp_3232 += tmp_3233;
   result += (-12) * tmp_3232;
   std::complex<double> tmp_3235;
   std::complex<double> tmp_3236;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3237;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3237 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_3236 += tmp_3237;
   }
   tmp_3235 += tmp_3236;
   result += (-4) * tmp_3235;
   std::complex<double> tmp_3238;
   std::complex<double> tmp_3239;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3240;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3240 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_3239 += tmp_3240;
   }
   tmp_3238 += tmp_3239;
   result += (-12) * tmp_3238;
   std::complex<double> tmp_3241;
   std::complex<double> tmp_3242;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3243;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3243 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_3242 += tmp_3243;
   }
   tmp_3241 += tmp_3242;
   result += (-4) * tmp_3241;
   std::complex<double> tmp_3244;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3244 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3244;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3245;
   std::complex<double> tmp_3246;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3246 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3245 += tmp_3246;
   result += (0.5) * tmp_3245;
   std::complex<double> tmp_3247;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3248;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3248 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_3248 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_3247 += tmp_3248;
   }
   result += tmp_3247;
   std::complex<double> tmp_3249;
   std::complex<double> tmp_3250;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3250 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3249 += tmp_3250;
   result += (3) * tmp_3249;
   std::complex<double> tmp_3251;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3251 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3251;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3253 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3252 += tmp_3253;
   result += (3) * tmp_3252;
   std::complex<double> tmp_3254;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3254 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3254;
   std::complex<double> tmp_3255;
   std::complex<double> tmp_3256;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3257;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3257 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3256 += tmp_3257;
   }
   tmp_3255 += tmp_3256;
   result += (-12) * tmp_3255;
   std::complex<double> tmp_3258;
   std::complex<double> tmp_3259;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3260;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3260 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3259 += tmp_3260;
   }
   tmp_3258 += tmp_3259;
   result += (-4) * tmp_3258;
   std::complex<double> tmp_3261;
   std::complex<double> tmp_3262;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3262 += AbsSqr(CpconjVWmHpmhh(1,gI2))*B00(p,Mhh(gI2),MHpm(1));
   }
   tmp_3261 += tmp_3262;
   result += (-4) * tmp_3261;
   std::complex<double> tmp_3263;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3263 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3263;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3264;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3265;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3265 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3264 += tmp_3265;
   }
   result += tmp_3264;
   std::complex<double> tmp_3266;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3267;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3267 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3266 += tmp_3267;
   }
   result += tmp_3266;
   std::complex<double> tmp_3268;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3269;
      std::complex<double> tmp_3270;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3270 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3269 += tmp_3270;
      tmp_3268 += (MFe(gI1)) * tmp_3269;
   }
   result += tmp_3268;
   std::complex<double> tmp_3271;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3272;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3272 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3271 += tmp_3272;
   }
   result += tmp_3271;
   std::complex<double> tmp_3273;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3274;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3274 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3273 += tmp_3274;
   }
   result += tmp_3273;
   std::complex<double> tmp_3275;
   std::complex<double> tmp_3276;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3276 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3275 += tmp_3276;
   result += (-4) * tmp_3275;
   std::complex<double> tmp_3277;
   std::complex<double> tmp_3278;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3278 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3277 += tmp_3278;
   result += (-4) * tmp_3277;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3279;
   std::complex<double> tmp_3280;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3281;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3281 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_3280 += tmp_3281;
   }
   tmp_3279 += tmp_3280;
   result += (-0.5) * tmp_3279;
   std::complex<double> tmp_3282;
   std::complex<double> tmp_3283;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3284;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3284 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3283 += tmp_3284;
   }
   tmp_3282 += tmp_3283;
   result += (-0.5) * tmp_3282;
   std::complex<double> tmp_3285;
   std::complex<double> tmp_3286;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3287;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3287 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3286 += tmp_3287;
   }
   tmp_3285 += tmp_3286;
   result += (-0.5) * tmp_3285;
   std::complex<double> tmp_3288;
   std::complex<double> tmp_3289;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3290;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3290 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPR(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3289 += tmp_3290;
   }
   tmp_3288 += tmp_3289;
   result += (-0.5) * tmp_3288;
   std::complex<double> tmp_3291;
   std::complex<double> tmp_3292;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3293;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3293 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3292 += tmp_3293;
   }
   tmp_3291 += tmp_3292;
   result += (-0.5) * tmp_3291;
   std::complex<double> tmp_3294;
   std::complex<double> tmp_3295;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3295 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3294 += tmp_3295;
   result += (-1) * tmp_3294;
   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3297 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3296 += tmp_3297;
   result += (-1) * tmp_3296;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3298;
   std::complex<double> tmp_3299;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3300;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3300 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3299 += tmp_3300;
   }
   tmp_3298 += tmp_3299;
   result += (-0.5) * tmp_3298;
   std::complex<double> tmp_3301;
   std::complex<double> tmp_3302;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3303;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3303 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3302 += tmp_3303;
   }
   tmp_3301 += tmp_3302;
   result += (-0.5) * tmp_3301;
   std::complex<double> tmp_3304;
   std::complex<double> tmp_3305;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3306;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3306 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3305 += tmp_3306;
   }
   tmp_3304 += tmp_3305;
   result += (-0.5) * tmp_3304;
   std::complex<double> tmp_3307;
   std::complex<double> tmp_3308;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3309;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3309 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3308 += tmp_3309;
   }
   tmp_3307 += tmp_3308;
   result += (-0.5) * tmp_3307;
   std::complex<double> tmp_3310;
   std::complex<double> tmp_3311;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3312;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3312 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3311 += tmp_3312;
   }
   tmp_3310 += tmp_3311;
   result += (-0.5) * tmp_3310;
   std::complex<double> tmp_3313;
   std::complex<double> tmp_3314;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3314 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3313 += tmp_3314;
   result += (-1) * tmp_3313;
   std::complex<double> tmp_3315;
   std::complex<double> tmp_3316;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3316 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3315 += tmp_3316;
   result += (-1) * tmp_3315;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3317;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3318;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3318 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3317 += tmp_3318;
   }
   result += tmp_3317;
   std::complex<double> tmp_3319;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3320;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3320 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3319 += tmp_3320;
   }
   result += tmp_3319;
   std::complex<double> tmp_3321;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3322;
      std::complex<double> tmp_3323;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3323 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3322 += tmp_3323;
      tmp_3321 += (MFd(gI1)) * tmp_3322;
   }
   result += tmp_3321;
   std::complex<double> tmp_3324;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3325;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3325 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3324 += tmp_3325;
   }
   result += tmp_3324;
   std::complex<double> tmp_3326;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3327;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3327 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3326 += tmp_3327;
   }
   result += tmp_3326;
   std::complex<double> tmp_3328;
   std::complex<double> tmp_3329;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3329 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3328 += tmp_3329;
   result += (-4) * tmp_3328;
   std::complex<double> tmp_3330;
   std::complex<double> tmp_3331;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3331 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3330 += tmp_3331;
   result += (-4) * tmp_3330;
   std::complex<double> tmp_3332;
   std::complex<double> tmp_3333;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3333 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3332 += tmp_3333;
   result += (1.3333333333333333*MGlu) * tmp_3332;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3334;
   std::complex<double> tmp_3335;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3336;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3336 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3335 += tmp_3336;
   }
   tmp_3334 += tmp_3335;
   result += (-0.5) * tmp_3334;
   std::complex<double> tmp_3337;
   std::complex<double> tmp_3338;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3339;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3339 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3338 += tmp_3339;
   }
   tmp_3337 += tmp_3338;
   result += (-0.5) * tmp_3337;
   std::complex<double> tmp_3340;
   std::complex<double> tmp_3341;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3342;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3342 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3341 += tmp_3342;
   }
   tmp_3340 += tmp_3341;
   result += (-0.5) * tmp_3340;
   std::complex<double> tmp_3343;
   std::complex<double> tmp_3344;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3344 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3343 += tmp_3344;
   result += (-0.6666666666666666) * tmp_3343;
   std::complex<double> tmp_3345;
   std::complex<double> tmp_3346;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3347;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3347 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPR(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3346 += tmp_3347;
   }
   tmp_3345 += tmp_3346;
   result += (-0.5) * tmp_3345;
   std::complex<double> tmp_3348;
   std::complex<double> tmp_3349;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3350;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3350 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3349 += tmp_3350;
   }
   tmp_3348 += tmp_3349;
   result += (-0.5) * tmp_3348;
   std::complex<double> tmp_3351;
   std::complex<double> tmp_3352;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3352 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3351 += tmp_3352;
   result += (-1) * tmp_3351;
   std::complex<double> tmp_3353;
   std::complex<double> tmp_3354;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3354 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3353 += tmp_3354;
   result += (-1) * tmp_3353;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3355;
   std::complex<double> tmp_3356;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3357;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3357 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3356 += tmp_3357;
   }
   tmp_3355 += tmp_3356;
   result += (-0.5) * tmp_3355;
   std::complex<double> tmp_3358;
   std::complex<double> tmp_3359;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3360;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3360 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3359 += tmp_3360;
   }
   tmp_3358 += tmp_3359;
   result += (-0.5) * tmp_3358;
   std::complex<double> tmp_3361;
   std::complex<double> tmp_3362;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3363;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3363 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3362 += tmp_3363;
   }
   tmp_3361 += tmp_3362;
   result += (-0.5) * tmp_3361;
   std::complex<double> tmp_3364;
   std::complex<double> tmp_3365;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3365 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPL(gO1,gI1,1);
   }
   tmp_3364 += tmp_3365;
   result += (-0.6666666666666666) * tmp_3364;
   std::complex<double> tmp_3366;
   std::complex<double> tmp_3367;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3368;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3368 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPL(gO1,gI1,gI2);
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
         tmp_3371 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3370 += tmp_3371;
   }
   tmp_3369 += tmp_3370;
   result += (-0.5) * tmp_3369;
   std::complex<double> tmp_3372;
   std::complex<double> tmp_3373;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3373 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3372 += tmp_3373;
   result += (-1) * tmp_3372;
   std::complex<double> tmp_3374;
   std::complex<double> tmp_3375;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3375 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3374 += tmp_3375;
   result += (-1) * tmp_3374;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3376;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3377;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3377 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3376 += tmp_3377;
   }
   result += tmp_3376;
   std::complex<double> tmp_3378;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3379;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3379 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3378 += tmp_3379;
   }
   result += tmp_3378;
   std::complex<double> tmp_3380;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3381;
      std::complex<double> tmp_3382;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3382 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3381 += tmp_3382;
      tmp_3380 += (MCha(gI1)) * tmp_3381;
   }
   result += tmp_3380;
   std::complex<double> tmp_3383;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3384;
      std::complex<double> tmp_3385;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3385 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3384 += tmp_3385;
      tmp_3383 += (MFu(gI1)) * tmp_3384;
   }
   result += tmp_3383;
   std::complex<double> tmp_3386;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3387;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3387 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3386 += tmp_3387;
   }
   result += tmp_3386;
   std::complex<double> tmp_3388;
   std::complex<double> tmp_3389;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3389 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3388 += tmp_3389;
   result += (-4) * tmp_3388;
   std::complex<double> tmp_3390;
   std::complex<double> tmp_3391;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3391 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3390 += tmp_3391;
   result += (-4) * tmp_3390;
   std::complex<double> tmp_3392;
   std::complex<double> tmp_3393;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3393 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3392 += tmp_3393;
   result += (-4) * tmp_3392;
   std::complex<double> tmp_3394;
   std::complex<double> tmp_3395;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3395 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3394 += tmp_3395;
   result += (1.3333333333333333*MGlu) * tmp_3394;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3396;
   std::complex<double> tmp_3397;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3398;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3398 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3397 += tmp_3398;
   }
   tmp_3396 += tmp_3397;
   result += (-0.5) * tmp_3396;
   std::complex<double> tmp_3399;
   std::complex<double> tmp_3400;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3401;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3401 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3400 += tmp_3401;
   }
   tmp_3399 += tmp_3400;
   result += (-0.5) * tmp_3399;
   std::complex<double> tmp_3402;
   std::complex<double> tmp_3403;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3404;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3404 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPR(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3403 += tmp_3404;
   }
   tmp_3402 += tmp_3403;
   result += (-0.5) * tmp_3402;
   std::complex<double> tmp_3405;
   std::complex<double> tmp_3406;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3407;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3407 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3406 += tmp_3407;
   }
   tmp_3405 += tmp_3406;
   result += (-0.5) * tmp_3405;
   std::complex<double> tmp_3408;
   std::complex<double> tmp_3409;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3409 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3408 += tmp_3409;
   result += (-0.6666666666666666) * tmp_3408;
   std::complex<double> tmp_3410;
   std::complex<double> tmp_3411;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3412;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3412 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3411 += tmp_3412;
   }
   tmp_3410 += tmp_3411;
   result += (-0.5) * tmp_3410;
   std::complex<double> tmp_3413;
   std::complex<double> tmp_3414;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3414 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3413 += tmp_3414;
   result += (-1) * tmp_3413;
   std::complex<double> tmp_3415;
   std::complex<double> tmp_3416;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3416 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3415 += tmp_3416;
   result += (-1) * tmp_3415;
   std::complex<double> tmp_3417;
   std::complex<double> tmp_3418;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3418 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3417 += tmp_3418;
   result += (-1) * tmp_3417;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3419;
   std::complex<double> tmp_3420;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3421;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3421 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3420 += tmp_3421;
   }
   tmp_3419 += tmp_3420;
   result += (-0.5) * tmp_3419;
   std::complex<double> tmp_3422;
   std::complex<double> tmp_3423;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3424;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3424 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3423 += tmp_3424;
   }
   tmp_3422 += tmp_3423;
   result += (-0.5) * tmp_3422;
   std::complex<double> tmp_3425;
   std::complex<double> tmp_3426;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3427;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3427 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3426 += tmp_3427;
   }
   tmp_3425 += tmp_3426;
   result += (-0.5) * tmp_3425;
   std::complex<double> tmp_3428;
   std::complex<double> tmp_3429;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3430;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3430 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3429 += tmp_3430;
   }
   tmp_3428 += tmp_3429;
   result += (-0.5) * tmp_3428;
   std::complex<double> tmp_3431;
   std::complex<double> tmp_3432;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3432 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPL(gO1,gI1,1);
   }
   tmp_3431 += tmp_3432;
   result += (-0.6666666666666666) * tmp_3431;
   std::complex<double> tmp_3433;
   std::complex<double> tmp_3434;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3435;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3435 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3434 += tmp_3435;
   }
   tmp_3433 += tmp_3434;
   result += (-0.5) * tmp_3433;
   std::complex<double> tmp_3436;
   std::complex<double> tmp_3437;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3437 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3436 += tmp_3437;
   result += (-1) * tmp_3436;
   std::complex<double> tmp_3438;
   std::complex<double> tmp_3439;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3439 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3438 += tmp_3439;
   result += (-1) * tmp_3438;
   std::complex<double> tmp_3440;
   std::complex<double> tmp_3441;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3441 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3440 += tmp_3441;
   result += (-1) * tmp_3440;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3442;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3443;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3443 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3442 += tmp_3443;
   }
   result += tmp_3442;
   std::complex<double> tmp_3444;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3445;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3445 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3444 += tmp_3445;
   }
   result += tmp_3444;
   std::complex<double> tmp_3446;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3447;
      std::complex<double> tmp_3448;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3448 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3447 += tmp_3448;
      tmp_3446 += (MCha(gI1)) * tmp_3447;
   }
   result += tmp_3446;
   std::complex<double> tmp_3449;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3450;
      std::complex<double> tmp_3451;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3451 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3450 += tmp_3451;
      tmp_3449 += (MFu(gI1)) * tmp_3450;
   }
   result += tmp_3449;
   std::complex<double> tmp_3452;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3453;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3453 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3452 += tmp_3453;
   }
   result += tmp_3452;
   std::complex<double> tmp_3454;
   std::complex<double> tmp_3455;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3455 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3454 += tmp_3455;
   result += (-4) * tmp_3454;
   std::complex<double> tmp_3456;
   std::complex<double> tmp_3457;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3457 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3456 += tmp_3457;
   result += (-4) * tmp_3456;
   std::complex<double> tmp_3458;
   std::complex<double> tmp_3459;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3459 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3458 += tmp_3459;
   result += (-4) * tmp_3458;
   std::complex<double> tmp_3460;
   std::complex<double> tmp_3461;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3461 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3460 += tmp_3461;
   result += (1.3333333333333333*MGlu) * tmp_3460;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3462;
   std::complex<double> tmp_3463;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3464;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3464 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3463 += tmp_3464;
   }
   tmp_3462 += tmp_3463;
   result += (-0.5) * tmp_3462;
   std::complex<double> tmp_3465;
   std::complex<double> tmp_3466;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3467;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3467 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3466 += tmp_3467;
   }
   tmp_3465 += tmp_3466;
   result += (-0.5) * tmp_3465;
   std::complex<double> tmp_3468;
   std::complex<double> tmp_3469;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3470;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3470 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3469 += tmp_3470;
   }
   tmp_3468 += tmp_3469;
   result += (-0.5) * tmp_3468;
   std::complex<double> tmp_3471;
   std::complex<double> tmp_3472;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3473;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3473 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3472 += tmp_3473;
   }
   tmp_3471 += tmp_3472;
   result += (-0.5) * tmp_3471;
   std::complex<double> tmp_3474;
   std::complex<double> tmp_3475;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3475 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3474 += tmp_3475;
   result += (-0.6666666666666666) * tmp_3474;
   std::complex<double> tmp_3476;
   std::complex<double> tmp_3477;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3478;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3478 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3477 += tmp_3478;
   }
   tmp_3476 += tmp_3477;
   result += (-0.5) * tmp_3476;
   std::complex<double> tmp_3479;
   std::complex<double> tmp_3480;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3480 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3479 += tmp_3480;
   result += (-1) * tmp_3479;
   std::complex<double> tmp_3481;
   std::complex<double> tmp_3482;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3482 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3481 += tmp_3482;
   result += (-1) * tmp_3481;
   std::complex<double> tmp_3483;
   std::complex<double> tmp_3484;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3484 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3483 += tmp_3484;
   result += (-1) * tmp_3483;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3485;
   std::complex<double> tmp_3486;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3487;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3487 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3486 += tmp_3487;
   }
   tmp_3485 += tmp_3486;
   result += (-0.5) * tmp_3485;
   std::complex<double> tmp_3488;
   std::complex<double> tmp_3489;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3490;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3490 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3489 += tmp_3490;
   }
   tmp_3488 += tmp_3489;
   result += (-0.5) * tmp_3488;
   std::complex<double> tmp_3491;
   std::complex<double> tmp_3492;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3493;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3493 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3492 += tmp_3493;
   }
   tmp_3491 += tmp_3492;
   result += (-0.5) * tmp_3491;
   std::complex<double> tmp_3494;
   std::complex<double> tmp_3495;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3496;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3496 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3495 += tmp_3496;
   }
   tmp_3494 += tmp_3495;
   result += (-0.5) * tmp_3494;
   std::complex<double> tmp_3497;
   std::complex<double> tmp_3498;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3498 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3497 += tmp_3498;
   result += (-0.6666666666666666) * tmp_3497;
   std::complex<double> tmp_3499;
   std::complex<double> tmp_3500;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3501;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3501 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3500 += tmp_3501;
   }
   tmp_3499 += tmp_3500;
   result += (-0.5) * tmp_3499;
   std::complex<double> tmp_3502;
   std::complex<double> tmp_3503;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3503 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3502 += tmp_3503;
   result += (-1) * tmp_3502;
   std::complex<double> tmp_3504;
   std::complex<double> tmp_3505;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3505 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3504 += tmp_3505;
   result += (-1) * tmp_3504;
   std::complex<double> tmp_3506;
   std::complex<double> tmp_3507;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3507 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3506 += tmp_3507;
   result += (-1) * tmp_3506;

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
   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3509 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3508 += tmp_3509;
   result += (-0.5) * tmp_3508;
   std::complex<double> tmp_3510;
   std::complex<double> tmp_3511;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3511 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3510 += tmp_3511;
   result += (-1) * tmp_3510;
   std::complex<double> tmp_3512;
   std::complex<double> tmp_3513;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3513 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3512 += tmp_3513;
   result += (-0.5) * tmp_3512;
   std::complex<double> tmp_3514;
   std::complex<double> tmp_3515;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3515 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_3514 += tmp_3515;
   result += (2) * tmp_3514;
   std::complex<double> tmp_3516;
   std::complex<double> tmp_3517;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3517 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3516 += tmp_3517;
   result += (6) * tmp_3516;
   std::complex<double> tmp_3518;
   std::complex<double> tmp_3519;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3519 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3518 += tmp_3519;
   result += (2) * tmp_3518;
   std::complex<double> tmp_3520;
   std::complex<double> tmp_3521;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3521 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3520 += tmp_3521;
   result += (6) * tmp_3520;
   std::complex<double> tmp_3522;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3522 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_3522;
   std::complex<double> tmp_3523;
   std::complex<double> tmp_3524;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3524 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3523 += tmp_3524;
   result += (-3) * tmp_3523;
   std::complex<double> tmp_3525;
   std::complex<double> tmp_3526;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3526 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3525 += tmp_3526;
   result += (-1) * tmp_3525;
   std::complex<double> tmp_3527;
   std::complex<double> tmp_3528;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3528 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3527 += tmp_3528;
   result += (-3) * tmp_3527;
   std::complex<double> tmp_3529;
   std::complex<double> tmp_3530;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3530 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3529 += tmp_3530;
   result += (-1) * tmp_3529;
   std::complex<double> tmp_3531;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3531 += A0(MFv(gI1))*(CpUhhFvFvPL(gO1,gI1,gI1) + CpUhhFvFvPR(gO1,
         gI1,gI1))*MFv(gI1);
   }
   result += tmp_3531;

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
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZD;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZD,
            eigenvalue_error);
         problems.flag_bad_mass(MSSMRHN_info::Sd, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZD);
      #endif

      PHYSICAL(MSd(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZD) = mix_ZD;
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
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZU;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZU,
            eigenvalue_error);
         problems.flag_bad_mass(MSSMRHN_info::Su, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZU);
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
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZE;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZE,
            eigenvalue_error);
         problems.flag_bad_mass(MSSMRHN_info::Se, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZE);
      #endif

      PHYSICAL(MSe(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZE) = mix_ZE;
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
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(MSSMRHN_info::Sv, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV);
      #endif

      PHYSICAL(MSv(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZV) = mix_ZV;
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
            problems.flag_bad_mass(MSSMRHN_info::hh,
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
      problems.flag_no_pole_mass_convergence(MSSMRHN_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(MSSMRHN_info::hh);
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
      Eigen::Matrix<double,2,2> self_energy;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Ah());

      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MAh(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = i1; i2 < 2; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,2,2> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(MSSMRHN_info::Ah,
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
      problems.flag_no_pole_mass_convergence(MSSMRHN_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(MSSMRHN_info::Ah);
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
         const Eigen::Matrix<double,2,2> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(MSSMRHN_info::Hpm,
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
      problems.flag_no_pole_mass_convergence(MSSMRHN_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(MSSMRHN_info::Hpm);
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
         problems.flag_bad_mass(MSSMRHN_info::Chi, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN);
      #endif
      if (es == 0)
         PHYSICAL(ZN) = mix_ZN;
      PHYSICAL(MChi(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,6,6> self_energy_1;
   Eigen::Matrix<double,6,6> self_energy_PL;
   Eigen::Matrix<double,6,6> self_energy_PR;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Fv());
   for (unsigned es = 0; es < 6; ++es) {
      const double p = Abs(MFv(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = 0; i2 < 6; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fv_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fv_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fv_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,6,6> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,6,6> M_loop(M_tree + 0.5 * (delta_M +
         delta_M.transpose()));
      Eigen::Array<double,6,1> eigen_values;
      decltype(UV) mix_UV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_UV,
            eigenvalue_error);
         problems.flag_bad_mass(MSSMRHN_info::Fv, eigenvalue_error
            > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_UV);
      #endif
      if (es == 0)
         PHYSICAL(UV) = mix_UV;
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
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(MSSMRHN_info::Cha, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
      decltype(ZEL) mix_ZEL;
      decltype(ZER) mix_ZER;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER, eigenvalue_error)
         ;
      problems.flag_bad_mass(MSSMRHN_info::Fe, eigenvalue_error >
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
      decltype(ZDL) mix_ZDL;
      decltype(ZDR) mix_ZDR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR, eigenvalue_error)
         ;
      problems.flag_bad_mass(MSSMRHN_info::Fd, eigenvalue_error >
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
      qcd_2l = -0.005191204615668296*Power(g3,4) -
         0.0032883224409535764*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(2))) -
         0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(currentScale,2)/Sqr(
         MFu(2))));
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
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error)
         ;
      problems.flag_bad_mass(MSSMRHN_info::Fu, eigenvalue_error >
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


double CLASSNAME::v() const
{
   return Sqrt(Sqr(vd) + Sqr(vu));
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


std::ostream& operator<<(std::ostream& ostr, const MSSMRHN_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
