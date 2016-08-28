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

// File generated at Sun 28 Aug 2016 15:02:45

/**
 * @file SplitMSSM_mass_eigenstates.cpp
 * @brief implementation of the SplitMSSM model class
 *
 * Contains the definition of the SplitMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sun 28 Aug 2016 15:02:45 with FlexibleSUSY
 * 1.6.1 (git commit: de78c399bea838d002fe54cbc91eb6f88fa62626) and SARAH 4.9.1 .
 */

#include "SplitMSSM_mass_eigenstates.hpp"
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

#include "sm_twoloophiggs.hpp"

#include "split_threeloophiggs.hpp"


#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace SplitMSSM_info;

#define CLASSNAME SplitMSSM_mass_eigenstates

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

CLASSNAME::SplitMSSM_mass_eigenstates(const SplitMSSM_input_parameters& input_)
   : SplitMSSM_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(SplitMSSM_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MHp(0), MFv(Eigen::Array<double,3,1>::Zero()), MGlu(0), MAh(0),
      Mhh(0), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,3,1>
      ::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MChi(Eigen::Array<double,4
      ,1>::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MVWp(0), MVP(0), MVZ(0
      )

   , Vd(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ud(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Vu(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Uu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ve(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ue(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZN(Eigen::Matrix<std::complex<double>,4,
      4>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZZ(Eigen::Matrix<double,2,
      2>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~SplitMSSM_mass_eigenstates()
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

const SplitMSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

SplitMSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const SplitMSSM_physical& physical_)
{
   physical = physical_;
}

const Problems<SplitMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<SplitMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
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

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh());

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
   SplitMSSM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mu2(gsl_vector_get(x, 0));


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
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision)),
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

   mu2 = solver->get_solution(0);


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

   const double old_mu2 = mu2;

   mu2 = Re(Lambdax*Sqr(v));

   const bool is_finite = IsFinite(mu2);

   if (!is_finite) {
      mu2 = old_mu2;
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

   x_init[0] = mu2;


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
      tadpole[0] += Re(tadpole_hh());

      if (ewsb_loop_order > 1) {

      }
   }

   double mu2;

   mu2 = Re((0.7071067811865475*(1.4142135623730951*Power(v,3)*Lambdax -
      tadpole[0]))/v);

   const bool is_finite = IsFinite(mu2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mu2;


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
   SplitMSSM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mu2 = gsl_vector_get(x, 0);

   model->set_mu2(mu2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = mu2;


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
           "SplitMSSM\n"
           "========================================\n";
   SplitMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
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
   const auto old_mu2 = mu2;

   solve_ewsb_tree_level();

   calculate_MVPVZ();
   calculate_MVWp();
   calculate_MCha();
   calculate_MChi();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_Mhh();
   calculate_MAh();
   calculate_MGlu();
   calculate_MFv();
   calculate_MHp();
   calculate_MVG();

   mu2 = old_mu2;

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
   typedef void (CLASSNAME::*Mem_fun_t)();
   typedef CLASSNAME* Obj_ptr_t;

   auto fut_MCha = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MCha_pole, this);
   auto fut_MChi = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MChi_pole, this);
   auto fut_MGlu = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MGlu_pole, this);

   if (calculate_sm_pole_masses) {
      auto fut_MVG = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MVG_pole, this);
      auto fut_MFv = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MFv_pole, this);
      auto fut_Mhh = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_Mhh_pole, this);
      auto fut_MVP = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MVP_pole, this);
      auto fut_MVZ = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MVZ_pole, this);
      auto fut_MFd = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MFd_pole, this);
      auto fut_MFu = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MFu_pole, this);
      auto fut_MFe = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MFe_pole, this);
      auto fut_MVWp = run_async<Mem_fun_t, Obj_ptr_t>(&CLASSNAME::calculate_MVWp_pole, this);
      fut_MVG.get();
      fut_MFv.get();
      fut_Mhh.get();
      fut_MVP.get();
      fut_MVZ.get();
      fut_MFd.get();
      fut_MFu.get();
      fut_MFe.get();
      fut_MVWp.get();
   }

   fut_MCha.get();
   fut_MChi.get();
   fut_MGlu.get();

#else
   calculate_MCha_pole();
   calculate_MChi_pole();
   calculate_MGlu_pole();

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_Mhh_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MVWp_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHp) = MHp;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
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
   PHYSICAL(MVWp) = MVWp;
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

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh) < 0.) problems.flag_tachyon(hh);

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
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MGlu = 0.;
   MAh = 0.;
   Mhh = 0.;
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
   MVWp = 0.;
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
   SplitMSSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MHp = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MGlu = pars(5);
   MAh = pars(6);
   Mhh = pars(7);
   MFd(0) = pars(8);
   MFd(1) = pars(9);
   MFd(2) = pars(10);
   MFu(0) = pars(11);
   MFu(1) = pars(12);
   MFu(2) = pars(13);
   MFe(0) = pars(14);
   MFe(1) = pars(15);
   MFe(2) = pars(16);
   MChi(0) = pars(17);
   MChi(1) = pars(18);
   MChi(2) = pars(19);
   MChi(3) = pars(20);
   MCha(0) = pars(21);
   MCha(1) = pars(22);
   MVWp = pars(23);
   MVP = pars(24);
   MVZ = pars(25);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(26);

   pars(0) = MVG;
   pars(1) = MHp;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MGlu;
   pars(6) = MAh;
   pars(7) = Mhh;
   pars(8) = MFd(0);
   pars(9) = MFd(1);
   pars(10) = MFd(2);
   pars(11) = MFu(0);
   pars(12) = MFu(1);
   pars(13) = MFu(2);
   pars(14) = MFe(0);
   pars(15) = MFe(1);
   pars(16) = MFe(2);
   pars(17) = MChi(0);
   pars(18) = MChi(1);
   pars(19) = MChi(2);
   pars(20) = MChi(3);
   pars(21) = MCha(0);
   pars(22) = MCha(1);
   pars(23) = MVWp;
   pars(24) = MVP;
   pars(25) = MVZ;

   return pars;
}

std::string CLASSNAME::name() const
{
   return "SplitMSSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   SplitMSSM_soft_parameters::run_to(scale, eps);
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

double CLASSNAME::get_mass_matrix_Hp() const
{
   const double mass_matrix_Hp = Re(-mu2 + Lambdax*Sqr(v) + 0.5*Sqr(g2)*
      Sqr(v));

   return mass_matrix_Hp;
}

void CLASSNAME::calculate_MHp()
{
   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = calculate_singlet_mass(mass_matrix_Hp);

   if (MHp < 0.) {
      problems.flag_tachyon(SplitMSSM_info::Hp);
   }

   MHp = AbsSqrt(MHp);
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

double CLASSNAME::get_mass_matrix_Ah() const
{
   const double mass_matrix_Ah = Re(-mu2 + Lambdax*Sqr(v) + 0.5*Sqr(v)*
      Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW())));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = calculate_singlet_mass(mass_matrix_Ah);

   if (MAh < 0.) {
      problems.flag_tachyon(SplitMSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double CLASSNAME::get_mass_matrix_hh() const
{
   const double mass_matrix_hh = Re(-mu2 + 3*Lambdax*Sqr(v));

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = calculate_singlet_mass(mass_matrix_hh);

   if (Mhh < 0.) {
      problems.flag_tachyon(SplitMSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = v*Yd(0,0);
   mass_matrix_Fd(0,1) = v*Yd(1,0);
   mass_matrix_Fd(0,2) = v*Yd(2,0);
   mass_matrix_Fd(1,0) = v*Yd(0,1);
   mass_matrix_Fd(1,1) = v*Yd(1,1);
   mass_matrix_Fd(1,2) = v*Yd(2,1);
   mass_matrix_Fd(2,0) = v*Yd(0,2);
   mass_matrix_Fd(2,1) = v*Yd(1,2);
   mass_matrix_Fd(2,2) = v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(SplitMSSM_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = v*Yu(0,0);
   mass_matrix_Fu(0,1) = v*Yu(1,0);
   mass_matrix_Fu(0,2) = v*Yu(2,0);
   mass_matrix_Fu(1,0) = v*Yu(0,1);
   mass_matrix_Fu(1,1) = v*Yu(1,1);
   mass_matrix_Fu(1,2) = v*Yu(2,1);
   mass_matrix_Fu(2,0) = v*Yu(0,2);
   mass_matrix_Fu(2,1) = v*Yu(1,2);
   mass_matrix_Fu(2,2) = v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(SplitMSSM_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = v*Ye(0,0);
   mass_matrix_Fe(0,1) = v*Ye(1,0);
   mass_matrix_Fe(0,2) = v*Ye(2,0);
   mass_matrix_Fe(1,0) = v*Ye(0,1);
   mass_matrix_Fe(1,1) = v*Ye(1,1);
   mass_matrix_Fe(1,2) = v*Ye(2,1);
   mass_matrix_Fe(2,0) = v*Ye(0,2);
   mass_matrix_Fe(2,1) = v*Ye(1,2);
   mass_matrix_Fe(2,2) = v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(SplitMSSM_info::Fe, eigenvalue_error >
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
   mass_matrix_Chi(0,2) = -0.7071067811865475*gYd*v;
   mass_matrix_Chi(0,3) = 0.7071067811865475*gYu*v;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.7071067811865475*g2d*v;
   mass_matrix_Chi(1,3) = -0.7071067811865475*g2u*v;
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
   problems.flag_bad_mass(SplitMSSM_info::Chi, eigenvalue_error >
      precision * Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = g2u*v;
   mass_matrix_Cha(1,0) = g2d*v;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(SplitMSSM_info::Cha, eigenvalue_error >
      precision * Abs(MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

double CLASSNAME::get_mass_matrix_VWp() const
{
   const double mass_matrix_VWp = Re(0.5*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void CLASSNAME::calculate_MVWp()
{
   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = calculate_singlet_mass(mass_matrix_VWp);

   if (MVWp < 0.) {
      problems.flag_tachyon(SplitMSSM_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.3*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.3872983346207417*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.5*Sqr(g2)*Sqr(v);

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
   double result = Re(-1.4142135623730951*mu2*v + 1.4142135623730951*Power(v,3)
      *Lambdax);

   return result;
}



double CLASSNAME::CpconjHpHphh() const
{
   double result = 0.0;

   result = -1.4142135623730951*v*Lambdax;

   return result;
}

double CLASSNAME::CpconjHpVWpVP() const
{
   double result = 0.0;

   result = 0.5477225575051661*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjHpVZVWp() const
{
   double result = 0.0;

   result = -0.5477225575051661*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpHpgWpCbargZ() const
{
   double result = 0.0;

   result = 0.35355339059327373*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpbargWpCgZ() const
{
   double result = 0.0;

   result = 0.07071067811865475*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpgZbargWp() const
{
   double result = 0.0;

   result = 0.07071067811865475*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpbargZgWp() const
{
   double result = 0.0;

   result = 0.35355339059327373*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjHpAhAh() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

double CLASSNAME::CpHpconjHphhhh() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

double CLASSNAME::CpHpconjHpconjHpHp() const
{
   double result = 0.0;

   result = -2*Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpVWpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjHpVWphh() const
{
   double result = 0.0;

   result = -0.5*g2;

   return result;
}

double CLASSNAME::CpconjHpVPHp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpVZHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjHpconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpconjHpVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarChaChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*UM(gI1,1)*(gYd*ZN(gI2,0) + g2d*ZN(gI2,1)) + g2d
      *UM(gI1,0)*ZN(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarChaChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*Conj(UP(gI1,1))*(gYu*Conj(ZN(gI2,0)) + g2u*Conj
      (ZN(gI2,1))) - g2u*Conj(UP(gI1,0))*Conj(ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1;
      std::complex<double> tmp_2;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_1 += tmp_2;
      tmp_0 += (Vd(gI1,j2)) * tmp_1;
   }
   result += tmp_0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3;
   std::complex<double> tmp_4;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5;
      std::complex<double> tmp_6;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_6 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_5 += tmp_6;
      tmp_4 += (Conj(Vu(gI2,j2))) * tmp_5;
   }
   tmp_3 += tmp_4;
   result += (-1) * tmp_3;

   return result;
}

double CLASSNAME::CpconjHpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_7;
   std::complex<double> tmp_8;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_8 += Conj(Ue(gI1,j1))*Ye(j1,gI2);
   }
   tmp_7 += tmp_8;
   result += (-1) * tmp_7;

   return result;
}

double CLASSNAME::CpAhhhAh() const
{
   double result = 0.0;

   result = -1.4142135623730951*v*Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpgWp() const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.35355339059327373)*v*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpCgWpC() const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.35355339059327373)*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpAhAhAhAh() const
{
   double result = 0.0;

   result = -3*Lambdax;

   return result;
}

double CLASSNAME::CpAhAhhhhh() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

double CLASSNAME::CpAhAhconjHpHp() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhVZhh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjVWpHp() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpAhAhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpAhbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2d*UM(gI1,1)*UP(gI2,0
      ) - g2u*UM(gI1,0)*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2d*Conj(UM(gI2,1))*
      Conj(UP(gI1,0)) - g2u*Conj(UM(gI2,0))*Conj(UP(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_9;
   std::complex<double> tmp_10;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_11;
      std::complex<double> tmp_12;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_12 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_11 += tmp_12;
      tmp_10 += (Vd(gI1,j2)) * tmp_11;
   }
   tmp_9 += tmp_10;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_9;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_13;
   std::complex<double> tmp_14;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_15;
      std::complex<double> tmp_16;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_16 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_15 += tmp_16;
      tmp_14 += (Conj(Vd(gI2,j2))) * tmp_15;
   }
   tmp_13 += tmp_14;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_13;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_17;
   std::complex<double> tmp_18;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_19;
      std::complex<double> tmp_20;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_20 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_19 += tmp_20;
      tmp_18 += (Ve(gI1,j2)) * tmp_19;
   }
   tmp_17 += tmp_18;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_17;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_21;
   std::complex<double> tmp_22;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_23;
      std::complex<double> tmp_24;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_24 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_23 += tmp_24;
      tmp_22 += (Conj(Ve(gI2,j2))) * tmp_23;
   }
   tmp_21 += tmp_22;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_21;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_25;
   std::complex<double> tmp_26;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_27;
      std::complex<double> tmp_28;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_28 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_27 += tmp_28;
      tmp_26 += (Vu(gI1,j2)) * tmp_27;
   }
   tmp_25 += tmp_26;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_25;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_29;
   std::complex<double> tmp_30;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_31;
      std::complex<double> tmp_32;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_32 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_31 += tmp_32;
      tmp_30 += (Conj(Vu(gI2,j2))) * tmp_31;
   }
   tmp_29 += tmp_30;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_29;

   return result;
}

std::complex<double> CLASSNAME::CpAhChiChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(ZN(gI1,2)*(-(gYd*ZN(gI2,0)) + g2d*ZN(
      gI2,1)) + ZN(gI1,3)*(-(gYu*ZN(gI2,0)) + g2u*ZN(gI2,1)) - gYd*ZN(gI1,0)*ZN(
      gI2,2) + g2d*ZN(gI1,1)*ZN(gI2,2) - gYu*ZN(gI1,0)*ZN(gI2,3) + g2u*ZN(gI1,1)*
      ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpAhChiChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(ZN(gI1,2))*(gYd*Conj(ZN(gI2,0)) -
      g2d*Conj(ZN(gI2,1))) + Conj(ZN(gI1,3))*(gYu*Conj(ZN(gI2,0)) - g2u*Conj(ZN(
      gI2,1))) + gYd*Conj(ZN(gI1,0))*Conj(ZN(gI2,2)) - g2d*Conj(ZN(gI1,1))*Conj(ZN
      (gI2,2)) + gYu*Conj(ZN(gI1,0))*Conj(ZN(gI2,3)) - g2u*Conj(ZN(gI1,1))*Conj(ZN
      (gI2,3)));

   return result;
}

double CLASSNAME::CphhAhAh() const
{
   double result = 0.0;

   result = -1.4142135623730951*v*Lambdax;

   return result;
}

double CLASSNAME::Cphhhhhh() const
{
   double result = 0.0;

   result = -4.242640687119286*v*Lambdax;

   return result;
}

double CLASSNAME::CphhVZVZ() const
{
   double result = 0.0;

   result = 0.7071067811865475*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhbargWpgWp() const
{
   double result = 0.0;

   result = -0.35355339059327373*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargWpCgWpC() const
{
   double result = 0.0;

   result = -0.35355339059327373*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargZgZ() const
{
   double result = 0.0;

   result = -0.35355339059327373*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1
      *Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhconjHpHp() const
{
   double result = 0.0;

   result = -1.4142135623730951*v*Lambdax;

   return result;
}

double CLASSNAME::CphhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.7071067811865475*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhhhAhAh() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

double CLASSNAME::Cphhhhhhhh() const
{
   double result = 0.0;

   result = -3*Lambdax;

   return result;
}

double CLASSNAME::CphhhhconjHpHp() const
{
   double result = 0.0;

   result = -Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CphhVZAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhconjVWpHp() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

double CLASSNAME::CphhhhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CphhbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2d*UM(gI1,1)*UP(gI2,0) + g2u*UM(gI1,0)*UP(gI2
      ,1));

   return result;
}

std::complex<double> CLASSNAME::CphhbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2d*Conj(UM(gI2,1))*Conj(UP(gI1,0)) + g2u*Conj
      (UM(gI2,0))*Conj(UP(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_33;
   std::complex<double> tmp_34;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_35;
      std::complex<double> tmp_36;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_36 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_35 += tmp_36;
      tmp_34 += (Vd(gI1,j2)) * tmp_35;
   }
   tmp_33 += tmp_34;
   result += (-0.7071067811865475) * tmp_33;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_37;
   std::complex<double> tmp_38;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_39;
      std::complex<double> tmp_40;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_40 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_39 += tmp_40;
      tmp_38 += (Conj(Vd(gI2,j2))) * tmp_39;
   }
   tmp_37 += tmp_38;
   result += (-0.7071067811865475) * tmp_37;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_41;
   std::complex<double> tmp_42;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_43;
      std::complex<double> tmp_44;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_44 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_43 += tmp_44;
      tmp_42 += (Ve(gI1,j2)) * tmp_43;
   }
   tmp_41 += tmp_42;
   result += (-0.7071067811865475) * tmp_41;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_45;
   std::complex<double> tmp_46;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_47;
      std::complex<double> tmp_48;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_48 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_47 += tmp_48;
      tmp_46 += (Conj(Ve(gI2,j2))) * tmp_47;
   }
   tmp_45 += tmp_46;
   result += (-0.7071067811865475) * tmp_45;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_49;
   std::complex<double> tmp_50;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_51;
      std::complex<double> tmp_52;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_52 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_51 += tmp_52;
      tmp_50 += (Vu(gI1,j2)) * tmp_51;
   }
   tmp_49 += tmp_50;
   result += (-0.7071067811865475) * tmp_49;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_53;
   std::complex<double> tmp_54;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_55;
      std::complex<double> tmp_56;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_56 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_55 += tmp_56;
      tmp_54 += (Conj(Vu(gI2,j2))) * tmp_55;
   }
   tmp_53 += tmp_54;
   result += (-0.7071067811865475) * tmp_53;

   return result;
}

std::complex<double> CLASSNAME::CphhChiChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(ZN(gI1,2)*(gYd*ZN(gI2,0) - g2d*ZN(gI2,1)) + ZN(gI1,3)*(-(gYu*
      ZN(gI2,0)) + g2u*ZN(gI2,1)) + gYd*ZN(gI1,0)*ZN(gI2,2) - g2d*ZN(gI1,1)*ZN(gI2
      ,2) - gYu*ZN(gI1,0)*ZN(gI2,3) + g2u*ZN(gI1,1)*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CphhChiChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI1,2))*(gYd*Conj(ZN(gI2,0)) - g2d*Conj(ZN(gI2,1))) +
      Conj(ZN(gI1,3))*(-(gYu*Conj(ZN(gI2,0))) + g2u*Conj(ZN(gI2,1))) + gYd*Conj(ZN
      (gI1,0))*Conj(ZN(gI2,2)) - g2d*Conj(ZN(gI1,1))*Conj(ZN(gI2,2)) - gYu*Conj(ZN
      (gI1,0))*Conj(ZN(gI2,3)) + g2u*Conj(ZN(gI1,1))*Conj(ZN(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZhh() const
{
   double result = 0.0;

   result = 0.7071067811865475*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWpgWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbargWpCgWpC() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjHpHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjVWpHp() const
{
   double result = 0.0;

   result = -0.5477225575051661*g1*g2*v*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHpHp() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpVZconjVWpVWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

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

double CLASSNAME::CpVZVZconjVWpVWp1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpHpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjVWpHphh() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

double CLASSNAME::CpconjVWpVPHp() const
{
   double result = 0.0;

   result = 0.5477225575051661*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVWphh() const
{
   double result = 0.0;

   result = 0.7071067811865475*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVZHp() const
{
   double result = 0.0;

   result = -0.5477225575051661*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargPgWp() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgZ() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargZgWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpAhAh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWphhhh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjHpHp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVWpVP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVZVWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarChaChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN(gI2,1))*UM(gI1,0) + 1.4142135623730951*Conj(ZN(
      gI2,2))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarChaChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UP(gI1,0))*ZN(gI2,1)) + 0.7071067811865475*g2*Conj(UP(gI1
      ,1))*ZN(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_57;
   std::complex<double> tmp_58;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_58 += Conj(Vu(gI2,j1))*Vd(gI1,j1);
   }
   tmp_57 += tmp_58;
   result += (-0.7071067811865475*g2) * tmp_57;

   return result;
}

double CLASSNAME::CpconjVWpbarFdFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gI1,gI2);
   }

   return result;
}

double CLASSNAME::CpconjVWpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp1() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp2() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp3() const
{
   double result = 0.0;

   result = 2*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_59;
      std::complex<double> tmp_60;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_60 += Conj(Vd(gI1,j2))*Yd(gO2,j2);
      }
      tmp_59 += tmp_60;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_59;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_61;
      std::complex<double> tmp_62;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_62 += Conj(Yd(j1,gO1))*Ud(gI1,j1);
      }
      tmp_61 += tmp_62;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_61;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_63;
      std::complex<double> tmp_64;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_64 += Conj(Vd(gI2,j2))*Yd(gO2,j2);
      }
      tmp_63 += tmp_64;
      result += (-0.7071067811865475) * tmp_63;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_65;
      std::complex<double> tmp_66;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_66 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_65 += tmp_66;
      result += (-0.7071067811865475) * tmp_65;
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

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_67;
      std::complex<double> tmp_68;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_68 += Conj(Vu(gI2,j2))*Yd(gO2,j2);
      }
      tmp_67 += tmp_68;
      result += (-1) * tmp_67;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_69;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_69 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      result += tmp_69;
   }

   return result;
}

double CLASSNAME::CpbarUFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vu(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_71 += Conj(Vu(gI1,j2))*Yu(gO2,j2);
      }
      tmp_70 += tmp_71;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_70;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_72;
      std::complex<double> tmp_73;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_73 += Conj(Yu(j1,gO1))*Uu(gI1,j1);
      }
      tmp_72 += tmp_73;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_72;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_74;
      std::complex<double> tmp_75;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_75 += Conj(Vu(gI2,j2))*Yu(gO2,j2);
      }
      tmp_74 += tmp_75;
      result += (-0.7071067811865475) * tmp_74;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_76;
      std::complex<double> tmp_77;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_77 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_76 += tmp_77;
      result += (-0.7071067811865475) * tmp_76;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_78;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_78 += Conj(Vd(gI2,j2))*Yu(gO2,j2);
      }
      result += tmp_78;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_79;
      std::complex<double> tmp_80;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_80 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_79 += tmp_80;
      result += (-1) * tmp_79;
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

double CLASSNAME::CpbarUFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vd(gI2,gO1));
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

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_81;
      std::complex<double> tmp_82;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_82 += Conj(Ve(gI1,j2))*Ye(gO2,j2);
      }
      tmp_81 += tmp_82;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_81;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_83;
      std::complex<double> tmp_84;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_84 += Conj(Ye(j1,gO1))*Ue(gI1,j1);
      }
      tmp_83 += tmp_84;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_83;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_85;
      std::complex<double> tmp_86;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_86 += Conj(Ve(gI2,j2))*Ye(gO2,j2);
      }
      tmp_85 += tmp_86;
      result += (-0.7071067811865475) * tmp_85;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_87;
      std::complex<double> tmp_88;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_88 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
      }
      tmp_87 += tmp_88;
      result += (-0.7071067811865475) * tmp_87;
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

std::complex<double> CLASSNAME::CpbarUFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -Ye(gO2,gI2);
   }

   return result;
}

double CLASSNAME::CpbarUFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(gI2,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(ZN(gI1,2))*(gYd*KroneckerDelta(0,
      gO2) - g2d*KroneckerDelta(1,gO2)) + Conj(ZN(gI1,3))*(gYu*KroneckerDelta(0,
      gO2) - g2u*KroneckerDelta(1,gO2)) + gYd*Conj(ZN(gI1,0))*KroneckerDelta(2,gO2
      ) - g2d*Conj(ZN(gI1,1))*KroneckerDelta(2,gO2) + gYu*Conj(ZN(gI1,0))*
      KroneckerDelta(3,gO2) - g2u*Conj(ZN(gI1,1))*KroneckerDelta(3,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(2,gO1)*(-(gYd*ZN(gI1,0)
      ) + g2d*ZN(gI1,1)) + KroneckerDelta(3,gO1)*(-(gYu*ZN(gI1,0)) + g2u*ZN(gI1,1)
      ) - gYd*KroneckerDelta(0,gO1)*ZN(gI1,2) + g2d*KroneckerDelta(1,gO1)*ZN(gI1,2
      ) - gYu*KroneckerDelta(0,gO1)*ZN(gI1,3) + g2u*KroneckerDelta(1,gO1)*ZN(gI1,3
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUChiHpChaPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*Conj(UM(gI2,1))*(gYd*KroneckerDelta(0,gO2) +
      g2d*KroneckerDelta(1,gO2)) + g2d*Conj(UM(gI2,0))*KroneckerDelta(2,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpUChiHpChaPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2u*KroneckerDelta(3,gO1)*UP(gI2,0)) - 0.7071067811865475*(gYu*
      KroneckerDelta(0,gO1) + g2u*KroneckerDelta(1,gO1))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiVWpChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(3,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiVWpChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1,gO1) +
      1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI2,2))*(gYd*KroneckerDelta(0,gO2) - g2d*
      KroneckerDelta(1,gO2)) + Conj(ZN(gI2,3))*(-(gYu*KroneckerDelta(0,gO2)) + g2u
      *KroneckerDelta(1,gO2)) + gYd*Conj(ZN(gI2,0))*KroneckerDelta(2,gO2) - g2d*
      Conj(ZN(gI2,1))*KroneckerDelta(2,gO2) - gYu*Conj(ZN(gI2,0))*KroneckerDelta(3
      ,gO2) + g2u*Conj(ZN(gI2,1))*KroneckerDelta(3,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(2,gO1)*(gYd*ZN(gI2,0) - g2d*ZN(gI2,1)) +
      KroneckerDelta(3,gO1)*(-(gYu*ZN(gI2,0)) + g2u*ZN(gI2,1)) + gYd*
      KroneckerDelta(0,gO1)*ZN(gI2,2) - g2d*KroneckerDelta(1,gO1)*ZN(gI2,2) - gYu*
      KroneckerDelta(0,gO1)*ZN(gI2,3) + g2u*KroneckerDelta(1,gO1)*ZN(gI2,3));

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

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2d*Conj(UM(gI1,1))*
      KroneckerDelta(0,gO2) - g2u*Conj(UM(gI1,0))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2d*KroneckerDelta(1,
      gO1)*UP(gI1,0) - g2u*KroneckerDelta(0,gO1)*UP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2d*Conj(UM(gI2,1))*KroneckerDelta(0,gO2) +
      g2u*Conj(UM(gI2,0))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2d*KroneckerDelta(1,gO1)*UP(gI2,0) + g2u*
      KroneckerDelta(0,gO1)*UP(gI2,1));

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

std::complex<double> CLASSNAME::CpbarUChaconjHpChiPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2u*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2)) - 0.7071067811865475*(
      gYu*Conj(ZN(gI2,0)) + g2u*Conj(ZN(gI2,1)))*KroneckerDelta(1,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjHpChiPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*KroneckerDelta(1,gO1)*(gYd*ZN(gI2,0) + g2d*ZN(
      gI2,1)) + g2d*KroneckerDelta(0,gO1)*ZN(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjVWpChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjVWpChiPL(unsigned gO1, unsigned gI2) const
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

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_89;
   std::complex<double> tmp_90;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_91;
      std::complex<double> tmp_92;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_92 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_91 += tmp_92;
      tmp_90 += (Conj(Vd(gI1,j2))) * tmp_91;
   }
   tmp_89 += tmp_90;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_89;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_93;
   std::complex<double> tmp_94;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_95;
      std::complex<double> tmp_96;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_96 += Conj(Yd(j1,j2))*Ud(gI1,j1);
      }
      tmp_95 += tmp_96;
      tmp_94 += (Vd(gO1,j2)) * tmp_95;
   }
   tmp_93 += tmp_94;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_93;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_97;
   std::complex<double> tmp_98;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_99;
      std::complex<double> tmp_100;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_100 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_99 += tmp_100;
      tmp_98 += (Conj(Vd(gI2,j2))) * tmp_99;
   }
   tmp_97 += tmp_98;
   result += (-0.7071067811865475) * tmp_97;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_101;
   std::complex<double> tmp_102;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_103;
      std::complex<double> tmp_104;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_104 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_103 += tmp_104;
      tmp_102 += (Vd(gO1,j2)) * tmp_103;
   }
   tmp_101 += tmp_102;
   result += (-0.7071067811865475) * tmp_101;

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

std::complex<double> CLASSNAME::CpbarFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_105;
   std::complex<double> tmp_106;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_107;
      std::complex<double> tmp_108;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_108 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_107 += tmp_108;
      tmp_106 += (Conj(Vu(gI2,j2))) * tmp_107;
   }
   tmp_105 += tmp_106;
   result += (-1) * tmp_105;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_109;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_111 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_110 += tmp_111;
      tmp_109 += (Vd(gO1,j2)) * tmp_110;
   }
   result += tmp_109;

   return result;
}

double CLASSNAME::CpbarFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_112;
   std::complex<double> tmp_113;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_113 += Conj(Vu(gI2,j1))*Vd(gO1,j1);
   }
   tmp_112 += tmp_113;
   result += (-0.7071067811865475*g2) * tmp_112;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_114;
   std::complex<double> tmp_115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_117 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_116 += tmp_117;
      tmp_115 += (Conj(Ve(gI1,j2))) * tmp_116;
   }
   tmp_114 += tmp_115;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_114;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_118;
   std::complex<double> tmp_119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_121 += Conj(Ye(j1,j2))*Ue(gI1,j1);
      }
      tmp_120 += tmp_121;
      tmp_119 += (Ve(gO1,j2)) * tmp_120;
   }
   tmp_118 += tmp_119;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_118;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_122;
   std::complex<double> tmp_123;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_125;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_125 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_124 += tmp_125;
      tmp_123 += (Conj(Ve(gI2,j2))) * tmp_124;
   }
   tmp_122 += tmp_123;
   result += (-0.7071067811865475) * tmp_122;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_126;
   std::complex<double> tmp_127;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_129;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_129 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_128 += tmp_129;
      tmp_127 += (Ve(gO1,j2)) * tmp_128;
   }
   tmp_126 += tmp_127;
   result += (-0.7071067811865475) * tmp_126;

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

std::complex<double> CLASSNAME::CpbarFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_130;
   std::complex<double> tmp_131;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_131 += Conj(Ue(gO2,j1))*Ye(j1,gI2);
   }
   tmp_130 += tmp_131;
   result += (-1) * tmp_130;

   return result;
}

double CLASSNAME::CpbarFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gO1,gI2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_132;
   std::complex<double> tmp_133;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_135 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_134 += tmp_135;
      tmp_133 += (Conj(Vu(gI1,j2))) * tmp_134;
   }
   tmp_132 += tmp_133;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_132;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_136;
   std::complex<double> tmp_137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_138;
      std::complex<double> tmp_139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_139 += Conj(Yu(j1,j2))*Uu(gI1,j1);
      }
      tmp_138 += tmp_139;
      tmp_137 += (Vu(gO1,j2)) * tmp_138;
   }
   tmp_136 += tmp_137;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_136;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_140;
   std::complex<double> tmp_141;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_142;
      std::complex<double> tmp_143;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_143 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_142 += tmp_143;
      tmp_141 += (Conj(Vu(gI2,j2))) * tmp_142;
   }
   tmp_140 += tmp_141;
   result += (-0.7071067811865475) * tmp_140;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_144;
   std::complex<double> tmp_145;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_146;
      std::complex<double> tmp_147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_147 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_146 += tmp_147;
      tmp_145 += (Vu(gO1,j2)) * tmp_146;
   }
   tmp_144 += tmp_145;
   result += (-0.7071067811865475) * tmp_144;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_148;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_149;
      std::complex<double> tmp_150;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_150 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_149 += tmp_150;
      tmp_148 += (Conj(Vd(gI2,j2))) * tmp_149;
   }
   result += tmp_148;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_151;
   std::complex<double> tmp_152;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_153;
      std::complex<double> tmp_154;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_154 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_153 += tmp_154;
      tmp_152 += (Vu(gO1,j2)) * tmp_153;
   }
   tmp_151 += tmp_152;
   result += (-1) * tmp_151;

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

double CLASSNAME::CpbarFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_155;
   std::complex<double> tmp_156;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_156 += Conj(Vd(gI2,j1))*Vu(gO1,j1);
   }
   tmp_155 += tmp_156;
   result += (-0.7071067811865475*g2) * tmp_155;

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


std::complex<double> CLASSNAME::self_energy_Hp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjHpHphh())*B0(p,MHp,Mhh);
   result += 4*AbsSqr(CpconjHpVWpVP())*(-0.5 + B0(p,0,MVWp));
   result += 4*AbsSqr(CpconjHpVZVWp())*(-0.5 + B0(p,MVWp,MVZ));
   result += -0.5*A0(MAh)*CpHpconjHpAhAh();
   result += -(A0(MHp)*CpHpconjHpconjHpHp());
   result += -0.5*A0(Mhh)*CpHpconjHphhhh();
   result += -(B0(p,MVZ,MVWp)*CpconjHpbargWpCgZ()*CpHpgWpCbargZ());
   result += -(B0(p,MVWp,MVZ)*CpconjHpbargZgWp()*CpHpgZbargWp());
   result += AbsSqr(CpconjHpVWpAh())*F0(p,MAh,MVWp);
   result += AbsSqr(CpconjHpVWphh())*F0(p,Mhh,MVWp);
   result += AbsSqr(CpconjHpVPHp())*F0(p,MHp,0);
   result += AbsSqr(CpconjHpVZHp())*F0(p,MHp,MVZ);
   std::complex<double> tmp_157;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_158;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_158 += (AbsSqr(CpconjHpbarChaChiPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarChaChiPR(gI1,gI2)))*G0(p,MCha(gI1),MChi(gI2));
      }
      tmp_157 += tmp_158;
   }
   result += tmp_157;
   std::complex<double> tmp_159;
   std::complex<double> tmp_160;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_161;
      std::complex<double> tmp_162;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_162 += B0(p,MCha(gI1),MChi(gI2))*(Conj(CpconjHpbarChaChiPR(
            gI1,gI2))*CpconjHpbarChaChiPL(gI1,gI2) + Conj(CpconjHpbarChaChiPL(gI1,
            gI2))*CpconjHpbarChaChiPR(gI1,gI2))*MChi(gI2);
      }
      tmp_161 += tmp_162;
      tmp_160 += (MCha(gI1)) * tmp_161;
   }
   tmp_159 += tmp_160;
   result += (-2) * tmp_159;
   std::complex<double> tmp_163;
   std::complex<double> tmp_164;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_165;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_165 += (AbsSqr(CpconjHpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFdFuPR(gI1,gI2)))*G0(p,MFd(gI1),MFu(gI2));
      }
      tmp_164 += tmp_165;
   }
   tmp_163 += tmp_164;
   result += (3) * tmp_163;
   std::complex<double> tmp_166;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_167;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_167 += (AbsSqr(CpconjHpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFeFvPR(gI1,gI2)))*G0(p,MFe(gI1),MFv(gI2));
      }
      tmp_166 += tmp_167;
   }
   result += tmp_166;
   std::complex<double> tmp_168;
   std::complex<double> tmp_169;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_171;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_171 += B0(p,MFd(gI1),MFu(gI2))*(Conj(CpconjHpbarFdFuPR(gI1,
            gI2))*CpconjHpbarFdFuPL(gI1,gI2) + Conj(CpconjHpbarFdFuPL(gI1,gI2))*
            CpconjHpbarFdFuPR(gI1,gI2))*MFu(gI2);
      }
      tmp_170 += tmp_171;
      tmp_169 += (MFd(gI1)) * tmp_170;
   }
   tmp_168 += tmp_169;
   result += (-6) * tmp_168;
   std::complex<double> tmp_172;
   std::complex<double> tmp_173;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_175;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_175 += B0(p,MFe(gI1),MFv(gI2))*(Conj(CpconjHpbarFeFvPR(gI1,
            gI2))*CpconjHpbarFeFvPL(gI1,gI2) + Conj(CpconjHpbarFeFvPL(gI1,gI2))*
            CpconjHpbarFeFvPR(gI1,gI2))*MFv(gI2);
      }
      tmp_174 += tmp_175;
      tmp_173 += (MFe(gI1)) * tmp_174;
   }
   tmp_172 += tmp_173;
   result += (-2) * tmp_172;
   result += 4*CpHpconjHpconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CpHpconjHpVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CpAhAhAhAh();
   result += -(A0(MHp)*CpAhAhconjHpHp());
   result += -0.5*A0(Mhh)*CpAhAhhhhh();
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpgWp()));
   result += AbsSqr(CpAhhhAh())*B0(p,Mhh,MAh);
   result += AbsSqr(CpAhVZhh())*F0(p,Mhh,MVZ);
   result += 2*AbsSqr(CpAhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_176;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_177;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_177 += (AbsSqr(CpAhbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpAhbarChaChaPR(gI1,gI2)))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_176 += tmp_177;
   }
   result += tmp_176;
   std::complex<double> tmp_178;
   std::complex<double> tmp_179;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_181;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_181 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpAhbarChaChaPR(gI1,
            gI2))*CpAhbarChaChaPL(gI1,gI2) + Conj(CpAhbarChaChaPL(gI1,gI2))*
            CpAhbarChaChaPR(gI1,gI2))*MCha(gI2);
      }
      tmp_180 += tmp_181;
      tmp_179 += (MCha(gI1)) * tmp_180;
   }
   tmp_178 += tmp_179;
   result += (-2) * tmp_178;
   std::complex<double> tmp_182;
   std::complex<double> tmp_183;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_184;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_184 += (AbsSqr(CpAhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpAhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_183 += tmp_184;
   }
   tmp_182 += tmp_183;
   result += (3) * tmp_182;
   std::complex<double> tmp_185;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_186;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_186 += (AbsSqr(CpAhbarFeFePL(gI1,gI2)) + AbsSqr(
            CpAhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_185 += tmp_186;
   }
   result += tmp_185;
   std::complex<double> tmp_187;
   std::complex<double> tmp_188;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_189;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_189 += (AbsSqr(CpAhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpAhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_188 += tmp_189;
   }
   tmp_187 += tmp_188;
   result += (3) * tmp_187;
   std::complex<double> tmp_190;
   std::complex<double> tmp_191;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_192;
      std::complex<double> tmp_193;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_193 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpAhbarFdFdPR(gI1,gI2))
            *CpAhbarFdFdPL(gI1,gI2) + Conj(CpAhbarFdFdPL(gI1,gI2))*CpAhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_192 += tmp_193;
      tmp_191 += (MFd(gI1)) * tmp_192;
   }
   tmp_190 += tmp_191;
   result += (-6) * tmp_190;
   std::complex<double> tmp_194;
   std::complex<double> tmp_195;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_196;
      std::complex<double> tmp_197;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_197 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpAhbarFeFePR(gI1,gI2))
            *CpAhbarFeFePL(gI1,gI2) + Conj(CpAhbarFeFePL(gI1,gI2))*CpAhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_196 += tmp_197;
      tmp_195 += (MFe(gI1)) * tmp_196;
   }
   tmp_194 += tmp_195;
   result += (-2) * tmp_194;
   std::complex<double> tmp_198;
   std::complex<double> tmp_199;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_200;
      std::complex<double> tmp_201;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_201 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpAhbarFuFuPR(gI1,gI2))
            *CpAhbarFuFuPL(gI1,gI2) + Conj(CpAhbarFuFuPL(gI1,gI2))*CpAhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_200 += tmp_201;
      tmp_199 += (MFu(gI1)) * tmp_200;
   }
   tmp_198 += tmp_199;
   result += (-6) * tmp_198;
   std::complex<double> tmp_202;
   std::complex<double> tmp_203;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_204;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_204 += (AbsSqr(CpAhChiChiPL(gI1,gI2)) + AbsSqr(CpAhChiChiPR(
            gI1,gI2)))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_203 += tmp_204;
   }
   tmp_202 += tmp_203;
   result += (0.5) * tmp_202;
   std::complex<double> tmp_205;
   std::complex<double> tmp_206;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_207;
      std::complex<double> tmp_208;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_208 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpAhChiChiPR(gI1,gI2)
            )*CpAhChiChiPL(gI1,gI2) + Conj(CpAhChiChiPL(gI1,gI2))*CpAhChiChiPR(gI1
            ,gI2))*MChi(gI2);
      }
      tmp_207 += tmp_208;
      tmp_206 += (MChi(gI1)) * tmp_207;
   }
   tmp_205 += tmp_206;
   result += (-1) * tmp_205;
   result += 4*CpAhAhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CpAhAhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CphhAhAh())*B0(p,MAh,MAh);
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpgWp()));
   result += -(B0(p,MVZ,MVZ)*Sqr(CphhbargZgZ()));
   result += AbsSqr(CphhconjHpHp())*B0(p,MHp,MHp);
   result += 4*AbsSqr(CphhconjVWpVWp())*(-0.5 + B0(p,MVWp,MVWp));
   result += -0.5*A0(MAh)*CphhhhAhAh();
   result += -(A0(MHp)*CphhhhconjHpHp());
   result += 0.5*AbsSqr(Cphhhhhh())*B0(p,Mhh,Mhh);
   result += -0.5*A0(Mhh)*Cphhhhhhhh();
   result += 2*AbsSqr(CphhVZVZ())*(-0.5 + B0(p,MVZ,MVZ));
   result += AbsSqr(CphhVZAh())*F0(p,MAh,MVZ);
   result += 2*AbsSqr(CphhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_209;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_210;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_210 += (AbsSqr(CphhbarChaChaPL(gI1,gI2)) + AbsSqr(
            CphhbarChaChaPR(gI1,gI2)))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_209 += tmp_210;
   }
   result += tmp_209;
   std::complex<double> tmp_211;
   std::complex<double> tmp_212;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_213;
      std::complex<double> tmp_214;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_214 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CphhbarChaChaPR(gI1,
            gI2))*CphhbarChaChaPL(gI1,gI2) + Conj(CphhbarChaChaPL(gI1,gI2))*
            CphhbarChaChaPR(gI1,gI2))*MCha(gI2);
      }
      tmp_213 += tmp_214;
      tmp_212 += (MCha(gI1)) * tmp_213;
   }
   tmp_211 += tmp_212;
   result += (-2) * tmp_211;
   std::complex<double> tmp_215;
   std::complex<double> tmp_216;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_217;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_217 += (AbsSqr(CphhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CphhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_216 += tmp_217;
   }
   tmp_215 += tmp_216;
   result += (3) * tmp_215;
   std::complex<double> tmp_218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_219;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_219 += (AbsSqr(CphhbarFeFePL(gI1,gI2)) + AbsSqr(
            CphhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_218 += tmp_219;
   }
   result += tmp_218;
   std::complex<double> tmp_220;
   std::complex<double> tmp_221;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_222;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_222 += (AbsSqr(CphhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CphhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_221 += tmp_222;
   }
   tmp_220 += tmp_221;
   result += (3) * tmp_220;
   std::complex<double> tmp_223;
   std::complex<double> tmp_224;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_225;
      std::complex<double> tmp_226;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_226 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CphhbarFdFdPR(gI1,gI2))
            *CphhbarFdFdPL(gI1,gI2) + Conj(CphhbarFdFdPL(gI1,gI2))*CphhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_225 += tmp_226;
      tmp_224 += (MFd(gI1)) * tmp_225;
   }
   tmp_223 += tmp_224;
   result += (-6) * tmp_223;
   std::complex<double> tmp_227;
   std::complex<double> tmp_228;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_229;
      std::complex<double> tmp_230;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_230 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CphhbarFeFePR(gI1,gI2))
            *CphhbarFeFePL(gI1,gI2) + Conj(CphhbarFeFePL(gI1,gI2))*CphhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_229 += tmp_230;
      tmp_228 += (MFe(gI1)) * tmp_229;
   }
   tmp_227 += tmp_228;
   result += (-2) * tmp_227;
   std::complex<double> tmp_231;
   std::complex<double> tmp_232;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_233;
      std::complex<double> tmp_234;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_234 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CphhbarFuFuPR(gI1,gI2))
            *CphhbarFuFuPL(gI1,gI2) + Conj(CphhbarFuFuPL(gI1,gI2))*CphhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_233 += tmp_234;
      tmp_232 += (MFu(gI1)) * tmp_233;
   }
   tmp_231 += tmp_232;
   result += (-6) * tmp_231;
   std::complex<double> tmp_235;
   std::complex<double> tmp_236;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_237;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_237 += (AbsSqr(CphhChiChiPL(gI1,gI2)) + AbsSqr(CphhChiChiPR(
            gI1,gI2)))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_236 += tmp_237;
   }
   tmp_235 += tmp_236;
   result += (0.5) * tmp_235;
   std::complex<double> tmp_238;
   std::complex<double> tmp_239;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_240;
      std::complex<double> tmp_241;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_241 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CphhChiChiPR(gI1,gI2)
            )*CphhChiChiPL(gI1,gI2) + Conj(CphhChiChiPL(gI1,gI2))*CphhChiChiPR(gI1
            ,gI2))*MChi(gI2);
      }
      tmp_240 += tmp_241;
      tmp_239 += (MChi(gI1)) * tmp_240;
   }
   tmp_238 += tmp_239;
   result += (-1) * tmp_238;
   result += 4*CphhhhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CphhhhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWpCgWpC())*B00(p,MVWp,MVWp);
   result += AbsSqr(CpVZbargWpgWp())*B00(p,MVWp,MVWp);
   result += -4*AbsSqr(CpVZconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVZconjVWpHp())*B0(p,MVWp,MHp);
   result += -4*AbsSqr(CpVZhhAh())*B00(p,MAh,Mhh);
   result += 0.5*A0(MAh)*CpVZVZAhAh();
   result += A0(MHp)*CpVZVZconjHpHp();
   result += -(A0(MVWp)*(4*CpVZVZconjVWpVWp1() + CpVZVZconjVWpVWp2() +
      CpVZVZconjVWpVWp3()));
   result += AbsSqr(CpVZVZhh())*B0(p,MVZ,Mhh);
   result += 0.5*A0(Mhh)*CpVZVZhhhh();
   std::complex<double> tmp_242;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_243;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_243 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_243 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_242 += tmp_243;
   }
   result += tmp_242;
   std::complex<double> tmp_244;
   std::complex<double> tmp_245;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_246;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_246 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_246 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_245 += tmp_246;
   }
   tmp_244 += tmp_245;
   result += (3) * tmp_244;
   std::complex<double> tmp_247;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_248;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_248 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_248 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_247 += tmp_248;
   }
   result += tmp_247;
   std::complex<double> tmp_249;
   std::complex<double> tmp_250;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_251;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_251 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_251 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_250 += tmp_251;
   }
   tmp_249 += tmp_250;
   result += (3) * tmp_249;
   std::complex<double> tmp_252;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_253;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_253 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_253 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_252 += tmp_253;
   }
   result += tmp_252;
   std::complex<double> tmp_254;
   std::complex<double> tmp_255;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_256;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_256 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_256 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_255 += tmp_256;
   }
   tmp_254 += tmp_255;
   result += (0.5) * tmp_254;
   result += 2*CpVZVZconjVWpVWp1()*Sqr(MVWp);
   result += -(AbsSqr(CpVZconjVWpVWp())*(2*A0(MVWp) + 10*B00(p,MVWp,MVWp) - 2*(
      2*Sqr(MVWp) - 0.3333333333333333*Sqr(p)) + B0(p,MVWp,MVWp)*(2*Sqr(MVWp) + 4*
      Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWpbargPgWp())*B00(p,MVWp,MVP);
   result += AbsSqr(CpconjVWpbargWpCgP())*B00(p,MVP,MVWp);
   result += AbsSqr(CpconjVWpbargWpCgZ())*B00(p,MVZ,MVWp);
   result += AbsSqr(CpconjVWpbargZgWp())*B00(p,MVWp,MVZ);
   result += -4*AbsSqr(CpconjVWpHpAh())*B00(p,MAh,MHp);
   result += -4*AbsSqr(CpconjVWpHphh())*B00(p,Mhh,MHp);
   result += AbsSqr(CpconjVWpVPHp())*B0(p,0,MHp);
   result += AbsSqr(CpconjVWpVWphh())*B0(p,MVWp,Mhh);
   result += AbsSqr(CpconjVWpVZHp())*B0(p,MVZ,MHp);
   result += 0.5*A0(MAh)*CpVWpconjVWpAhAh();
   result += A0(MHp)*CpVWpconjVWpconjHpHp();
   result += -(A0(MVWp)*(4*CpVWpconjVWpconjVWpVWp1() + CpVWpconjVWpconjVWpVWp2(
      ) + CpVWpconjVWpconjVWpVWp3()));
   result += 0.5*A0(Mhh)*CpVWpconjVWphhhh();
   result += 0;
   std::complex<double> tmp_257;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_258;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_258 += (AbsSqr(CpconjVWpbarChaChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarChaChiPR(gI1,gI2)))*H0(p,MCha(gI1),MChi(gI2));
         tmp_258 += 4*B0(p,MCha(gI1),MChi(gI2))*MCha(gI1)*MChi(gI2)*Re(
            Conj(CpconjVWpbarChaChiPL(gI1,gI2))*CpconjVWpbarChaChiPR(gI1,gI2));
      }
      tmp_257 += tmp_258;
   }
   result += tmp_257;
   std::complex<double> tmp_259;
   std::complex<double> tmp_260;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_261;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_261 += (AbsSqr(CpconjVWpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFdFuPR(gI1,gI2)))*H0(p,MFd(gI1),MFu(gI2));
         tmp_261 += 4*B0(p,MFd(gI1),MFu(gI2))*MFd(gI1)*MFu(gI2)*Re(Conj(
            CpconjVWpbarFdFuPL(gI1,gI2))*CpconjVWpbarFdFuPR(gI1,gI2));
      }
      tmp_260 += tmp_261;
   }
   tmp_259 += tmp_260;
   result += (3) * tmp_259;
   std::complex<double> tmp_262;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_263;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_263 += (AbsSqr(CpconjVWpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFeFvPR(gI1,gI2)))*H0(p,MFe(gI1),MFv(gI2));
         tmp_263 += 4*B0(p,MFe(gI1),MFv(gI2))*MFe(gI1)*MFv(gI2)*Re(Conj(
            CpconjVWpbarFeFvPL(gI1,gI2))*CpconjVWpbarFeFvPR(gI1,gI2));
      }
      tmp_262 += tmp_263;
   }
   result += tmp_262;
   result += 2*CpVWpconjVWpconjVWpVWp1()*Sqr(MVWp);
   result += -(AbsSqr(CpconjVWpVWpVP())*(A0(MVWp) + 10*B00(p,MVWp,0) - 2*(Sqr(
      MVWp) - 0.3333333333333333*Sqr(p)) + B0(p,MVWp,0)*(Sqr(MVWp) + 4*Sqr(p))));
   result += 0.5*(-(A0(MVZ)*(4*CpVWpconjVWpVZVZ1() + CpVWpconjVWpVZVZ2() +
      CpVWpconjVWpVZVZ3())) + 2*CpVWpconjVWpVZVZ1()*Sqr(MVZ));
   result += -(AbsSqr(CpconjVWpVZVWp())*(A0(MVWp) + A0(MVZ) + 10*B00(p,MVZ,MVWp
      ) - 2*(Sqr(MVWp) + Sqr(MVZ) - 0.3333333333333333*Sqr(p)) + B0(p,MVZ,MVWp)*(
      Sqr(MVWp) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_264;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_264 += B0(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_264;
   std::complex<double> tmp_265;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_265 += B0(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_265;
   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_267 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_266 += tmp_267;
   result += (-5.333333333333333) * tmp_266;
   std::complex<double> tmp_268;
   std::complex<double> tmp_269;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_269 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_268 += tmp_269;
   result += (-4) * tmp_268;
   std::complex<double> tmp_270;
   std::complex<double> tmp_271;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_271 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_270 += tmp_271;
   result += (-4) * tmp_270;
   std::complex<double> tmp_272;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_272 += B0(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_272;
   std::complex<double> tmp_273;
   std::complex<double> tmp_274;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_274 += (-0.5 + B0(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPR(gO2,
         gI2))*CpbarUFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_273 += tmp_274;
   result += (-4) * tmp_273;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_275;
   std::complex<double> tmp_276;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_276 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPR(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1);
   }
   tmp_275 += tmp_276;
   result += (-0.5) * tmp_275;
   std::complex<double> tmp_277;
   std::complex<double> tmp_278;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_278 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPR(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2);
   }
   tmp_277 += tmp_278;
   result += (-0.5) * tmp_277;
   std::complex<double> tmp_279;
   std::complex<double> tmp_280;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_280 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPL(gO2,
         gI2))*CpbarUFdconjVWpFuPL(gO1,gI2);
   }
   tmp_279 += tmp_280;
   result += (-1) * tmp_279;
   std::complex<double> tmp_281;
   std::complex<double> tmp_282;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_282 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPR(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2);
   }
   tmp_281 += tmp_282;
   result += (-0.5) * tmp_281;
   std::complex<double> tmp_283;
   std::complex<double> tmp_284;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_284 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_283 += tmp_284;
   result += (-1.3333333333333333) * tmp_283;
   std::complex<double> tmp_285;
   std::complex<double> tmp_286;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_286 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_285 += tmp_286;
   result += (-1) * tmp_285;
   std::complex<double> tmp_287;
   std::complex<double> tmp_288;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_288 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_287 += tmp_288;
   result += (-1) * tmp_287;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_289;
   std::complex<double> tmp_290;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_290 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPL(gO1,gI1);
   }
   tmp_289 += tmp_290;
   result += (-0.5) * tmp_289;
   std::complex<double> tmp_291;
   std::complex<double> tmp_292;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_292 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPL(gO1,gI2);
   }
   tmp_291 += tmp_292;
   result += (-0.5) * tmp_291;
   std::complex<double> tmp_293;
   std::complex<double> tmp_294;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_294 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPR(gO2,
         gI2))*CpbarUFdconjVWpFuPR(gO1,gI2);
   }
   tmp_293 += tmp_294;
   result += (-1) * tmp_293;
   std::complex<double> tmp_295;
   std::complex<double> tmp_296;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_296 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPL(gO1,gI2);
   }
   tmp_295 += tmp_296;
   result += (-0.5) * tmp_295;
   std::complex<double> tmp_297;
   std::complex<double> tmp_298;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_298 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_297 += tmp_298;
   result += (-1.3333333333333333) * tmp_297;
   std::complex<double> tmp_299;
   std::complex<double> tmp_300;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_300 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_299 += tmp_300;
   result += (-1) * tmp_299;
   std::complex<double> tmp_301;
   std::complex<double> tmp_302;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_302 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_301 += tmp_302;
   result += (-1) * tmp_301;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_303;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_303 += B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_303;
   std::complex<double> tmp_304;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_304 += B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_304;
   std::complex<double> tmp_305;
   std::complex<double> tmp_306;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_306 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))
         *CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_305 += tmp_306;
   result += (-4) * tmp_305;
   std::complex<double> tmp_307;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_307 += B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_307;
   std::complex<double> tmp_308;
   std::complex<double> tmp_309;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_309 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_308 += tmp_309;
   result += (-5.333333333333333) * tmp_308;
   std::complex<double> tmp_310;
   std::complex<double> tmp_311;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_311 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_310 += tmp_311;
   result += (-4) * tmp_310;
   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_313 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_312 += tmp_313;
   result += (-4) * tmp_312;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_314;
   std::complex<double> tmp_315;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_315 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1);
   }
   tmp_314 += tmp_315;
   result += (-0.5) * tmp_314;
   std::complex<double> tmp_316;
   std::complex<double> tmp_317;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_317 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2);
   }
   tmp_316 += tmp_317;
   result += (-0.5) * tmp_316;
   std::complex<double> tmp_318;
   std::complex<double> tmp_319;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_319 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2);
   }
   tmp_318 += tmp_319;
   result += (-0.5) * tmp_318;
   std::complex<double> tmp_320;
   std::complex<double> tmp_321;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_321 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_320 += tmp_321;
   result += (-1.3333333333333333) * tmp_320;
   std::complex<double> tmp_322;
   std::complex<double> tmp_323;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_323 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_322 += tmp_323;
   result += (-1) * tmp_322;
   std::complex<double> tmp_324;
   std::complex<double> tmp_325;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_325 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2);
   }
   tmp_324 += tmp_325;
   result += (-1) * tmp_324;
   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_327 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_326 += tmp_327;
   result += (-1) * tmp_326;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_328;
   std::complex<double> tmp_329;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_329 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPL(gO1,gI1);
   }
   tmp_328 += tmp_329;
   result += (-0.5) * tmp_328;
   std::complex<double> tmp_330;
   std::complex<double> tmp_331;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_331 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPL(gO1,gI2);
   }
   tmp_330 += tmp_331;
   result += (-0.5) * tmp_330;
   std::complex<double> tmp_332;
   std::complex<double> tmp_333;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_333 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPL(gO1,gI2);
   }
   tmp_332 += tmp_333;
   result += (-0.5) * tmp_332;
   std::complex<double> tmp_334;
   std::complex<double> tmp_335;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_335 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_334 += tmp_335;
   result += (-1.3333333333333333) * tmp_334;
   std::complex<double> tmp_336;
   std::complex<double> tmp_337;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_337 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_336 += tmp_337;
   result += (-1) * tmp_336;
   std::complex<double> tmp_338;
   std::complex<double> tmp_339;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_339 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPR(gO1,gI2);
   }
   tmp_338 += tmp_339;
   result += (-1) * tmp_338;
   std::complex<double> tmp_340;
   std::complex<double> tmp_341;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_341 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_340 += tmp_341;
   result += (-1) * tmp_340;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_342;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_342 += B0(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_342;
   std::complex<double> tmp_343;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_343 += B0(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_343;
   std::complex<double> tmp_344;
   std::complex<double> tmp_345;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_345 += (-0.5 + B0(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_344 += tmp_345;
   result += (-4) * tmp_344;
   std::complex<double> tmp_346;
   std::complex<double> tmp_347;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_347 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_346 += tmp_347;
   result += (-4) * tmp_346;
   std::complex<double> tmp_348;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_348 += B0(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_348;
   std::complex<double> tmp_349;
   std::complex<double> tmp_350;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_350 += (-0.5 + B0(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPR(gO2,
         gI2))*CpbarUFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_349 += tmp_350;
   result += (-4) * tmp_349;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_351;
   std::complex<double> tmp_352;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_352 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPR(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1);
   }
   tmp_351 += tmp_352;
   result += (-0.5) * tmp_351;
   std::complex<double> tmp_353;
   std::complex<double> tmp_354;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_354 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPR(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2);
   }
   tmp_353 += tmp_354;
   result += (-0.5) * tmp_353;
   std::complex<double> tmp_355;
   std::complex<double> tmp_356;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_356 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPL(gO2,
         gI2))*CpbarUFeconjVWpFvPL(gO1,gI2);
   }
   tmp_355 += tmp_356;
   result += (-1) * tmp_355;
   std::complex<double> tmp_357;
   std::complex<double> tmp_358;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_358 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePR(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2);
   }
   tmp_357 += tmp_358;
   result += (-0.5) * tmp_357;
   std::complex<double> tmp_359;
   std::complex<double> tmp_360;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_360 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_359 += tmp_360;
   result += (-1) * tmp_359;
   std::complex<double> tmp_361;
   std::complex<double> tmp_362;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_362 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_361 += tmp_362;
   result += (-1) * tmp_361;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_363;
   std::complex<double> tmp_364;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_364 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPL(gO1,gI1);
   }
   tmp_363 += tmp_364;
   result += (-0.5) * tmp_363;
   std::complex<double> tmp_365;
   std::complex<double> tmp_366;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_366 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPL(gO1,gI2);
   }
   tmp_365 += tmp_366;
   result += (-0.5) * tmp_365;
   std::complex<double> tmp_367;
   std::complex<double> tmp_368;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_368 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPR(gO2,
         gI2))*CpbarUFeconjVWpFvPR(gO1,gI2);
   }
   tmp_367 += tmp_368;
   result += (-1) * tmp_367;
   std::complex<double> tmp_369;
   std::complex<double> tmp_370;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_370 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePL(gO1,gI2);
   }
   tmp_369 += tmp_370;
   result += (-0.5) * tmp_369;
   std::complex<double> tmp_371;
   std::complex<double> tmp_372;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_372 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_371 += tmp_372;
   result += (-1) * tmp_371;
   std::complex<double> tmp_373;
   std::complex<double> tmp_374;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_374 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_373 += tmp_374;
   result += (-1) * tmp_373;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_375;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_375 += B0(p,MChi(gI1),MAh)*Conj(CpUChiChiAhPL(gO2,gI1))*
         CpUChiChiAhPR(gO1,gI1)*MChi(gI1);
   }
   result += tmp_375;
   std::complex<double> tmp_376;
   std::complex<double> tmp_377;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_377 += B0(p,MCha(gI2),MHp)*Conj(CpUChiHpChaPL(gO2,gI2))*
         CpUChiHpChaPR(gO1,gI2)*MCha(gI2);
   }
   tmp_376 += tmp_377;
   result += (2) * tmp_376;
   std::complex<double> tmp_378;
   std::complex<double> tmp_379;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_379 += (-0.5 + B0(p,MCha(gI2),MVWp))*Conj(CpUChiVWpChaPR(gO2,gI2))
         *CpUChiVWpChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_378 += tmp_379;
   result += (-8) * tmp_378;
   std::complex<double> tmp_380;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_380 += B0(p,MChi(gI2),Mhh)*Conj(CpUChihhChiPL(gO2,gI2))*
         CpUChihhChiPR(gO1,gI2)*MChi(gI2);
   }
   result += tmp_380;
   std::complex<double> tmp_381;
   std::complex<double> tmp_382;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_382 += (-0.5 + B0(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_381 += tmp_382;
   result += (-4) * tmp_381;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_383;
   std::complex<double> tmp_384;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_384 += B1(p,MChi(gI1),MAh)*Conj(CpUChiChiAhPR(gO2,gI1))*
         CpUChiChiAhPR(gO1,gI1);
   }
   tmp_383 += tmp_384;
   result += (-0.5) * tmp_383;
   std::complex<double> tmp_385;
   std::complex<double> tmp_386;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_386 += B1(p,MCha(gI2),MHp)*Conj(CpUChiHpChaPR(gO2,gI2))*
         CpUChiHpChaPR(gO1,gI2);
   }
   tmp_385 += tmp_386;
   result += (-1) * tmp_385;
   std::complex<double> tmp_387;
   std::complex<double> tmp_388;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_388 += (0.5 + B1(p,MCha(gI2),MVWp))*Conj(CpUChiVWpChaPL(gO2,gI2))*
         CpUChiVWpChaPL(gO1,gI2);
   }
   tmp_387 += tmp_388;
   result += (-2) * tmp_387;
   std::complex<double> tmp_389;
   std::complex<double> tmp_390;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_390 += B1(p,MChi(gI2),Mhh)*Conj(CpUChihhChiPR(gO2,gI2))*
         CpUChihhChiPR(gO1,gI2);
   }
   tmp_389 += tmp_390;
   result += (-0.5) * tmp_389;
   std::complex<double> tmp_391;
   std::complex<double> tmp_392;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_392 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_391 += tmp_392;
   result += (-1) * tmp_391;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_393;
   std::complex<double> tmp_394;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_394 += B1(p,MChi(gI1),MAh)*Conj(CpUChiChiAhPL(gO2,gI1))*
         CpUChiChiAhPL(gO1,gI1);
   }
   tmp_393 += tmp_394;
   result += (-0.5) * tmp_393;
   std::complex<double> tmp_395;
   std::complex<double> tmp_396;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_396 += B1(p,MCha(gI2),MHp)*Conj(CpUChiHpChaPL(gO2,gI2))*
         CpUChiHpChaPL(gO1,gI2);
   }
   tmp_395 += tmp_396;
   result += (-1) * tmp_395;
   std::complex<double> tmp_397;
   std::complex<double> tmp_398;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_398 += (0.5 + B1(p,MCha(gI2),MVWp))*Conj(CpUChiVWpChaPR(gO2,gI2))*
         CpUChiVWpChaPR(gO1,gI2);
   }
   tmp_397 += tmp_398;
   result += (-2) * tmp_397;
   std::complex<double> tmp_399;
   std::complex<double> tmp_400;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_400 += B1(p,MChi(gI2),Mhh)*Conj(CpUChihhChiPL(gO2,gI2))*
         CpUChihhChiPL(gO1,gI2);
   }
   tmp_399 += tmp_400;
   result += (-0.5) * tmp_399;
   std::complex<double> tmp_401;
   std::complex<double> tmp_402;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_402 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_401 += tmp_402;
   result += (-1) * tmp_401;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_403;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_403 += B0(p,MCha(gI1),MAh)*Conj(CpbarUChaChaAhPL(gO2,gI1))*
         CpbarUChaChaAhPR(gO1,gI1)*MCha(gI1);
   }
   result += tmp_403;
   std::complex<double> tmp_404;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_404 += B0(p,MCha(gI2),Mhh)*Conj(CpbarUChahhChaPL(gO2,gI2))*
         CpbarUChahhChaPR(gO1,gI2)*MCha(gI2);
   }
   result += tmp_404;
   std::complex<double> tmp_405;
   std::complex<double> tmp_406;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_406 += (-0.5 + B0(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_405 += tmp_406;
   result += (-4) * tmp_405;
   std::complex<double> tmp_407;
   std::complex<double> tmp_408;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_408 += (-0.5 + B0(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2)
         )*CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_407 += tmp_408;
   result += (-4) * tmp_407;
   std::complex<double> tmp_409;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_409 += B0(p,MChi(gI2),MHp)*Conj(CpbarUChaconjHpChiPL(gO2,gI2))*
         CpbarUChaconjHpChiPR(gO1,gI2)*MChi(gI2);
   }
   result += tmp_409;
   std::complex<double> tmp_410;
   std::complex<double> tmp_411;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_411 += (-0.5 + B0(p,MChi(gI2),MVWp))*Conj(CpbarUChaconjVWpChiPR(
         gO2,gI2))*CpbarUChaconjVWpChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_410 += tmp_411;
   result += (-4) * tmp_410;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_413 += B1(p,MCha(gI1),MAh)*Conj(CpbarUChaChaAhPR(gO2,gI1))*
         CpbarUChaChaAhPR(gO1,gI1);
   }
   tmp_412 += tmp_413;
   result += (-0.5) * tmp_412;
   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_415 += B1(p,MCha(gI2),Mhh)*Conj(CpbarUChahhChaPR(gO2,gI2))*
         CpbarUChahhChaPR(gO1,gI2);
   }
   tmp_414 += tmp_415;
   result += (-0.5) * tmp_414;
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_417 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_416 += tmp_417;
   result += (-1) * tmp_416;
   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_419 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPL(gO2,gI2))
         *CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_418 += tmp_419;
   result += (-1) * tmp_418;
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_421 += B1(p,MChi(gI2),MHp)*Conj(CpbarUChaconjHpChiPR(gO2,gI2))*
         CpbarUChaconjHpChiPR(gO1,gI2);
   }
   tmp_420 += tmp_421;
   result += (-0.5) * tmp_420;
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_423 += (0.5 + B1(p,MChi(gI2),MVWp))*Conj(CpbarUChaconjVWpChiPL(gO2
         ,gI2))*CpbarUChaconjVWpChiPL(gO1,gI2);
   }
   tmp_422 += tmp_423;
   result += (-1) * tmp_422;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_424;
   std::complex<double> tmp_425;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_425 += B1(p,MCha(gI1),MAh)*Conj(CpbarUChaChaAhPL(gO2,gI1))*
         CpbarUChaChaAhPL(gO1,gI1);
   }
   tmp_424 += tmp_425;
   result += (-0.5) * tmp_424;
   std::complex<double> tmp_426;
   std::complex<double> tmp_427;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_427 += B1(p,MCha(gI2),Mhh)*Conj(CpbarUChahhChaPL(gO2,gI2))*
         CpbarUChahhChaPL(gO1,gI2);
   }
   tmp_426 += tmp_427;
   result += (-0.5) * tmp_426;
   std::complex<double> tmp_428;
   std::complex<double> tmp_429;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_429 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_428 += tmp_429;
   result += (-1) * tmp_428;
   std::complex<double> tmp_430;
   std::complex<double> tmp_431;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_431 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2))
         *CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_430 += tmp_431;
   result += (-1) * tmp_430;
   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_433 += B1(p,MChi(gI2),MHp)*Conj(CpbarUChaconjHpChiPL(gO2,gI2))*
         CpbarUChaconjHpChiPL(gO1,gI2);
   }
   tmp_432 += tmp_433;
   result += (-0.5) * tmp_432;
   std::complex<double> tmp_434;
   std::complex<double> tmp_435;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_435 += (0.5 + B1(p,MChi(gI2),MVWp))*Conj(CpbarUChaconjVWpChiPR(gO2
         ,gI2))*CpbarUChaconjVWpChiPR(gO1,gI2);
   }
   tmp_434 += tmp_435;
   result += (-1) * tmp_434;

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

   std::complex<double> tmp_436;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_437;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_437 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_437 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_436 += tmp_437;
   }
   result += tmp_436;
   std::complex<double> tmp_438;
   std::complex<double> tmp_439;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_440;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_440 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_440 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_439 += tmp_440;
   }
   tmp_438 += tmp_439;
   result += (0.5) * tmp_438;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_441;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_442;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_442 += (AbsSqr(CpconjVWpbarChaChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarChaChiPR(gI1,gI2)))*H0(p,MCha(gI1),MChi(gI2));
         tmp_442 += 4*B0(p,MCha(gI1),MChi(gI2))*MCha(gI1)*MChi(gI2)*Re(
            Conj(CpconjVWpbarChaChiPL(gI1,gI2))*CpconjVWpbarChaChiPR(gI1,gI2));
      }
      tmp_441 += tmp_442;
   }
   result += tmp_441;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_443;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_443 += B0(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_443;
   std::complex<double> tmp_444;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_444 += B0(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_444;
   std::complex<double> tmp_445;
   std::complex<double> tmp_446;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_446 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_445 += tmp_446;
   result += (-4) * tmp_445;
   std::complex<double> tmp_447;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_447 += B0(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_447;
   std::complex<double> tmp_448;
   std::complex<double> tmp_449;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_449 += (-0.5 + B0(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPR(gO2,
         gI2))*CpbarFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_448 += tmp_449;
   result += (-4) * tmp_448;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_450;
   std::complex<double> tmp_451;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_451 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPR(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1);
   }
   tmp_450 += tmp_451;
   result += (-0.5) * tmp_450;
   std::complex<double> tmp_452;
   std::complex<double> tmp_453;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_453 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPR(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2);
   }
   tmp_452 += tmp_453;
   result += (-0.5) * tmp_452;
   std::complex<double> tmp_454;
   std::complex<double> tmp_455;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_455 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPL(gO2,gI2
         ))*CpbarFdconjVWpFuPL(gO1,gI2);
   }
   tmp_454 += tmp_455;
   result += (-1) * tmp_454;
   std::complex<double> tmp_456;
   std::complex<double> tmp_457;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_457 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPR(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2);
   }
   tmp_456 += tmp_457;
   result += (-0.5) * tmp_456;
   std::complex<double> tmp_458;
   std::complex<double> tmp_459;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_459 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_458 += tmp_459;
   result += (-1) * tmp_458;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_460;
   std::complex<double> tmp_461;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_461 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPL(gO1,gI1);
   }
   tmp_460 += tmp_461;
   result += (-0.5) * tmp_460;
   std::complex<double> tmp_462;
   std::complex<double> tmp_463;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_463 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPL(gO1,gI2);
   }
   tmp_462 += tmp_463;
   result += (-0.5) * tmp_462;
   std::complex<double> tmp_464;
   std::complex<double> tmp_465;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_465 += (0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPR(gO2,gI2
         ))*CpbarFdconjVWpFuPR(gO1,gI2);
   }
   tmp_464 += tmp_465;
   result += (-1) * tmp_464;
   std::complex<double> tmp_466;
   std::complex<double> tmp_467;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_467 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPL(gO1,gI2);
   }
   tmp_466 += tmp_467;
   result += (-0.5) * tmp_466;
   std::complex<double> tmp_468;
   std::complex<double> tmp_469;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_469 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_468 += tmp_469;
   result += (-1) * tmp_468;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_470;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_470 += B0(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_470;
   std::complex<double> tmp_471;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_471 += B0(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_471;
   std::complex<double> tmp_472;
   std::complex<double> tmp_473;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_473 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_472 += tmp_473;
   result += (-4) * tmp_472;
   std::complex<double> tmp_474;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_474 += B0(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_474;
   std::complex<double> tmp_475;
   std::complex<double> tmp_476;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_476 += (-0.5 + B0(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPR(gO2,
         gI2))*CpbarFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_475 += tmp_476;
   result += (-4) * tmp_475;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_477;
   std::complex<double> tmp_478;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_478 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPR(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1);
   }
   tmp_477 += tmp_478;
   result += (-0.5) * tmp_477;
   std::complex<double> tmp_479;
   std::complex<double> tmp_480;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_480 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPR(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2);
   }
   tmp_479 += tmp_480;
   result += (-0.5) * tmp_479;
   std::complex<double> tmp_481;
   std::complex<double> tmp_482;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_482 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPL(gO2,gI2
         ))*CpbarFeconjVWpFvPL(gO1,gI2);
   }
   tmp_481 += tmp_482;
   result += (-1) * tmp_481;
   std::complex<double> tmp_483;
   std::complex<double> tmp_484;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_484 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePR(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2);
   }
   tmp_483 += tmp_484;
   result += (-0.5) * tmp_483;
   std::complex<double> tmp_485;
   std::complex<double> tmp_486;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_486 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_485 += tmp_486;
   result += (-1) * tmp_485;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_487;
   std::complex<double> tmp_488;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_488 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPL(gO1,gI1);
   }
   tmp_487 += tmp_488;
   result += (-0.5) * tmp_487;
   std::complex<double> tmp_489;
   std::complex<double> tmp_490;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_490 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPL(gO1,gI2);
   }
   tmp_489 += tmp_490;
   result += (-0.5) * tmp_489;
   std::complex<double> tmp_491;
   std::complex<double> tmp_492;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_492 += (0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPR(gO2,gI2
         ))*CpbarFeconjVWpFvPR(gO1,gI2);
   }
   tmp_491 += tmp_492;
   result += (-1) * tmp_491;
   std::complex<double> tmp_493;
   std::complex<double> tmp_494;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_494 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePL(gO1,gI2);
   }
   tmp_493 += tmp_494;
   result += (-0.5) * tmp_493;
   std::complex<double> tmp_495;
   std::complex<double> tmp_496;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_496 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_495 += tmp_496;
   result += (-1) * tmp_495;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_497;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_497 += B0(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_497;
   std::complex<double> tmp_498;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_498 += B0(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_498;
   std::complex<double> tmp_499;
   std::complex<double> tmp_500;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_500 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_499 += tmp_500;
   result += (-4) * tmp_499;
   std::complex<double> tmp_501;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_501 += B0(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_501;
   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_503 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_502 += tmp_503;
   result += (-4) * tmp_502;
   std::complex<double> tmp_504;
   std::complex<double> tmp_505;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_505 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_504 += tmp_505;
   result += (-4) * tmp_504;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_506;
   std::complex<double> tmp_507;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_507 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPR(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1);
   }
   tmp_506 += tmp_507;
   result += (-0.5) * tmp_506;
   std::complex<double> tmp_508;
   std::complex<double> tmp_509;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_509 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPR(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2);
   }
   tmp_508 += tmp_509;
   result += (-0.5) * tmp_508;
   std::complex<double> tmp_510;
   std::complex<double> tmp_511;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_511 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPR(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2);
   }
   tmp_510 += tmp_511;
   result += (-0.5) * tmp_510;
   std::complex<double> tmp_512;
   std::complex<double> tmp_513;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_513 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_512 += tmp_513;
   result += (-1) * tmp_512;
   std::complex<double> tmp_514;
   std::complex<double> tmp_515;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_515 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPL(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2);
   }
   tmp_514 += tmp_515;
   result += (-1) * tmp_514;
   std::complex<double> tmp_516;
   std::complex<double> tmp_517;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_517 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_516 += tmp_517;
   result += (-1) * tmp_516;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_518;
   std::complex<double> tmp_519;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_519 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPL(gO1,gI1);
   }
   tmp_518 += tmp_519;
   result += (-0.5) * tmp_518;
   std::complex<double> tmp_520;
   std::complex<double> tmp_521;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_521 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPL(gO1,gI2);
   }
   tmp_520 += tmp_521;
   result += (-0.5) * tmp_520;
   std::complex<double> tmp_522;
   std::complex<double> tmp_523;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_523 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPL(gO1,gI2);
   }
   tmp_522 += tmp_523;
   result += (-0.5) * tmp_522;
   std::complex<double> tmp_524;
   std::complex<double> tmp_525;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_525 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_524 += tmp_525;
   result += (-1) * tmp_524;
   std::complex<double> tmp_526;
   std::complex<double> tmp_527;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_527 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPR(gO1,gI2);
   }
   tmp_526 += tmp_527;
   result += (-1) * tmp_526;
   std::complex<double> tmp_528;
   std::complex<double> tmp_529;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_529 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_528 += tmp_529;
   result += (-1) * tmp_528;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_530;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_530 += B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_530;
   std::complex<double> tmp_531;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_531 += B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_531;
   std::complex<double> tmp_532;
   std::complex<double> tmp_533;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_533 += (-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))
         *CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_532 += tmp_533;
   result += (-4) * tmp_532;
   std::complex<double> tmp_534;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_534 += B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_534;
   std::complex<double> tmp_535;
   std::complex<double> tmp_536;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_536 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_535 += tmp_536;
   result += (-4) * tmp_535;
   std::complex<double> tmp_537;
   std::complex<double> tmp_538;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_538 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_537 += tmp_538;
   result += (-4) * tmp_537;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_539;
   std::complex<double> tmp_540;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_540 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1);
   }
   tmp_539 += tmp_540;
   result += (-0.5) * tmp_539;
   std::complex<double> tmp_541;
   std::complex<double> tmp_542;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_542 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2);
   }
   tmp_541 += tmp_542;
   result += (-0.5) * tmp_541;
   std::complex<double> tmp_543;
   std::complex<double> tmp_544;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_544 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2);
   }
   tmp_543 += tmp_544;
   result += (-0.5) * tmp_543;
   std::complex<double> tmp_545;
   std::complex<double> tmp_546;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_546 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_545 += tmp_546;
   result += (-1) * tmp_545;
   std::complex<double> tmp_547;
   std::complex<double> tmp_548;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_548 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2);
   }
   tmp_547 += tmp_548;
   result += (-1) * tmp_547;
   std::complex<double> tmp_549;
   std::complex<double> tmp_550;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_550 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_549 += tmp_550;
   result += (-1) * tmp_549;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_551;
   std::complex<double> tmp_552;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_552 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPL(gO1,gI1);
   }
   tmp_551 += tmp_552;
   result += (-0.5) * tmp_551;
   std::complex<double> tmp_553;
   std::complex<double> tmp_554;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_554 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPL(gO1,gI2);
   }
   tmp_553 += tmp_554;
   result += (-0.5) * tmp_553;
   std::complex<double> tmp_555;
   std::complex<double> tmp_556;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_556 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPL(gO1,gI2);
   }
   tmp_555 += tmp_556;
   result += (-0.5) * tmp_555;
   std::complex<double> tmp_557;
   std::complex<double> tmp_558;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_558 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_557 += tmp_558;
   result += (-1) * tmp_557;
   std::complex<double> tmp_559;
   std::complex<double> tmp_560;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_560 += (0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPR(gO1,gI2);
   }
   tmp_559 += tmp_560;
   result += (-1) * tmp_559;
   std::complex<double> tmp_561;
   std::complex<double> tmp_562;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_562 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_561 += tmp_562;
   result += (-1) * tmp_561;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh() const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CphhAhAh();
   result += A0(MVWp)*CphhbargWpCgWpC();
   result += A0(MVWp)*CphhbargWpgWp();
   result += A0(MVZ)*CphhbargZgZ();
   result += -(A0(MHp)*CphhconjHpHp());
   result += -0.5*A0(Mhh)*Cphhhhhh();
   std::complex<double> tmp_563;
   std::complex<double> tmp_564;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_564 += A0(MCha(gI1))*(CphhbarChaChaPL(gI1,gI1) + CphhbarChaChaPR(
         gI1,gI1))*MCha(gI1);
   }
   tmp_563 += tmp_564;
   result += (2) * tmp_563;
   std::complex<double> tmp_565;
   std::complex<double> tmp_566;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_566 += A0(MFd(gI1))*(CphhbarFdFdPL(gI1,gI1) + CphhbarFdFdPR(gI1,
         gI1))*MFd(gI1);
   }
   tmp_565 += tmp_566;
   result += (6) * tmp_565;
   std::complex<double> tmp_567;
   std::complex<double> tmp_568;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_568 += A0(MFe(gI1))*(CphhbarFeFePL(gI1,gI1) + CphhbarFeFePR(gI1,
         gI1))*MFe(gI1);
   }
   tmp_567 += tmp_568;
   result += (2) * tmp_567;
   std::complex<double> tmp_569;
   std::complex<double> tmp_570;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_570 += A0(MFu(gI1))*(CphhbarFuFuPL(gI1,gI1) + CphhbarFuFuPR(gI1,
         gI1))*MFu(gI1);
   }
   tmp_569 += tmp_570;
   result += (6) * tmp_569;
   std::complex<double> tmp_571;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_571 += A0(MChi(gI1))*(CphhChiChiPL(gI1,gI1) + CphhChiChiPR(gI1,gI1
         ))*MChi(gI1);
   }
   result += tmp_571;
   result += 4*CphhconjVWpVWp()*(A0(MVWp) - 0.5*Sqr(MVWp));
   result += 2*CphhVZVZ()*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}




double CLASSNAME::self_energy_hh_2loop() const
{
   const double mt = MFu(2);
   const double yt = Yu(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy += self_energy_higgs_2loop_at_at_sm(scale, mt, yt);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy += self_energy_higgs_2loop_at_as_sm(scale, mt, yt, gs);
   }

   return self_energy;
}





double CLASSNAME::self_energy_hh_3loop() const
{
   const double mt = MFu(2);
   const double mg = MGlu;
   const double yt = Yu(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
      self_energy += self_energy_higgs_3loop_gluino_split(scale, mt, yt, gs, mg);
   }

   return self_energy;
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

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(hh))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const double M_tree(get_mass_matrix_hh());
      const double p = old_Mhh;
      double self_energy = Re(self_energy_hh(p));
      if (pole_mass_loop_order > 1)
         self_energy += self_energy_hh_2loop();
      if (pole_mass_loop_order > 2)
         self_energy += self_energy_hh_3loop();
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(Mhh) = SignedAbsSqrt(mass_sqr);

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(SplitMSSM_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(SplitMSSM_info::hh);
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
      problems.flag_bad_mass(SplitMSSM_info::Fd, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
      qcd_2l = -0.006995771808874528*Power(g3,4) -
         0.004518101565212638*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(2))) -
         0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(currentScale,2)/Sqr(
         MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = Power(g3,6)*(-0.0017408026847411467 -
         0.000984413176263005*Log(Sqr(currentScale)/Sqr(MFu(2))) -
         0.00003352082872926087*Power(Log(Sqr(currentScale)/Sqr(MFu(2))),3) -
         0.00029813221915266867*Sqr(Log(Power(currentScale,2)/Sqr(MFu(2)))));
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
      problems.flag_bad_mass(SplitMSSM_info::Fu, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
      problems.flag_bad_mass(SplitMSSM_info::Fe, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
         problems.flag_bad_mass(SplitMSSM_info::Chi,
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
      problems.flag_bad_mass(SplitMSSM_info::Cha, eigenvalue_error >
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

void CLASSNAME::calculate_MVWp_pole()
{
   if (!force_output && problems.is_tachyon(VWp))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWp));
   const double p = MVWp;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWp);

   PHYSICAL(MVWp) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_pole(double p)
{
   if (!force_output && problems.is_tachyon(VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(MVWp) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWp);

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
      qcd_2l = -0.005855107113909601*Power(g3,4) -
         0.0028071045227652486*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(currentScale,2)
         /Sqr(MFu(idx))));
   }

   if (get_thresholds() > 2) {
      qcd_3l = -0.0013067805969741943*Power(g3,6) -
         0.0004114970933517977*Power(g3,6)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 5.078913443827405e-6*Power(g3,6)*Power(Log(Sqr(currentScale)/Sqr(
         MFu(idx))),3) - 0.00007466002762426286*Power(g3,6)*Sqr(Log(Power(
         currentScale,2)/Sqr(MFu(idx))));
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

double CLASSNAME::calculate_MVWp_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(VWp);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}


std::ostream& operator<<(std::ostream& ostr, const SplitMSSM_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
