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


/**
 * @file NUHMSSM_two_scale_ewsb_solver.cpp
 *
 * @brief implementation of EWSB solver for two-scale iteration
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "NUHMSSM_two_scale_ewsb_solver.hpp"
#include "NUHMSSM_mass_eigenstates.hpp"
#include "logger.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "raii.hpp"
#include "wrappers.hpp"

#include <memory>

namespace flexiblesusy {

#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) INPUT(parameter)
#define INPUTPARAMETER(parameter) LOCALINPUT(parameter)
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.parameter()
#define PHASE(p) MODELPARAMETER(p)
#define LowEnergyConstant(p) Electroweak_constants::p
#define CLASSNAME NUHMSSM_ewsb_solver<Two_scale>

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_iteratively(NUHMSSM_mass_eigenstates& model_to_solve)
{
   auto model = model_to_solve;
   model.set_ewsb_loop_order(loop_order);

   auto ewsb_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double mHd2 = ewsb_pars(0);
      const double mHu2 = ewsb_pars(1);

      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);


      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->ewsb_step(model);
   };

   auto tadpole_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double mHd2 = ewsb_pars(0);
      const double mHu2 = ewsb_pars(1);

      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);


      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->tadpole_equations(model);
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations> >(ewsb_stepper, number_of_iterations, fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations>(precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, number_of_iterations, precision, Root_finder<number_of_ewsb_equations>::GSLHybridS)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, number_of_iterations, precision, Root_finder<number_of_ewsb_equations>::GSLBroyden))
   };

   const auto x_init(initial_guess(model_to_solve));

   VERBOSE_MSG("\t\tSolving EWSB equations ...");
   VERBOSE_MSG("\t\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (auto& solver: solvers) {
      VERBOSE_MSG("\t\t\tStarting EWSB iteration using " << solver->name());
      status = solve_iteratively_with(model_to_solve, solver.get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\t\t\t" << solver->name() << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\t\t\t" << solver->name() << " could not find a solution!"
                 " (requested precision: " << precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      model_to_solve.get_problems().unflag_no_ewsb();
   } else {
      set_best_ewsb_solution(model_to_solve, std::begin(solvers), std::end(solvers));
      model_to_solve.get_problems().flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\t\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param model_to_solve model to solve EWSB conditions for
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_iteratively_with(
   NUHMSSM_mass_eigenstates& model_to_solve, EWSB_solver* solver, const EWSB_vector_t& x_init)
{
   const int status = solver->solve(x_init);

   if (status == EWSB_solver::SUCCESS)
      set_ewsb_solution(model_to_solve, solver);

   return status;
}

/**
 * Sets EWSB output parameters from given solver.
 *
 * @param model model to set EWSB output parameters in
 * @param solver solver
 */
void CLASSNAME::set_ewsb_solution(NUHMSSM_mass_eigenstates& model, const EWSB_solver* solver)
{
   const auto solution = solver->get_solution();

   const double mHd2 = solution(0);
   const double mHu2 = solution(1);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);


   model.calculate_DRbar_masses();
}

/**
 * Sets EWSB output parameters from the solver from the range [first,
 * last), which minimizes the tadpole equations at most.
 *
 * @param model model to set EWSB output parameters in
 * @param first iterator to first solver
 * @param last iterator to last solver
 */
template <typename It>
void CLASSNAME::set_best_ewsb_solution(NUHMSSM_mass_eigenstates& model, It first, It last)
{
   auto ma(model), mb(model);

   const auto best_solver =
      std::min_element(first, last,
                       [this, &ma, &mb](const std::unique_ptr<EWSB_solver>& a, const std::unique_ptr<EWSB_solver>& b) {
                          this->set_ewsb_solution(ma, a.get());
                          this->set_ewsb_solution(mb, b.get());
                          return Total(Abs(Re(ma.tadpole_equations()))) < Total(Abs(Re(mb.tadpole_equations())));
                       });

   VERBOSE_MSG("\t\tUsing best solution from " << (*best_solver)->name());

   set_ewsb_solution(model, best_solver->get());
}

int CLASSNAME::solve_iteratively_at(NUHMSSM_mass_eigenstates& model_to_solve, int l)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(loop_order);
   loop_order = l;

   return solve_iteratively(model_to_solve);
}

int CLASSNAME::solve(NUHMSSM_mass_eigenstates& model_to_solve)
{
   if (loop_order == 0) {
      return solve_tree_level(model_to_solve);
   }
   return solve_iteratively_at(model_to_solve, loop_order);
}

int CLASSNAME::solve_tree_level(NUHMSSM_mass_eigenstates& model)
{
   int error = EWSB_solver::SUCCESS;

   const auto BMu = MODELPARAMETER(BMu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);

   double mHd2;
   double mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) - 3*Cube(vd)*
      Sqr(g1) - 5*Cube(vd)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu)))
      /vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) - 3*Cube(vu)*
      Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd)))
      /vu);

   
   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);

   if (is_finite) {
      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);
      model.get_problems().unflag_no_ewsb_tree_level();
   } else {
      error = EWSB_solver::FAIL;
      model.get_problems().flag_no_ewsb_tree_level();
   }
   return error;
}

CLASSNAME::EWSB_vector_t CLASSNAME::initial_guess(const NUHMSSM_mass_eigenstates& model) const
{
   EWSB_vector_t x_init(EWSB_vector_t::Zero());

   const auto mHd2 = MODELPARAMETER(mHd2);
   const auto mHu2 = MODELPARAMETER(mHu2);
   x_init[0] = mHd2;
   x_init[1] = mHu2;


   return x_init;
}

CLASSNAME::EWSB_vector_t CLASSNAME::tadpole_equations(const NUHMSSM_mass_eigenstates& model) const
{
   return model.tadpole_equations();
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @param model model to use for calculating the EWSB output parameters
 *
 * @return new set of EWSB output parameters
 */
CLASSNAME::EWSB_vector_t CLASSNAME::ewsb_step(const NUHMSSM_mass_eigenstates& model) const
{
   std::array<double, number_of_ewsb_equations> tadpole{};
   EWSB_vector_t ewsb_parameters(EWSB_vector_t::Zero());

   if (loop_order > 0) {
   tadpole[0] += Re(model.tadpole_hh_1loop(0));
   tadpole[1] += Re(model.tadpole_hh_1loop(1));

      if (loop_order > 1) {
   const auto tadpole_2l(model.tadpole_hh_2loop());
   tadpole[0] += tadpole_2l(0);
   tadpole[1] += tadpole_2l(1);

      }
   }

   const auto BMu = MODELPARAMETER(BMu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   double mHd2;
   double mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) + 40*tadpole[
      0] - 3*Cube(vd)*Sqr(g1) - 5*Cube(vd)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*
      Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) + 40*tadpole[
      1] - 3*Cube(vu)*Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*
      Sqr(g2)*Sqr(vd)))/vu);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mHd2;
   ewsb_parameters[1] = mHu2;


   return ewsb_parameters;
}

} // namespace flexiblesusy
