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
 * @file MRSSMEFTHiggs_two_scale_ewsb_solver.cpp
 *
 * @brief implementation of EWSB solver for two-scale iteration
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#include "MRSSMEFTHiggs_two_scale_ewsb_solver.hpp"
#include "MRSSMEFTHiggs_mass_eigenstates.hpp"
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
#define CLASSNAME MRSSMEFTHiggs_ewsb_solver<Two_scale>

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_iteratively(MRSSMEFTHiggs_mass_eigenstates& model_to_solve)
{
   auto model = model_to_solve;
   model.set_ewsb_loop_order(loop_order);

   auto ewsb_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double mHd2 = ewsb_pars(0);
      const double mHu2 = ewsb_pars(1);
      const double vS = ewsb_pars(2);
      const double vT = ewsb_pars(3);

      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);
      model.set_vS(vS);
      model.set_vT(vT);


      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->ewsb_step(model);
   };

   auto tadpole_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double mHd2 = ewsb_pars(0);
      const double mHu2 = ewsb_pars(1);
      const double vS = ewsb_pars(2);
      const double vT = ewsb_pars(3);

      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);
      model.set_vS(vS);
      model.set_vT(vT);


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
   MRSSMEFTHiggs_mass_eigenstates& model_to_solve, EWSB_solver* solver, const EWSB_vector_t& x_init)
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
void CLASSNAME::set_ewsb_solution(MRSSMEFTHiggs_mass_eigenstates& model, const EWSB_solver* solver)
{
   const auto solution = solver->get_solution();

   const double mHd2 = solution(0);
   const double mHu2 = solution(1);
   const double vS = solution(2);
   const double vT = solution(3);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);
   model.set_vS(vS);
   model.set_vT(vT);


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
void CLASSNAME::set_best_ewsb_solution(MRSSMEFTHiggs_mass_eigenstates& model, It first, It last)
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

int CLASSNAME::solve_iteratively_at(MRSSMEFTHiggs_mass_eigenstates& model_to_solve, int l)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(loop_order);
   loop_order = l;

   return solve_iteratively(model_to_solve);
}

int CLASSNAME::solve(MRSSMEFTHiggs_mass_eigenstates& model_to_solve)
{
   if (loop_order == 0) {
      return solve_tree_level(model_to_solve);
   }
   return solve_iteratively_at(model_to_solve, loop_order);
}

int CLASSNAME::solve_tree_level(MRSSMEFTHiggs_mass_eigenstates& model)
{
   int error = EWSB_solver::SUCCESS;

   const auto BMu = MODELPARAMETER(BMu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto mS2 = MODELPARAMETER(mS2);
   const auto mT2 = MODELPARAMETER(mT2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vd = MODELPARAMETER(vd);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto MuD = MODELPARAMETER(MuD);
   const auto MuU = MODELPARAMETER(MuU);
   const auto Mu = MODELPARAMETER(Mu);

   double vT;
   double vS;
   double mHd2;
   double mHu2;

   vT = Re((0.2*(-20*g2*MDWBT*AbsSqr(LamSD)*Quad(vd) - 5.477225575051661*g1*LamTD*
      MDBS*Conj(LamSD)*Quad(vd) - 5.477225575051661*g1*LamSD*MDBS*Conj(LamTD)*Quad
      (vd) - 10*MuD*AbsSqr(LamSD)*Conj(LamTD)*Quad(vd) - 5.477225575051661*g1*
      LamTD*Conj(LamSD)*Conj(MDBS)*Quad(vd) - 5.477225575051661*g1*LamSD*Conj(
      LamTD)*Conj(MDBS)*Quad(vd) - 20*g2*AbsSqr(LamSD)*Conj(MDWBT)*Quad(vd) - 10*
      LamTD*AbsSqr(LamSD)*Conj(MuD)*Quad(vd) + 20*g2*MDWBT*AbsSqr(LamSU)*Quad(vu)
      - 5.477225575051661*g1*LamTU*MDBS*Conj(LamSU)*Quad(vu) - 5.477225575051661*
      g1*LamSU*MDBS*Conj(LamTU)*Quad(vu) + 10*MuU*AbsSqr(LamSU)*Conj(LamTU)*Quad(
      vu) - 5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(MDBS)*Quad(vu) -
      5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Quad(vu) + 20*g2*AbsSqr(
      LamSU)*Conj(MDWBT)*Quad(vu) + 10*LamTU*AbsSqr(LamSU)*Conj(MuU)*Quad(vu) + 10
      *Conj(LamTD)*Conj(MuD)*Quad(vd)*Sqr(LamSD) - 10*Conj(LamTU)*Conj(MuU)*Quad(
      vu)*Sqr(LamSU) - 40*g2*MDWBT*mS2*Sqr(vd) - 40*mS2*MuD*Conj(LamTD)*Sqr(vd) -
      40*g2*mS2*Conj(MDWBT)*Sqr(vd) - 40*LamTD*mS2*Conj(MuD)*Sqr(vd) - 160*g2*
      MDWBT*Sqr(MDBS)*Sqr(vd) - 160*MuD*Conj(LamTD)*Sqr(MDBS)*Sqr(vd) - 160*g2*
      Conj(MDWBT)*Sqr(MDBS)*Sqr(vd) - 160*LamTD*Conj(MuD)*Sqr(MDBS)*Sqr(vd) + 40*
      g2*MDWBT*mS2*Sqr(vu) + 40*mS2*MuU*Conj(LamTU)*Sqr(vu) + 40*g2*mS2*Conj(MDWBT
      )*Sqr(vu) + 40*LamTU*mS2*Conj(MuU)*Sqr(vu) + 160*g2*MDWBT*Sqr(MDBS)*Sqr(vu)
      + 160*MuU*Conj(LamTU)*Sqr(MDBS)*Sqr(vu) + 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(
      vu) + 160*LamTU*Conj(MuU)*Sqr(MDBS)*Sqr(vu) + 20*g2*MDWBT*AbsSqr(LamSD)*Sqr(
      vd)*Sqr(vu) - 20*g2*MDWBT*AbsSqr(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*
      g1*LamTD*MDBS*Conj(LamSD)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*MDBS*
      Conj(LamSU)*Sqr(vd)*Sqr(vu) - 10*LamTU*MuD*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*
      Sqr(vu) + 10*LamTD*MuU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSD*MDBS*Conj(LamTD)*Sqr(vd)*Sqr(vu) - 20*MuD*AbsSqr(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 10*LamSD*MuU*Conj(LamSU)*Conj(LamTD)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSU*MDBS*Conj(LamTU)*Sqr(vd)*Sqr(vu
      ) + 20*MuU*AbsSqr(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) - 10*LamSU*MuD*Conj(
      LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTD*Conj(LamSD)*
      Conj(MDBS)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(
      MDBS)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Sqr(vd)*
      Sqr(vu) + 20*g2*AbsSqr(LamSD)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*g2*AbsSqr(
      LamSU)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*LamTD*AbsSqr(LamSU)*Conj(MuD)*Sqr(vd
      )*Sqr(vu) - 10*LamSD*LamTU*Conj(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*LamSD*
      LamSU*Conj(LamTU)*Conj(MuD)*Sqr(vd)*Sqr(vu) + 20*LamTU*AbsSqr(LamSD)*Conj(
      MuU)*Sqr(vd)*Sqr(vu) + 10*LamSU*LamTD*Conj(LamSD)*Conj(MuU)*Sqr(vd)*Sqr(vu)
      + 10*LamSD*LamSU*Conj(LamTD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuD*Quad(
      vd)*Sqr(Conj(LamSD)) - 10*LamTU*MuU*Quad(vu)*Sqr(Conj(LamSU))))/(32*mS2*mT2
      + 2*AbsSqr(LamSD)*AbsSqr(LamTD)*Quad(vd) + 2*AbsSqr(LamSU)*AbsSqr(LamTU)*
      Quad(vu) + 128*mT2*Sqr(MDBS) + 128*mS2*Sqr(MDWBT) + 512*Sqr(MDBS)*Sqr(MDWBT)
      + 16*mT2*AbsSqr(LamSD)*Sqr(vd) + 8*mS2*AbsSqr(LamTD)*Sqr(vd) + 32*AbsSqr(
      LamTD)*Sqr(MDBS)*Sqr(vd) + 64*AbsSqr(LamSD)*Sqr(MDWBT)*Sqr(vd) + 16*mT2*
      AbsSqr(LamSU)*Sqr(vu) + 8*mS2*AbsSqr(LamTU)*Sqr(vu) + 32*AbsSqr(LamTU)*Sqr(
      MDBS)*Sqr(vu) + 64*AbsSqr(LamSU)*Sqr(MDWBT)*Sqr(vu) + 4*AbsSqr(LamSU)*AbsSqr
      (LamTD)*Sqr(vd)*Sqr(vu) + 4*AbsSqr(LamSD)*AbsSqr(LamTU)*Sqr(vd)*Sqr(vu) + 2*
      LamTD*LamTU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 2*LamSD*LamTU*Conj(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 2*LamSU*LamTD*Conj(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) + 2*LamSD*LamSU*Conj(LamTD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) -
      Quad(vd)*Sqr(LamTD)*Sqr(Conj(LamSD)) - Quad(vu)*Sqr(LamTU)*Sqr(Conj(LamSU))
      - Quad(vd)*Sqr(LamSD)*Sqr(Conj(LamTD)) - Quad(vu)*Sqr(LamSU)*Sqr(Conj(LamTU)
      )));
   vS = Re((-1.4142135623730951*(4*mT2*vT + 16*vT*Sqr(MDWBT) + g2*MDWBT*Sqr(vd) +
      vT*AbsSqr(LamTD)*Sqr(vd) + MuD*Conj(LamTD)*Sqr(vd) + g2*Conj(MDWBT)*Sqr(vd)
      + LamTD*Conj(MuD)*Sqr(vd) - g2*MDWBT*Sqr(vu) + vT*AbsSqr(LamTU)*Sqr(vu) -
      MuU*Conj(LamTU)*Sqr(vu) - g2*Conj(MDWBT)*Sqr(vu) - LamTU*Conj(MuU)*Sqr(vu)))
      /(LamTD*Conj(LamSD)*Sqr(vd) + LamSD*Conj(LamTD)*Sqr(vd) - LamTU*Conj(LamSU)*
      Sqr(vu) - LamSU*Conj(LamTU)*Sqr(vu)));
   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*vd*
      AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS*
      Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) - 3*Cube(vd)*Sqr(g1) - 5*Cube(vd)*Sqr(g2) - 20*vd*AbsSqr(LamSD)*
      Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*
      Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40*vu
      *AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*vu*
      Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*vu*
      Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) - 3*Cube(vu)*Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) +
      5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu*AbsSqr(LamTU)*Sqr
      (vT)))/vu);

   
   const bool is_finite = IsFinite(vT) && IsFinite(vS) && IsFinite(mHd2) &&
      IsFinite(mHu2);

   if (is_finite) {
      model.set_vT(vT);
      model.set_vS(vS);
      model.set_mHd2(mHd2);
      model.set_mHu2(mHu2);
      model.get_problems().unflag_no_ewsb_tree_level();
   } else {
      error = EWSB_solver::FAIL;
      model.get_problems().flag_no_ewsb_tree_level();
   }
   return error;
}

CLASSNAME::EWSB_vector_t CLASSNAME::initial_guess(const MRSSMEFTHiggs_mass_eigenstates& model) const
{
   EWSB_vector_t x_init(EWSB_vector_t::Zero());

   const auto mHd2 = MODELPARAMETER(mHd2);
   const auto mHu2 = MODELPARAMETER(mHu2);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   x_init[0] = mHd2;
   x_init[1] = mHu2;
   x_init[2] = vS;
   x_init[3] = vT;


   return x_init;
}

CLASSNAME::EWSB_vector_t CLASSNAME::tadpole_equations(const MRSSMEFTHiggs_mass_eigenstates& model) const
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
CLASSNAME::EWSB_vector_t CLASSNAME::ewsb_step(const MRSSMEFTHiggs_mass_eigenstates& model) const
{
   std::array<double, number_of_ewsb_equations> tadpole{};
   EWSB_vector_t ewsb_parameters(EWSB_vector_t::Zero());

   if (loop_order > 0) {
   tadpole[0] += Re(model.tadpole_hh_1loop(0));
   tadpole[1] += Re(model.tadpole_hh_1loop(1));
   tadpole[2] += Re(model.tadpole_hh_1loop(3));
   tadpole[3] += Re(model.tadpole_hh_1loop(2));

      if (loop_order > 1) {

      }
   }

   const auto BMu = MODELPARAMETER(BMu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto mS2 = MODELPARAMETER(mS2);
   const auto mT2 = MODELPARAMETER(mT2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vd = MODELPARAMETER(vd);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto MuD = MODELPARAMETER(MuD);
   const auto MuU = MODELPARAMETER(MuU);
   const auto Mu = MODELPARAMETER(Mu);
   double vT;
   double vS;
   double mHd2;
   double mHu2;

   vT = Re((0.2*(160*mS2*tadpole[2] - 20*g2*MDWBT*AbsSqr(LamSD)*Quad(vd) -
      5.477225575051661*g1*LamTD*MDBS*Conj(LamSD)*Quad(vd) - 5.477225575051661*g1*
      LamSD*MDBS*Conj(LamTD)*Quad(vd) - 10*MuD*AbsSqr(LamSD)*Conj(LamTD)*Quad(vd)
      - 5.477225575051661*g1*LamTD*Conj(LamSD)*Conj(MDBS)*Quad(vd) -
      5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*Quad(vd) - 20*g2*AbsSqr(
      LamSD)*Conj(MDWBT)*Quad(vd) - 10*LamTD*AbsSqr(LamSD)*Conj(MuD)*Quad(vd) + 20
      *g2*MDWBT*AbsSqr(LamSU)*Quad(vu) - 5.477225575051661*g1*LamTU*MDBS*Conj(
      LamSU)*Quad(vu) - 5.477225575051661*g1*LamSU*MDBS*Conj(LamTU)*Quad(vu) + 10*
      MuU*AbsSqr(LamSU)*Conj(LamTU)*Quad(vu) - 5.477225575051661*g1*LamTU*Conj(
      LamSU)*Conj(MDBS)*Quad(vu) - 5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(
      MDBS)*Quad(vu) + 20*g2*AbsSqr(LamSU)*Conj(MDWBT)*Quad(vu) + 10*LamTU*AbsSqr(
      LamSU)*Conj(MuU)*Quad(vu) + 10*Conj(LamTD)*Conj(MuD)*Quad(vd)*Sqr(LamSD) -
      10*Conj(LamTU)*Conj(MuU)*Quad(vu)*Sqr(LamSU) + 640*tadpole[2]*Sqr(MDBS) - 40
      *g2*MDWBT*mS2*Sqr(vd) - 40*mS2*MuD*Conj(LamTD)*Sqr(vd) - 40*g2*mS2*Conj(
      MDWBT)*Sqr(vd) - 40*LamTD*mS2*Conj(MuD)*Sqr(vd) + 80*AbsSqr(LamSD)*tadpole[2
      ]*Sqr(vd) - 28.284271247461902*LamTD*Conj(LamSD)*tadpole[3]*Sqr(vd) -
      28.284271247461902*LamSD*Conj(LamTD)*tadpole[3]*Sqr(vd) - 160*g2*MDWBT*Sqr(
      MDBS)*Sqr(vd) - 160*MuD*Conj(LamTD)*Sqr(MDBS)*Sqr(vd) - 160*g2*Conj(MDWBT)*
      Sqr(MDBS)*Sqr(vd) - 160*LamTD*Conj(MuD)*Sqr(MDBS)*Sqr(vd) + 40*g2*MDWBT*mS2*
      Sqr(vu) + 40*mS2*MuU*Conj(LamTU)*Sqr(vu) + 40*g2*mS2*Conj(MDWBT)*Sqr(vu) +
      40*LamTU*mS2*Conj(MuU)*Sqr(vu) + 80*AbsSqr(LamSU)*tadpole[2]*Sqr(vu) +
      28.284271247461902*LamTU*Conj(LamSU)*tadpole[3]*Sqr(vu) + 28.284271247461902
      *LamSU*Conj(LamTU)*tadpole[3]*Sqr(vu) + 160*g2*MDWBT*Sqr(MDBS)*Sqr(vu) + 160
      *MuU*Conj(LamTU)*Sqr(MDBS)*Sqr(vu) + 160*g2*Conj(MDWBT)*Sqr(MDBS)*Sqr(vu) +
      160*LamTU*Conj(MuU)*Sqr(MDBS)*Sqr(vu) + 20*g2*MDWBT*AbsSqr(LamSD)*Sqr(vd)*
      Sqr(vu) - 20*g2*MDWBT*AbsSqr(LamSU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*
      LamTD*MDBS*Conj(LamSD)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*MDBS*
      Conj(LamSU)*Sqr(vd)*Sqr(vu) - 10*LamTU*MuD*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*
      Sqr(vu) + 10*LamTD*MuU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) +
      5.477225575051661*g1*LamSD*MDBS*Conj(LamTD)*Sqr(vd)*Sqr(vu) - 20*MuD*AbsSqr(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 10*LamSD*MuU*Conj(LamSU)*Conj(LamTD)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSU*MDBS*Conj(LamTU)*Sqr(vd)*Sqr(vu
      ) + 20*MuU*AbsSqr(LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) - 10*LamSU*MuD*Conj(
      LamSD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTD*Conj(LamSD)*
      Conj(MDBS)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamTU*Conj(LamSU)*Conj(
      MDBS)*Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSD*Conj(LamTD)*Conj(MDBS)*
      Sqr(vd)*Sqr(vu) + 5.477225575051661*g1*LamSU*Conj(LamTU)*Conj(MDBS)*Sqr(vd)*
      Sqr(vu) + 20*g2*AbsSqr(LamSD)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*g2*AbsSqr(
      LamSU)*Conj(MDWBT)*Sqr(vd)*Sqr(vu) - 20*LamTD*AbsSqr(LamSU)*Conj(MuD)*Sqr(vd
      )*Sqr(vu) - 10*LamSD*LamTU*Conj(LamSU)*Conj(MuD)*Sqr(vd)*Sqr(vu) - 10*LamSD*
      LamSU*Conj(LamTU)*Conj(MuD)*Sqr(vd)*Sqr(vu) + 20*LamTU*AbsSqr(LamSD)*Conj(
      MuU)*Sqr(vd)*Sqr(vu) + 10*LamSU*LamTD*Conj(LamSD)*Conj(MuU)*Sqr(vd)*Sqr(vu)
      + 10*LamSD*LamSU*Conj(LamTD)*Conj(MuU)*Sqr(vd)*Sqr(vu) + 10*LamTD*MuD*Quad(
      vd)*Sqr(Conj(LamSD)) - 10*LamTU*MuU*Quad(vu)*Sqr(Conj(LamSU))))/(32*mS2*mT2
      + 2*AbsSqr(LamSD)*AbsSqr(LamTD)*Quad(vd) + 2*AbsSqr(LamSU)*AbsSqr(LamTU)*
      Quad(vu) + 128*mT2*Sqr(MDBS) + 128*mS2*Sqr(MDWBT) + 512*Sqr(MDBS)*Sqr(MDWBT)
      + 16*mT2*AbsSqr(LamSD)*Sqr(vd) + 8*mS2*AbsSqr(LamTD)*Sqr(vd) + 32*AbsSqr(
      LamTD)*Sqr(MDBS)*Sqr(vd) + 64*AbsSqr(LamSD)*Sqr(MDWBT)*Sqr(vd) + 16*mT2*
      AbsSqr(LamSU)*Sqr(vu) + 8*mS2*AbsSqr(LamTU)*Sqr(vu) + 32*AbsSqr(LamTU)*Sqr(
      MDBS)*Sqr(vu) + 64*AbsSqr(LamSU)*Sqr(MDWBT)*Sqr(vu) + 4*AbsSqr(LamSU)*AbsSqr
      (LamTD)*Sqr(vd)*Sqr(vu) + 4*AbsSqr(LamSD)*AbsSqr(LamTU)*Sqr(vd)*Sqr(vu) + 2*
      LamTD*LamTU*Conj(LamSD)*Conj(LamSU)*Sqr(vd)*Sqr(vu) + 2*LamSD*LamTU*Conj(
      LamSU)*Conj(LamTD)*Sqr(vd)*Sqr(vu) + 2*LamSU*LamTD*Conj(LamSD)*Conj(LamTU)*
      Sqr(vd)*Sqr(vu) + 2*LamSD*LamSU*Conj(LamTD)*Conj(LamTU)*Sqr(vd)*Sqr(vu) -
      Quad(vd)*Sqr(LamTD)*Sqr(Conj(LamSD)) - Quad(vu)*Sqr(LamTU)*Sqr(Conj(LamSU))
      - Quad(vd)*Sqr(LamSD)*Sqr(Conj(LamTD)) - Quad(vu)*Sqr(LamSU)*Sqr(Conj(LamTU)
      )));
   vS = Re((-1.4142135623730951*(4*mT2*vT - 4*tadpole[2] + 16*vT*Sqr(MDWBT) + g2*
      MDWBT*Sqr(vd) + vT*AbsSqr(LamTD)*Sqr(vd) + MuD*Conj(LamTD)*Sqr(vd) + g2*Conj
      (MDWBT)*Sqr(vd) + LamTD*Conj(MuD)*Sqr(vd) - g2*MDWBT*Sqr(vu) + vT*AbsSqr(
      LamTU)*Sqr(vu) - MuU*Conj(LamTU)*Sqr(vu) - g2*Conj(MDWBT)*Sqr(vu) - LamTU*
      Conj(MuU)*Sqr(vu)))/(LamTD*Conj(LamSD)*Sqr(vd) + LamSD*Conj(LamTD)*Sqr(vd) -
      LamTU*Conj(LamSU)*Sqr(vu) - LamSU*Conj(LamTU)*Sqr(vu)));
   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*vd*
      AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS*
      Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) + 40*tadpole[0] - 3*Cube(vd)*Sqr(g1) - 5*Cube(vd)*Sqr(g2) - 20*vd*
      AbsSqr(LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) +
      5*vd*Sqr(g2)*Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40*vu
      *AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*vu*
      Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*vu*
      Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) + 40*tadpole[1] - 3*Cube(vu)*Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*
      Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu
      *AbsSqr(LamTU)*Sqr(vT)))/vu);

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

} // namespace flexiblesusy
