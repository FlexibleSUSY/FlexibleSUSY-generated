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


#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vS.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_1_loop(const Susy_traces& susy_traces) const
{
   const auto Qs = INPUT(Qs);


   double beta_vS;

   beta_vS = Re(2*vS*(-AbsSqr(Lambdax) + Sqr(gp)*Sqr(Qs)));


   return oneLoop * beta_vS;
}

/**
 * Calculates the 2-loop beta function of vS.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qs = INPUT(Qs);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vS;

   beta_vS = Re(-0.2*vS*(-30*traceYdAdjYd*AbsSqr(Lambdax) - 10*traceYeAdjYe*
      AbsSqr(Lambdax) - 30*traceYuAdjYu*AbsSqr(Lambdax) - 10*traceYvAdjYv*
      AbsSqr(Lambdax) + 15*Quad(gp)*Quad(Qs) + 6*AbsSqr(Lambdax)*Sqr(g1) + 30*
      AbsSqr(Lambdax)*Sqr(g2) + 20*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*AbsSqr
      (Lambdax)*Sqr(gp)*Sqr(QHu) + 45*Quad(gp)*Sqr(Qd)*Sqr(Qs) + 15*Quad(gp)*
      Sqr(Qe)*Sqr(Qs) + 10*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 10*Quad(gp)*Sqr(QHu)*Sqr
      (Qs) + 30*Quad(gp)*Sqr(Ql)*Sqr(Qs) + 90*Quad(gp)*Sqr(Qq)*Sqr(Qs) + 45*
      Quad(gp)*Sqr(Qs)*Sqr(Qu) + 15*Quad(gp)*Sqr(Qs)*Sqr(Qv) - 20*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)));


   return twoLoop * beta_vS;
}

/**
 * Calculates the 3-loop beta function of vS.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return threeLoop * beta_vS;
}

/**
 * Calculates the 4-loop beta function of vS.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return fourLoop * beta_vS;
}

/**
 * Calculates the 5-loop beta function of vS.
 *
 * @return 5-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return fiveLoop * beta_vS;
}

} // namespace flexiblesusy
