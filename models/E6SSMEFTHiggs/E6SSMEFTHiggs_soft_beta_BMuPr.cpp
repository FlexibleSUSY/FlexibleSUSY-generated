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


#include "E6SSMEFTHiggs_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of BMuPr.
 *
 * @return 1-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_BMuPr_1_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(0.2*(-3*BMuPr*Sqr(g1) + 6*MassB*MuPr*Sqr(g1) - 15*BMuPr*Sqr(
      g2) + 30*MassWB*MuPr*Sqr(g2) - 2*BMuPr*Sqr(gN) + 4*MassBp*MuPr*Sqr(gN)));


   return oneLoop * beta_BMuPr;
}

/**
 * Calculates the 2-loop beta function of BMuPr.
 *
 * @return 2-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_BMuPr_2_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(-0.06*(-99*BMuPr*Quad(g1) + 396*MassB*MuPr*Quad(g1) - 275*
      BMuPr*Quad(g2) + 1100*MassWB*MuPr*Quad(g2) - 64*BMuPr*Quad(gN) + 256*
      MassBp*MuPr*Quad(gN) - 30*BMuPr*Sqr(g1)*Sqr(g2) + 60*MassB*MuPr*Sqr(g1)*
      Sqr(g2) + 60*MassWB*MuPr*Sqr(g1)*Sqr(g2) - 12*BMuPr*Sqr(g1)*Sqr(gN) + 24*
      MassB*MuPr*Sqr(g1)*Sqr(gN) + 24*MassBp*MuPr*Sqr(g1)*Sqr(gN) - 20*BMuPr*
      Sqr(g2)*Sqr(gN) + 40*MassBp*MuPr*Sqr(g2)*Sqr(gN) + 40*MassWB*MuPr*Sqr(g2)
      *Sqr(gN)));


   return twoLoop * beta_BMuPr;
}

/**
 * Calculates the 3-loop beta function of BMuPr.
 *
 * @return 3-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_BMuPr_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return threeLoop * beta_BMuPr;
}

/**
 * Calculates the 4-loop beta function of BMuPr.
 *
 * @return 4-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_BMuPr_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return fourLoop * beta_BMuPr;
}

/**
 * Calculates the 5-loop beta function of BMuPr.
 *
 * @return 5-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_BMuPr_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return fiveLoop * beta_BMuPr;
}

} // namespace flexiblesusy
