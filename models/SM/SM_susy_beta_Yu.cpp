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


#include "SM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (-0.05*Yu*(-60*traceYdAdjYd - 20*traceYeAdjYe - 60*traceYuAdjYu +
      17*Sqr(g1) + 45*Sqr(g2) + 160*Sqr(g3)) - 1.5*(Yu*Yd.adjoint()*Yd) + 1.5*(
      Yu*Yu.adjoint()*Yu)).real();


   return oneLoop * beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.0016666666666666668*Yu*(-4050*traceYdAdjYdYdAdjYd + 900*
      traceYdAdjYuYuAdjYd - 1350*traceYeAdjYeYeAdjYe - 4050*traceYuAdjYuYuAdjYu
       + 1187*Quad(g1) - 3450*Quad(g2) - 64800*Quad(g3) + 375*traceYdAdjYd*Sqr(
      g1) + 1125*traceYeAdjYe*Sqr(g1) + 1275*traceYuAdjYu*Sqr(g1) + 3375*
      traceYdAdjYd*Sqr(g2) + 1125*traceYeAdjYe*Sqr(g2) + 3375*traceYuAdjYu*Sqr(
      g2) - 270*Sqr(g1)*Sqr(g2) + 12000*traceYdAdjYd*Sqr(g3) + 12000*
      traceYuAdjYu*Sqr(g3) + 760*Sqr(g1)*Sqr(g3) + 5400*Sqr(g2)*Sqr(g3) + 900*
      Sqr(Lambdax)) + 0.0125*(300*traceYdAdjYd + 100*traceYeAdjYe + 300*
      traceYuAdjYu - 43*Sqr(g1) + 45*Sqr(g2) - 1280*Sqr(g3))*(Yu*Yd.adjoint()*
      Yd) + 0.0125*(-540*traceYdAdjYd - 180*traceYeAdjYe - 540*traceYuAdjYu -
      480*Lambdax + 223*Sqr(g1) + 675*Sqr(g2) + 1280*Sqr(g3))*(Yu*Yu.adjoint()*
      Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*
      Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*Yu.
      adjoint()*Yu*Yu.adjoint()*Yu)).real();


   return twoLoop * beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (16.09896843832796*(6.1494623322761255*PROJECTOR*Lambdax*Power5(Yu
      (2,2)) + 3.6401567229074585*PROJECTOR*Power7(Yu(2,2)) -
      1.5170193131852032*PROJECTOR*Cube(Yu(2,2))*Quad(g1) + 1.0553414347426693*
      PROJECTOR*Cube(Yu(2,2))*Quad(g2) + 22.59550733368146*PROJECTOR*Cube(Yu(2,
      2))*Quad(g3) + 2.821688161099187*PROJECTOR*Cube(Yu(2,2))*Quad(Yd(2,2)) +
      1.6072458368448965*PROJECTOR*Cube(Yu(2,2))*Quad(Ye(2,2)) -
      0.3944352102015494*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g1) -
      1.8922019827975904*PROJECTOR*Power5(Yu(2,2))*Sqr(g1) - 2.096407613275952*
      PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g2) - 6.184402459164058*PROJECTOR*
      Power5(Yu(2,2))*Sqr(g2) + 2.163462366965912*PROJECTOR*Cube(Yu(2,2))*Sqr(
      g1)*Sqr(g2) + 0.49692624907281824*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g3)
      - 9.752177638054057*PROJECTOR*Power5(Yu(2,2))*Sqr(g3) +
      1.1226836416869554*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(g3) +
      3.0045554008024813*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(g3) +
      0.058233544813220885*PROJECTOR*Cube(Yu(2,2))*Sqr(Lambdax) +
      2.8883838227357557*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(Yd(2,2)) +
      2.8689726411313488*PROJECTOR*Power5(Yu(2,2))*Sqr(Yd(2,2)) -
      0.49957993142419754*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(Yd(2,2)) -
      4.814159761919977*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(Yd(2,2)) -
      0.7122084216159787*PROJECTOR*Cube(Yu(2,2))*Sqr(g3)*Sqr(Yd(2,2)) +
      0.9317367170115342*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(Ye(2,2)) +
      0.6522157019080739*PROJECTOR*Power5(Yu(2,2))*Sqr(Ye(2,2)) -
      0.4242586654541776*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(Ye(2,2)) -
      1.9298448995321542*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(Ye(2,2)) +
      0.1552894528352557*PROJECTOR*Cube(Yu(2,2))*Sqr(g3)*Sqr(Ye(2,2)) +
      0.21740523396935796*PROJECTOR*Cube(Yu(2,2))*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      0.2795210151034602*PROJECTOR*Cube(Lambdax)*Yu(2,2) + 1.*PROJECTOR*Power6(
      g1)*Yu(2,2) + 10.549072334227892*PROJECTOR*Power6(g2)*Yu(2,2) -
      38.47142449015523*PROJECTOR*Power6(g3)*Yu(2,2) + 2.187826890843716*
      PROJECTOR*Power6(Yd(2,2))*Yu(2,2) + 0.4996388893047732*PROJECTOR*Power6(
      Ye(2,2))*Yu(2,2) - 0.08455510706879672*PROJECTOR*Lambdax*Quad(g1)*Yu(2,2)
      - 0.33193120543535903*PROJECTOR*Lambdax*Quad(g2)*Yu(2,2) +
      0.4658683585057671*PROJECTOR*Lambdax*Quad(Yd(2,2))*Yu(2,2) +
      0.4658683585057671*PROJECTOR*Lambdax*Quad(Ye(2,2))*Yu(2,2) -
      0.29462264646735065*PROJECTOR*Quad(g2)*Sqr(g1)*Yu(2,2) -
      0.9376834535939306*PROJECTOR*Quad(g3)*Sqr(g1)*Yu(2,2) -
      0.23043972650846906*PROJECTOR*Quad(Yd(2,2))*Sqr(g1)*Yu(2,2) -
      0.5080205796037476*PROJECTOR*Quad(Ye(2,2))*Sqr(g1)*Yu(2,2) -
      0.2759050553623628*PROJECTOR*Quad(g1)*Sqr(g2)*Yu(2,2) + 4.575063689774406
      *PROJECTOR*Quad(g3)*Sqr(g2)*Yu(2,2) - 2.0795715998030686*PROJECTOR*Quad(
      Yd(2,2))*Sqr(g2)*Yu(2,2) - 0.5509041095110554*PROJECTOR*Quad(Ye(2,2))*Sqr
      (g2)*Yu(2,2) + 0.09084432990862458*PROJECTOR*Lambdax*Sqr(g1)*Sqr(g2)*Yu(2
      ,2) - 1.3863544090646365*PROJECTOR*Quad(g1)*Sqr(g3)*Yu(2,2) -
      1.3089127804653757*PROJECTOR*Quad(g2)*Sqr(g3)*Yu(2,2) +
      0.3148250285229063*PROJECTOR*Quad(Yd(2,2))*Sqr(g3)*Yu(2,2) -
      0.9969582872023416*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(g3)*Yu(2,2) +
      0.1397605075517301*PROJECTOR*Sqr(g1)*Sqr(Lambdax)*Yu(2,2) +
      0.6988025377586506*PROJECTOR*Sqr(g2)*Sqr(Lambdax)*Yu(2,2) -
      0.20587853488920188*PROJECTOR*Quad(g1)*Sqr(Yd(2,2))*Yu(2,2) +
      0.40913682599408235*PROJECTOR*Quad(g2)*Sqr(Yd(2,2))*Yu(2,2) -
      12.757991577276114*PROJECTOR*Quad(g3)*Sqr(Yd(2,2))*Yu(2,2) +
      0.8230341000268552*PROJECTOR*Quad(Ye(2,2))*Sqr(Yd(2,2))*Yu(2,2) +
      0.5641039159322749*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(Yd(2,2))*Yu(2,2) -
      1.3643639389175877*PROJECTOR*Sqr(g1)*Sqr(g3)*Sqr(Yd(2,2))*Yu(2,2) -
      8.902567024109379*PROJECTOR*Sqr(g2)*Sqr(g3)*Sqr(Yd(2,2))*Yu(2,2) -
      1.129730769376485*PROJECTOR*Sqr(Lambdax)*Sqr(Yd(2,2))*Yu(2,2) -
      1.5400233128895373*PROJECTOR*Quad(g1)*Sqr(Ye(2,2))*Yu(2,2) -
      0.9718846178821886*PROJECTOR*Quad(g2)*Sqr(Ye(2,2))*Yu(2,2) +
      1.3665471849502502*PROJECTOR*Quad(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.3364704668151931*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2))*Yu(2,2) -
      0.3494012688793253*PROJECTOR*Sqr(Lambdax)*Sqr(Ye(2,2))*Yu(2,2) -
      0.14904312779958145*PROJECTOR*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.5159639826231229*PROJECTOR*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.445163098127733*PROJECTOR*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2))).
      real();


   return threeLoop * beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (2308.1827208150457*PROJECTOR*Power8(g3)*Yu(2,2)).real();


   return fourLoop * beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yu;
}

} // namespace flexiblesusy
