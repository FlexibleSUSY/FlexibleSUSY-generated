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


#include "HSSUSY_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.01*(-1200*traceYdAdjYdYdAdjYd - 400*traceYeAdjYeYeAdjYe
      - 1200*traceYuAdjYuYuAdjYu + 1200*traceYdAdjYd*Lambdax + 400*traceYeAdjYe
      *Lambdax + 1200*traceYuAdjYu*Lambdax + 27*Quad(g1) + 225*Quad(g2) - 180*
      Lambdax*Sqr(g1) - 900*Lambdax*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 1200*Sqr(
      Lambdax)));


   return oneLoop * beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe = TRACE_STRUCT.
      traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu = TRACE_STRUCT.
      traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.001*(60000*traceYdAdjYdYdAdjYdYdAdjYd - 24000*
      traceYdAdjYdYdAdjYuYuAdjYd + 12000*traceYdAdjYuYuAdjYdYdAdjYd - 12000*
      traceYdAdjYuYuAdjYuYuAdjYd + 20000*traceYeAdjYeYeAdjYeYeAdjYe + 60000*
      traceYuAdjYuYuAdjYuYuAdjYu - 78000*Cube(Lambdax) - 3000*
      traceYdAdjYdYdAdjYd*Lambdax - 42000*traceYdAdjYuYuAdjYd*Lambdax - 1000*
      traceYeAdjYeYeAdjYe*Lambdax - 3000*traceYuAdjYuYuAdjYu*Lambdax - 3411*
      Power6(g1) + 38125*Power6(g2) + 900*traceYdAdjYd*Quad(g1) - 4500*
      traceYeAdjYe*Quad(g1) - 3420*traceYuAdjYu*Quad(g1) + 9435*Lambdax*Quad(g1
      ) - 4500*traceYdAdjYd*Quad(g2) - 1500*traceYeAdjYe*Quad(g2) - 4500*
      traceYuAdjYu*Quad(g2) - 9125*Lambdax*Quad(g2) + 1600*traceYdAdjYdYdAdjYd*
      Sqr(g1) - 4800*traceYeAdjYeYeAdjYe*Sqr(g1) - 3200*traceYuAdjYuYuAdjYu*Sqr
      (g1) + 2500*traceYdAdjYd*Lambdax*Sqr(g1) + 7500*traceYeAdjYe*Lambdax*Sqr(
      g1) + 8500*traceYuAdjYu*Lambdax*Sqr(g1) - 7225*Quad(g2)*Sqr(g1) + 22500*
      traceYdAdjYd*Lambdax*Sqr(g2) + 7500*traceYeAdjYe*Lambdax*Sqr(g2) + 22500*
      traceYuAdjYu*Lambdax*Sqr(g2) - 8385*Quad(g1)*Sqr(g2) + 5400*traceYdAdjYd*
      Sqr(g1)*Sqr(g2) + 6600*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 12600*traceYuAdjYu*
      Sqr(g1)*Sqr(g2) + 5850*Lambdax*Sqr(g1)*Sqr(g2) - 64000*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 64000*traceYuAdjYuYuAdjYu*Sqr(g3) + 80000*
      traceYdAdjYd*Lambdax*Sqr(g3) + 80000*traceYuAdjYu*Lambdax*Sqr(g3) - 72000
      *traceYdAdjYd*Sqr(Lambdax) - 24000*traceYeAdjYe*Sqr(Lambdax) - 72000*
      traceYuAdjYu*Sqr(Lambdax) + 10800*Sqr(g1)*Sqr(Lambdax) + 54000*Sqr(g2)*
      Sqr(Lambdax)));


   return twoLoop * beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   const double beta_Lambdax_1 = Re(-6.031600412410262*(-9.339962685289576*
      Lambdax*Power6(g1) - 286.9827742168705*Lambdax*Power6(g2) +
      74.07064511965402*Lambdax*Power6(Yd(2,2)) + 77.74412745450469*Lambdax*
      Power6(Ye(2,2)) + 1.*Power8(g1) + 75.66229767543491*Power8(g2) +
      161.25010404433053*Power8(Yd(2,2)) + 21.420306803784438*Power8(Ye(2,2)) -
      4.310550825785817*Quad(g1)*Quad(g2) - 249.16051734798086*Quad(Lambdax) +
      11.266458319444492*Quad(g1)*Quad(Yd(2,2)) - 10.53370494294467*Quad(g2)*
      Quad(Yd(2,2)) + 33.291890929823126*Quad(g3)*Quad(Yd(2,2)) -
      21.23162373590703*Quad(g1)*Quad(Ye(2,2)) - 3.262544766538025*Quad(g2)*
      Quad(Ye(2,2)) + 47.74852117315802*Quad(Yd(2,2))*Quad(Ye(2,2)) +
      6.423699789119175*Cube(Lambdax)*Sqr(g1) + 25.12726423496584*Power6(g2)*
      Sqr(g1) - 9.226669872878011*Power6(Yd(2,2))*Sqr(g1) - 17.462679807754416*
      Power6(Ye(2,2))*Sqr(g1) - 26.40696244250939*Lambdax*Quad(g2)*Sqr(g1) -
      11.990947297383288*Lambdax*Quad(Yd(2,2))*Sqr(g1) + 15.37217037999308*
      Lambdax*Quad(Ye(2,2))*Sqr(g1) + 32.118498945595874*Cube(Lambdax)*Sqr(g2)
      + 1.0231901443099605*Power6(g1)*Sqr(g2) - 49.16652864612754*Power6(Yd(2,2
      ))*Sqr(g2) - 16.388842882042518*Power6(Ye(2,2))*Sqr(g2) -
      20.47661959122464*Lambdax*Quad(g1)*Sqr(g2) + 1.8137171911699768*Lambdax*
      Quad(Yd(2,2))*Sqr(g2) - 2.3797101762657173*Lambdax*Quad(Ye(2,2))*Sqr(g2)
      + 25.306618248571045*Quad(Yd(2,2))*Sqr(g1)*Sqr(g2) + 22.872470753746928*
      Quad(Ye(2,2))*Sqr(g1)*Sqr(g2) - 0.4397999591948324*Power6(g1)*Sqr(g3) -
      4.997726809032187*Power6(g2)*Sqr(g3) - 166.12085657591095*Power6(Yd(2,2))
      *Sqr(g3) + 2.7790841962446198*Lambdax*Quad(g1)*Sqr(g3) +
      18.948301338031502*Lambdax*Quad(g2)*Sqr(g3) + 219.79763285742905*Lambdax*
      Quad(Yd(2,2))*Sqr(g3) - 0.9995453618064374*Quad(g2)*Sqr(g1)*Sqr(g3) -
      7.513283589518858*Quad(Yd(2,2))*Sqr(g1)*Sqr(g3) - 0.7329999319913874*Quad
      (g1)*Sqr(g2)*Sqr(g3) - 8.852950966952918*Quad(Yd(2,2))*Sqr(g2)*Sqr(g3) +
      30.759991811525552*Quad(g1)*Sqr(Lambdax) + 131.0233001667084*Quad(g2)*Sqr
      (Lambdax) - 293.1651465422671*Quad(Yd(2,2))*Sqr(Lambdax) -
      109.65884580737854*Quad(Ye(2,2))*Sqr(Lambdax) + 52.49683489237998*Sqr(g1)
      *Sqr(g2)*Sqr(Lambdax) - 72.36885240306762*Cube(Lambdax)*Sqr(Yd(2,2)) -
      3.23203546695533*Power6(g1)*Sqr(Yd(2,2)) - 41.447585739271474*Power6(g2)*
      Sqr(Yd(2,2)) + 24.620331229909603*Power6(Ye(2,2))*Sqr(Yd(2,2)) +
      11.460062483294552*Lambdax*Quad(g1)*Sqr(Yd(2,2)) + 105.99644858660962*
      Lambdax*Quad(g2)*Sqr(Yd(2,2)) - 118.365909843938*Lambdax*Quad(g3)*Sqr(Yd(
      2,2)) - 79.58086862193002*Lambdax*Quad(Ye(2,2))*Sqr(Yd(2,2)) -
      13.784493073708608*Quad(g2)*Sqr(g1)*Sqr(Yd(2,2)) - 8.381888765119909*Quad
      (g1)*Sqr(g2)*Sqr(Yd(2,2)) + 15.837899406921814*Lambdax*Sqr(g1)*Sqr(g2)*
      Sqr(Yd(2,2)) - 1.6285371043301784*Quad(g1)*Sqr(g3)*Sqr(Yd(2,2)) -
      10.91844691535382*Quad(g2)*Sqr(g3)*Sqr(Yd(2,2)) + 1.3873490852923183*
      Lambdax*Sqr(g1)*Sqr(g3)*Sqr(Yd(2,2)) - 5.021649329193744*Lambdax*Sqr(g2)*
      Sqr(g3)*Sqr(Yd(2,2)) - 5.958996138580129*Sqr(g1)*Sqr(g2)*Sqr(g3)*Sqr(Yd(2
      ,2)) + 12.5881938544474*Sqr(g1)*Sqr(Lambdax)*Sqr(Yd(2,2)) +
      59.60915139954888*Sqr(g2)*Sqr(Lambdax)*Sqr(Yd(2,2)) - 26.654542981504992*
      Sqr(g3)*Sqr(Lambdax)*Sqr(Yd(2,2)) - 24.12295080102254*Cube(Lambdax)*Sqr(
      Ye(2,2)) - 5.492026720736577*Power6(g1)*Sqr(Ye(2,2)) - 13.81586191309049*
      Power6(g2)*Sqr(Ye(2,2)) + 24.620331229909603*Power6(Yd(2,2))*Sqr(Ye(2,2))
      + 17.969663482497715*Lambdax*Quad(g1)*Sqr(Ye(2,2)) + 35.33214952886986*
      Lambdax*Quad(g2)*Sqr(Ye(2,2)) - 79.58086862193002*Lambdax*Quad(Yd(2,2))*
      Sqr(Ye(2,2)) - 1.5347452609060228*Quad(g2)*Sqr(g1)*Sqr(Ye(2,2)) -
      3.459017345205606*Quad(g1)*Sqr(g2)*Sqr(Ye(2,2)) - 6.687860766560821*
      Lambdax*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2)) + 1.9748527030237164*Sqr(g1)*Sqr(
      Lambdax)*Sqr(Ye(2,2)) + 19.86971713318296*Sqr(g2)*Sqr(Lambdax)*Sqr(Ye(2,2
      )) - 0.4078519516873914*Quad(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      0.746070643330594*Quad(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 1.7905695439934257
      *Lambdax*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 8.952847719967128*Lambdax*
      Sqr(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 0.4973804288870627*Sqr(g1)*Sqr(g2)*
      Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 35.81139087986851*Sqr(Lambdax)*Sqr(Yd(2,2))*
      Sqr(Ye(2,2)) - 7.372368699637704*Power6(g1)*Sqr(Yu(2,2)) -
      41.447585739271474*Power6(g2)*Sqr(Yu(2,2)) + 24.822556597226036*Lambdax*
      Quad(g1)*Sqr(Yu(2,2)) + 105.99644858660962*Lambdax*Quad(g2)*Sqr(Yu(2,2))
      - 118.365909843938*Lambdax*Quad(g3)*Sqr(Yu(2,2)) - 8.648442082445648*Quad
      (g2)*Sqr(g1)*Sqr(Yu(2,2)) - 7.047713089694501*Quad(g1)*Sqr(g2)*Sqr(Yu(2,2
      )) - 1.8617589799206342*Lambdax*Sqr(g1)*Sqr(g2)*Sqr(Yu(2,2)) -
      0.6735666808670179*Quad(g1)*Sqr(g3)*Sqr(Yu(2,2)) - 10.91844691535382*Quad
      (g2)*Sqr(g3)*Sqr(Yu(2,2)) - 5.787687768100881*Lambdax*Sqr(g1)*Sqr(g3)*Sqr
      (Yu(2,2)) - 5.021649329193744*Lambdax*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) -
      7.55061351101873*Sqr(g1)*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) +
      10.589103130834523*Sqr(g1)*Sqr(Lambdax)*Sqr(Yu(2,2)) + 59.60915139954888*
      Sqr(g2)*Sqr(Lambdax)*Sqr(Yu(2,2)) - 26.654542981504992*Sqr(g3)*Sqr(
      Lambdax)*Sqr(Yu(2,2))));
   const double beta_Lambdax_2 = Re(436.5*(-1.023515541010766*Lambdax*Power6(Yu
      (2,2)) - 2.228169974925502*Power8(Yu(2,2)) + 0.14614411445535522*Quad(g1)
      *Quad(Yu(2,2)) + 0.14555578253854087*Quad(g2)*Quad(Yu(2,2)) -
      0.4600306599364001*Quad(g3)*Quad(Yu(2,2)) + 0.7931097093011755*Quad(Yd(2,
      2))*Quad(Yu(2,2)) - 0.6597938144329897*Quad(Ye(2,2))*Quad(Yu(2,2)) +
      0.31092536460231723*Power6(Yu(2,2))*Sqr(g1) - 0.09628707405242558*Lambdax
      *Quad(Yu(2,2))*Sqr(g1) + 0.6793879827233994*Power6(Yu(2,2))*Sqr(g2) -
      0.025062124532088004*Lambdax*Quad(Yu(2,2))*Sqr(g2) - 0.6447330669284618*
      Quad(Yu(2,2))*Sqr(g1)*Sqr(g2) + 2.2954745178309515*Power6(Yu(2,2))*Sqr(g3
      ) - 3.0371855509499848*Lambdax*Quad(Yu(2,2))*Sqr(g3) +
      0.16101073097774002*Quad(Yu(2,2))*Sqr(g1)*Sqr(g3) + 0.12233095693773434*
      Quad(Yu(2,2))*Sqr(g2)*Sqr(g3) + 4.050985151864039*Quad(Yu(2,2))*Sqr(
      Lambdax) - 1.2178606965749863*Power6(Yu(2,2))*Sqr(Yd(2,2)) +
      4.458058162909423*Lambdax*Quad(Yu(2,2))*Sqr(Yd(2,2)) -
      0.10848390353339003*Quad(Yu(2,2))*Sqr(g1)*Sqr(Yd(2,2)) +
      0.13659793814432988*Quad(Yu(2,2))*Sqr(g2)*Sqr(Yd(2,2)) -
      0.5470674121572557*Quad(Yu(2,2))*Sqr(g3)*Sqr(Yd(2,2)) -
      0.3402061855670103*Power6(Yu(2,2))*Sqr(Ye(2,2)) + 1.099656357388316*
      Lambdax*Quad(Yu(2,2))*Sqr(Ye(2,2)) + 0.05154639175257732*Quad(Yu(2,2))*
      Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 1.*Cube(Lambdax)*Sqr(Yu(2,2)) -
      1.2178606965749863*Power6(Yd(2,2))*Sqr(Yu(2,2)) - 0.3402061855670103*
      Power6(Ye(2,2))*Sqr(Yu(2,2)) + 4.458058162909423*Lambdax*Quad(Yd(2,2))*
      Sqr(Yu(2,2)) + 1.099656357388316*Lambdax*Quad(Ye(2,2))*Sqr(Yu(2,2)) +
      0.040168042971726645*Quad(Yd(2,2))*Sqr(g1)*Sqr(Yu(2,2)) +
      0.13659793814432988*Quad(Yd(2,2))*Sqr(g2)*Sqr(Yu(2,2)) -
      0.5470674121572557*Quad(Yd(2,2))*Sqr(g3)*Sqr(Yu(2,2)) -
      0.0405119402990832*Quad(g1)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.2642469996319474*Quad(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      1.7594501718213058*Quad(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.10996563573883161*Quad(Ye(2,2))*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.2194385717470402*Lambdax*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.31083128169247154*Lambdax*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.15977437983399062*Sqr(g1)*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.15302388409310907*Lambdax*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.9841691885756799*Sqr(g2)*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      2.1112878907901247*Sqr(Lambdax)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.0963573883161512*Quad(g1)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.010309278350515464*Quad(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.05154639175257732*Quad(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.024742268041237116*Lambdax*Sqr(g1)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.12371134020618557*Lambdax*Sqr(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.03986254295532646*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.4948453608247423*Sqr(Lambdax)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.09621993127147767*Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2))));

   beta_Lambdax = beta_Lambdax_1 + beta_Lambdax_2;


   return threeLoop * beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = Re(16616.342695047893*Power6(g3)*Quad(Yu(2,2)));


   return fourLoop * beta_Lambdax;
}

/**
 * Calculates the 5-loop beta function of Lambdax.
 *
 * @return 5-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return fiveLoop * beta_Lambdax;
}

} // namespace flexiblesusy
