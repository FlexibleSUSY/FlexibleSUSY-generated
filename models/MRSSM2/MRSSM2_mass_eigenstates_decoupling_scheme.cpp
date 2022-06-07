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
 * @file MRSSM2_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the MRSSM2 model class in the decoupling scheme
 *
 * Contains the definition of the MRSSM2 model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "MRSSM2_mass_eigenstates_decoupling_scheme.hpp"
#include "MRSSM2_mass_eigenstates.hpp"
#include "MRSSM2_info.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "error.hpp"
#include "ewsb_solver.hpp"
#include "ew_input.hpp"
#include "functors.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "numerics2.hpp"
#include "raii.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace flexiblesusy {

#define CLASSNAME MRSSM2_mass_eigenstates_decoupling_scheme

#define DERIVEDPARAMETER(p) model.p()
#define EXTRAPARAMETER(parameter) model.get_##parameter()
#define INPUT(parameter) model.get_input().parameter
#define INPUTPARAMETER(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define LowEnergyGaugeCoupling(i) new_g##i
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODELPARAMETER(parameter) model.get_##parameter()
#define PHASE(p) model.get_##p()
#define PHYSICAL(parameter) physical.parameter

CLASSNAME::CLASSNAME(const MRSSM2_input_parameters& input_)
   : MRSSM2_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const MRSSM2_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<MRSSM2_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<MRSSM2_mass_eigenstates_decoupling_scheme>(*this);
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
}

void CLASSNAME::fill_from(const standard_model::Standard_model& sm_input)
{
   using RM33 = Eigen::Matrix<double, 3, 3>;
   using CM33 = Eigen::Matrix<std::complex<double>, 3, 3>;

   // make a local copy and ensure that the tree-level masses are calculated
   auto sm = sm_input;
   sm.calculate_DRbar_masses();

   const auto sm_gY = sm.get_g1()*standard_model_info::normalization_g1;
   const auto sm_g2 = sm.get_g2()*standard_model_info::normalization_g2;
   const auto sm_g3 = sm.get_g3()*standard_model_info::normalization_g3;
   const auto VEV   = sm.get_v();

   // SM interface parameters for the low-scale constraint
   const CM33 CKM = sm.get_Vu().adjoint() * sm.get_Vd();
   const double MZMSbar = sm.get_MVZ();
   const double MZDRbar = sm.get_MVZ();
   const double MWMSbar = sm.get_MVWp();
   const double MWDRbar = sm.get_MVWp();
   const double EDRbar = sm_gY * sm_g2 / std::sqrt(sm_gY*sm_gY + sm_g2*sm_g2);
   const double EMSbar = EDRbar;
   const double THETAW = sm.ThetaW();
   const double ThetaWDRbar = THETAW;
   const double AlphaS = sm_g3*sm_g3*0.07957747154594767; // g3^2/(4 Pi)

   RM33 upQuarksDRbar(RM33::Zero());
   RM33 downQuarksDRbar(RM33::Zero());
   RM33 downLeptonsDRbar(RM33::Zero());

   upQuarksDRbar.diagonal()    = sm.get_MFu();
   downQuarksDRbar.diagonal()  = sm.get_MFd();
   downLeptonsDRbar.diagonal() = sm.get_MFe();

   // new gauge couplings
   double new_g1 = 0., new_g2 = 0., new_g3 = 0.;

   // calculate new gauge couplings
   {
      auto model = this;

      new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
      new_g2 = EDRbar*Csc(ThetaWDRbar);
      new_g3 = 3.5449077018110318*Sqrt(AlphaS);

      if (IsFinite(new_g1)) {
         model->get_problems().unflag_non_perturbative_parameter(MRSSM2_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            MRSSM2_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(MRSSM2_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            MRSSM2_info::g2, new_g2, get_scale());
         new_g2 = Electroweak_constants::g2;
      }
   }

   // set new gauge couplings
   {
      auto& model = *this;
      auto MODEL = this;
      
      MODEL->set_g1(new_g1);
      MODEL->set_g2(new_g2);
      MODEL->set_g3(new_g3);

   }

   // apply user-defined low-energy constraint for the VEV(s)
   {
      auto& model = *this;
      auto MODEL = this;
      const auto TanBeta = INPUTPARAMETER(TanBeta);
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);

      MODEL->set_vd(Re((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
         ))));
      MODEL->set_vu(Re((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
         TanBeta)))));

   }

   // apply user-defined low-energy constraint for the Yukawa couplings
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vu = MODELPARAMETER(vu);
      MODEL->set_Yu(((1.4142135623730951*(upQuarksDRbar*CKM))/vu).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Yd(((1.4142135623730951*downQuarksDRbar)/vd).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Ye(((1.4142135623730951*downLeptonsDRbar)/vd).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const MRSSM2_mass_eigenstates& model)
{
   set(model.get());
   set_scale(model.get_scale());
   set_loops(model.get_loops());
   set_thresholds(model.get_thresholds());
   set_zero_threshold(model.get_zero_threshold());
   set_input_parameters(model.get_input());
   force_output = model.do_force_output();
   precision = model.get_precision();
   physical = model.get_physical();

#define OTHER(p) model.get_##p()
   MVG = OTHER(MVG);
   MGlu = OTHER(MGlu);
   MFv = OTHER(MFv);
   MSRdp = OTHER(MSRdp);
   MSRum = OTHER(MSRum);
   MsigmaO = OTHER(MsigmaO);
   MphiO = OTHER(MphiO);
   MSd = OTHER(MSd);
   ZD = OTHER(ZD);
   MSv = OTHER(MSv);
   ZV = OTHER(ZV);
   MSu = OTHER(MSu);
   ZU = OTHER(ZU);
   MSe = OTHER(MSe);
   ZE = OTHER(ZE);
   Mhh = OTHER(Mhh);
   ZH = OTHER(ZH);
   MAh = OTHER(MAh);
   ZA = OTHER(ZA);
   MRh = OTHER(MRh);
   ZHR = OTHER(ZHR);
   MHpm = OTHER(MHpm);
   ZP = OTHER(ZP);
   MChi = OTHER(MChi);
   ZN1 = OTHER(ZN1);
   ZN2 = OTHER(ZN2);
   MCha1 = OTHER(MCha1);
   UM1 = OTHER(UM1);
   UP1 = OTHER(UP1);
   MCha2 = OTHER(MCha2);
   UM2 = OTHER(UM2);
   UP2 = OTHER(UP2);
   MFe = OTHER(MFe);
   ZEL = OTHER(ZEL);
   ZER = OTHER(ZER);
   MFd = OTHER(MFd);
   ZDL = OTHER(ZDL);
   ZDR = OTHER(ZDR);
   MFu = OTHER(MFu);
   ZUL = OTHER(ZUL);
   ZUR = OTHER(ZUR);
   MVWm = OTHER(MVWm);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const MRSSM2_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MRSSM2_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MRSSM2_physical& physical_)
{
   physical = physical_;
}

const Problems& CLASSNAME::get_problems() const
{
   return problems;
}

Problems& CLASSNAME::get_problems()
{
   return problems;
}

int CLASSNAME::solve_ewsb_equations_tree_level()
{
   int error = EWSB_solver::SUCCESS;

   
   const double old_vT = vT;
   const double old_vS = vS;
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

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

   if (!is_finite) {
      vT = old_vT;
      vS = old_vS;
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      error = EWSB_solver::FAIL;
   }

   return error;
}

int CLASSNAME::solve_ewsb_equations()
{
   return solve_ewsb_equations_tree_level();
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "MRSSM2\n"
           "========================================\n";
   MRSSM2_soft_parameters::print(ostr);
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
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

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
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_tree_level_mass_spectrum()
{
   const auto save_mHd2_raii = make_raii_save(mHd2);
   const auto save_mHu2_raii = make_raii_save(mHu2);
   const auto save_vS_raii = make_raii_save(vS);
   const auto save_vT_raii = make_raii_save(vT);

   solve_ewsb_equations_tree_level();

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

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 *
 * @note Does currently nothing, because it is not clear how to
 * calculate the pole masses in this scheme.
 */
void CLASSNAME::calculate_pole_mass_spectrum()
{
   calculate_tree_level_mass_spectrum();
   // move goldstone bosons to the front
   reorder_tree_level_masses();
   copy_tree_level_masses_to_pole_masses();
   check_pole_masses_for_tachyons();
}

void CLASSNAME::copy_tree_level_masses_to_pole_masses()
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
   PHYSICAL(ZZ) = ZZ;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_tree_level_masses()
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
   if (PHYSICAL(MSRdp) < 0.) { problems.flag_pole_tachyon(MRSSM2_info::SRdp); }
   if (PHYSICAL(MSRum) < 0.) { problems.flag_pole_tachyon(MRSSM2_info::SRum); }
   if (PHYSICAL(MsigmaO) < 0.) { problems.flag_pole_tachyon(MRSSM2_info::sigmaO); }
   if (PHYSICAL(MphiO) < 0.) { problems.flag_pole_tachyon(MRSSM2_info::phiO); }
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Se); }
   if (PHYSICAL(Mhh).tail<4>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::hh); }
   if (PHYSICAL(MAh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Ah); }
   if (PHYSICAL(MRh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Rh); }
   if (PHYSICAL(MHpm).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(MRSSM2_info::Hpm); }
}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_mass_spectrum()
{
   calculate_tree_level_mass_spectrum();
   calculate_pole_mass_spectrum();
}

void CLASSNAME::clear_tree_level_parameters()
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
   problems.clear();
}

void CLASSNAME::clear()
{
   MRSSM2_soft_parameters::clear();
   clear_tree_level_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& pars)
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

const MRSSM2_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

MRSSM2_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
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

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZD(0,0) = pars(64);
   ZD(0,1) = pars(65);
   ZD(0,2) = pars(66);
   ZD(0,3) = pars(67);
   ZD(0,4) = pars(68);
   ZD(0,5) = pars(69);
   ZD(1,0) = pars(70);
   ZD(1,1) = pars(71);
   ZD(1,2) = pars(72);
   ZD(1,3) = pars(73);
   ZD(1,4) = pars(74);
   ZD(1,5) = pars(75);
   ZD(2,0) = pars(76);
   ZD(2,1) = pars(77);
   ZD(2,2) = pars(78);
   ZD(2,3) = pars(79);
   ZD(2,4) = pars(80);
   ZD(2,5) = pars(81);
   ZD(3,0) = pars(82);
   ZD(3,1) = pars(83);
   ZD(3,2) = pars(84);
   ZD(3,3) = pars(85);
   ZD(3,4) = pars(86);
   ZD(3,5) = pars(87);
   ZD(4,0) = pars(88);
   ZD(4,1) = pars(89);
   ZD(4,2) = pars(90);
   ZD(4,3) = pars(91);
   ZD(4,4) = pars(92);
   ZD(4,5) = pars(93);
   ZD(5,0) = pars(94);
   ZD(5,1) = pars(95);
   ZD(5,2) = pars(96);
   ZD(5,3) = pars(97);
   ZD(5,4) = pars(98);
   ZD(5,5) = pars(99);
   ZV(0,0) = pars(100);
   ZV(0,1) = pars(101);
   ZV(0,2) = pars(102);
   ZV(1,0) = pars(103);
   ZV(1,1) = pars(104);
   ZV(1,2) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZU(0,0) = pars(109);
   ZU(0,1) = pars(110);
   ZU(0,2) = pars(111);
   ZU(0,3) = pars(112);
   ZU(0,4) = pars(113);
   ZU(0,5) = pars(114);
   ZU(1,0) = pars(115);
   ZU(1,1) = pars(116);
   ZU(1,2) = pars(117);
   ZU(1,3) = pars(118);
   ZU(1,4) = pars(119);
   ZU(1,5) = pars(120);
   ZU(2,0) = pars(121);
   ZU(2,1) = pars(122);
   ZU(2,2) = pars(123);
   ZU(2,3) = pars(124);
   ZU(2,4) = pars(125);
   ZU(2,5) = pars(126);
   ZU(3,0) = pars(127);
   ZU(3,1) = pars(128);
   ZU(3,2) = pars(129);
   ZU(3,3) = pars(130);
   ZU(3,4) = pars(131);
   ZU(3,5) = pars(132);
   ZU(4,0) = pars(133);
   ZU(4,1) = pars(134);
   ZU(4,2) = pars(135);
   ZU(4,3) = pars(136);
   ZU(4,4) = pars(137);
   ZU(4,5) = pars(138);
   ZU(5,0) = pars(139);
   ZU(5,1) = pars(140);
   ZU(5,2) = pars(141);
   ZU(5,3) = pars(142);
   ZU(5,4) = pars(143);
   ZU(5,5) = pars(144);
   ZE(0,0) = pars(145);
   ZE(0,1) = pars(146);
   ZE(0,2) = pars(147);
   ZE(0,3) = pars(148);
   ZE(0,4) = pars(149);
   ZE(0,5) = pars(150);
   ZE(1,0) = pars(151);
   ZE(1,1) = pars(152);
   ZE(1,2) = pars(153);
   ZE(1,3) = pars(154);
   ZE(1,4) = pars(155);
   ZE(1,5) = pars(156);
   ZE(2,0) = pars(157);
   ZE(2,1) = pars(158);
   ZE(2,2) = pars(159);
   ZE(2,3) = pars(160);
   ZE(2,4) = pars(161);
   ZE(2,5) = pars(162);
   ZE(3,0) = pars(163);
   ZE(3,1) = pars(164);
   ZE(3,2) = pars(165);
   ZE(3,3) = pars(166);
   ZE(3,4) = pars(167);
   ZE(3,5) = pars(168);
   ZE(4,0) = pars(169);
   ZE(4,1) = pars(170);
   ZE(4,2) = pars(171);
   ZE(4,3) = pars(172);
   ZE(4,4) = pars(173);
   ZE(4,5) = pars(174);
   ZE(5,0) = pars(175);
   ZE(5,1) = pars(176);
   ZE(5,2) = pars(177);
   ZE(5,3) = pars(178);
   ZE(5,4) = pars(179);
   ZE(5,5) = pars(180);
   ZH(0,0) = pars(181);
   ZH(0,1) = pars(182);
   ZH(0,2) = pars(183);
   ZH(0,3) = pars(184);
   ZH(1,0) = pars(185);
   ZH(1,1) = pars(186);
   ZH(1,2) = pars(187);
   ZH(1,3) = pars(188);
   ZH(2,0) = pars(189);
   ZH(2,1) = pars(190);
   ZH(2,2) = pars(191);
   ZH(2,3) = pars(192);
   ZH(3,0) = pars(193);
   ZH(3,1) = pars(194);
   ZH(3,2) = pars(195);
   ZH(3,3) = pars(196);
   ZA(0,0) = pars(197);
   ZA(0,1) = pars(198);
   ZA(0,2) = pars(199);
   ZA(0,3) = pars(200);
   ZA(1,0) = pars(201);
   ZA(1,1) = pars(202);
   ZA(1,2) = pars(203);
   ZA(1,3) = pars(204);
   ZA(2,0) = pars(205);
   ZA(2,1) = pars(206);
   ZA(2,2) = pars(207);
   ZA(2,3) = pars(208);
   ZA(3,0) = pars(209);
   ZA(3,1) = pars(210);
   ZA(3,2) = pars(211);
   ZA(3,3) = pars(212);
   ZHR(0,0) = pars(213);
   ZHR(0,1) = pars(214);
   ZHR(1,0) = pars(215);
   ZHR(1,1) = pars(216);
   ZP(0,0) = pars(217);
   ZP(0,1) = pars(218);
   ZP(0,2) = pars(219);
   ZP(0,3) = pars(220);
   ZP(1,0) = pars(221);
   ZP(1,1) = pars(222);
   ZP(1,2) = pars(223);
   ZP(1,3) = pars(224);
   ZP(2,0) = pars(225);
   ZP(2,1) = pars(226);
   ZP(2,2) = pars(227);
   ZP(2,3) = pars(228);
   ZP(3,0) = pars(229);
   ZP(3,1) = pars(230);
   ZP(3,2) = pars(231);
   ZP(3,3) = pars(232);
   ZN1(0,0) = std::complex<double>(pars(233), pars(234));
   ZN1(0,1) = std::complex<double>(pars(235), pars(236));
   ZN1(0,2) = std::complex<double>(pars(237), pars(238));
   ZN1(0,3) = std::complex<double>(pars(239), pars(240));
   ZN1(1,0) = std::complex<double>(pars(241), pars(242));
   ZN1(1,1) = std::complex<double>(pars(243), pars(244));
   ZN1(1,2) = std::complex<double>(pars(245), pars(246));
   ZN1(1,3) = std::complex<double>(pars(247), pars(248));
   ZN1(2,0) = std::complex<double>(pars(249), pars(250));
   ZN1(2,1) = std::complex<double>(pars(251), pars(252));
   ZN1(2,2) = std::complex<double>(pars(253), pars(254));
   ZN1(2,3) = std::complex<double>(pars(255), pars(256));
   ZN1(3,0) = std::complex<double>(pars(257), pars(258));
   ZN1(3,1) = std::complex<double>(pars(259), pars(260));
   ZN1(3,2) = std::complex<double>(pars(261), pars(262));
   ZN1(3,3) = std::complex<double>(pars(263), pars(264));
   ZN2(0,0) = std::complex<double>(pars(265), pars(266));
   ZN2(0,1) = std::complex<double>(pars(267), pars(268));
   ZN2(0,2) = std::complex<double>(pars(269), pars(270));
   ZN2(0,3) = std::complex<double>(pars(271), pars(272));
   ZN2(1,0) = std::complex<double>(pars(273), pars(274));
   ZN2(1,1) = std::complex<double>(pars(275), pars(276));
   ZN2(1,2) = std::complex<double>(pars(277), pars(278));
   ZN2(1,3) = std::complex<double>(pars(279), pars(280));
   ZN2(2,0) = std::complex<double>(pars(281), pars(282));
   ZN2(2,1) = std::complex<double>(pars(283), pars(284));
   ZN2(2,2) = std::complex<double>(pars(285), pars(286));
   ZN2(2,3) = std::complex<double>(pars(287), pars(288));
   ZN2(3,0) = std::complex<double>(pars(289), pars(290));
   ZN2(3,1) = std::complex<double>(pars(291), pars(292));
   ZN2(3,2) = std::complex<double>(pars(293), pars(294));
   ZN2(3,3) = std::complex<double>(pars(295), pars(296));
   UM1(0,0) = std::complex<double>(pars(297), pars(298));
   UM1(0,1) = std::complex<double>(pars(299), pars(300));
   UM1(1,0) = std::complex<double>(pars(301), pars(302));
   UM1(1,1) = std::complex<double>(pars(303), pars(304));
   UP1(0,0) = std::complex<double>(pars(305), pars(306));
   UP1(0,1) = std::complex<double>(pars(307), pars(308));
   UP1(1,0) = std::complex<double>(pars(309), pars(310));
   UP1(1,1) = std::complex<double>(pars(311), pars(312));
   UM2(0,0) = std::complex<double>(pars(313), pars(314));
   UM2(0,1) = std::complex<double>(pars(315), pars(316));
   UM2(1,0) = std::complex<double>(pars(317), pars(318));
   UM2(1,1) = std::complex<double>(pars(319), pars(320));
   UP2(0,0) = std::complex<double>(pars(321), pars(322));
   UP2(0,1) = std::complex<double>(pars(323), pars(324));
   UP2(1,0) = std::complex<double>(pars(325), pars(326));
   UP2(1,1) = std::complex<double>(pars(327), pars(328));
   ZEL(0,0) = std::complex<double>(pars(329), pars(330));
   ZEL(0,1) = std::complex<double>(pars(331), pars(332));
   ZEL(0,2) = std::complex<double>(pars(333), pars(334));
   ZEL(1,0) = std::complex<double>(pars(335), pars(336));
   ZEL(1,1) = std::complex<double>(pars(337), pars(338));
   ZEL(1,2) = std::complex<double>(pars(339), pars(340));
   ZEL(2,0) = std::complex<double>(pars(341), pars(342));
   ZEL(2,1) = std::complex<double>(pars(343), pars(344));
   ZEL(2,2) = std::complex<double>(pars(345), pars(346));
   ZER(0,0) = std::complex<double>(pars(347), pars(348));
   ZER(0,1) = std::complex<double>(pars(349), pars(350));
   ZER(0,2) = std::complex<double>(pars(351), pars(352));
   ZER(1,0) = std::complex<double>(pars(353), pars(354));
   ZER(1,1) = std::complex<double>(pars(355), pars(356));
   ZER(1,2) = std::complex<double>(pars(357), pars(358));
   ZER(2,0) = std::complex<double>(pars(359), pars(360));
   ZER(2,1) = std::complex<double>(pars(361), pars(362));
   ZER(2,2) = std::complex<double>(pars(363), pars(364));
   ZDL(0,0) = std::complex<double>(pars(365), pars(366));
   ZDL(0,1) = std::complex<double>(pars(367), pars(368));
   ZDL(0,2) = std::complex<double>(pars(369), pars(370));
   ZDL(1,0) = std::complex<double>(pars(371), pars(372));
   ZDL(1,1) = std::complex<double>(pars(373), pars(374));
   ZDL(1,2) = std::complex<double>(pars(375), pars(376));
   ZDL(2,0) = std::complex<double>(pars(377), pars(378));
   ZDL(2,1) = std::complex<double>(pars(379), pars(380));
   ZDL(2,2) = std::complex<double>(pars(381), pars(382));
   ZDR(0,0) = std::complex<double>(pars(383), pars(384));
   ZDR(0,1) = std::complex<double>(pars(385), pars(386));
   ZDR(0,2) = std::complex<double>(pars(387), pars(388));
   ZDR(1,0) = std::complex<double>(pars(389), pars(390));
   ZDR(1,1) = std::complex<double>(pars(391), pars(392));
   ZDR(1,2) = std::complex<double>(pars(393), pars(394));
   ZDR(2,0) = std::complex<double>(pars(395), pars(396));
   ZDR(2,1) = std::complex<double>(pars(397), pars(398));
   ZDR(2,2) = std::complex<double>(pars(399), pars(400));
   ZUL(0,0) = std::complex<double>(pars(401), pars(402));
   ZUL(0,1) = std::complex<double>(pars(403), pars(404));
   ZUL(0,2) = std::complex<double>(pars(405), pars(406));
   ZUL(1,0) = std::complex<double>(pars(407), pars(408));
   ZUL(1,1) = std::complex<double>(pars(409), pars(410));
   ZUL(1,2) = std::complex<double>(pars(411), pars(412));
   ZUL(2,0) = std::complex<double>(pars(413), pars(414));
   ZUL(2,1) = std::complex<double>(pars(415), pars(416));
   ZUL(2,2) = std::complex<double>(pars(417), pars(418));
   ZUR(0,0) = std::complex<double>(pars(419), pars(420));
   ZUR(0,1) = std::complex<double>(pars(421), pars(422));
   ZUR(0,2) = std::complex<double>(pars(423), pars(424));
   ZUR(1,0) = std::complex<double>(pars(425), pars(426));
   ZUR(1,1) = std::complex<double>(pars(427), pars(428));
   ZUR(1,2) = std::complex<double>(pars(429), pars(430));
   ZUR(2,0) = std::complex<double>(pars(431), pars(432));
   ZUR(2,1) = std::complex<double>(pars(433), pars(434));
   ZUR(2,2) = std::complex<double>(pars(435), pars(436));
   ZZ(0,0) = pars(437);
   ZZ(0,1) = pars(438);
   ZZ(1,0) = pars(439);
   ZZ(1,1) = pars(440);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(441);

   pars(64) = ZD(0,0);
   pars(65) = ZD(0,1);
   pars(66) = ZD(0,2);
   pars(67) = ZD(0,3);
   pars(68) = ZD(0,4);
   pars(69) = ZD(0,5);
   pars(70) = ZD(1,0);
   pars(71) = ZD(1,1);
   pars(72) = ZD(1,2);
   pars(73) = ZD(1,3);
   pars(74) = ZD(1,4);
   pars(75) = ZD(1,5);
   pars(76) = ZD(2,0);
   pars(77) = ZD(2,1);
   pars(78) = ZD(2,2);
   pars(79) = ZD(2,3);
   pars(80) = ZD(2,4);
   pars(81) = ZD(2,5);
   pars(82) = ZD(3,0);
   pars(83) = ZD(3,1);
   pars(84) = ZD(3,2);
   pars(85) = ZD(3,3);
   pars(86) = ZD(3,4);
   pars(87) = ZD(3,5);
   pars(88) = ZD(4,0);
   pars(89) = ZD(4,1);
   pars(90) = ZD(4,2);
   pars(91) = ZD(4,3);
   pars(92) = ZD(4,4);
   pars(93) = ZD(4,5);
   pars(94) = ZD(5,0);
   pars(95) = ZD(5,1);
   pars(96) = ZD(5,2);
   pars(97) = ZD(5,3);
   pars(98) = ZD(5,4);
   pars(99) = ZD(5,5);
   pars(100) = ZV(0,0);
   pars(101) = ZV(0,1);
   pars(102) = ZV(0,2);
   pars(103) = ZV(1,0);
   pars(104) = ZV(1,1);
   pars(105) = ZV(1,2);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZU(0,0);
   pars(110) = ZU(0,1);
   pars(111) = ZU(0,2);
   pars(112) = ZU(0,3);
   pars(113) = ZU(0,4);
   pars(114) = ZU(0,5);
   pars(115) = ZU(1,0);
   pars(116) = ZU(1,1);
   pars(117) = ZU(1,2);
   pars(118) = ZU(1,3);
   pars(119) = ZU(1,4);
   pars(120) = ZU(1,5);
   pars(121) = ZU(2,0);
   pars(122) = ZU(2,1);
   pars(123) = ZU(2,2);
   pars(124) = ZU(2,3);
   pars(125) = ZU(2,4);
   pars(126) = ZU(2,5);
   pars(127) = ZU(3,0);
   pars(128) = ZU(3,1);
   pars(129) = ZU(3,2);
   pars(130) = ZU(3,3);
   pars(131) = ZU(3,4);
   pars(132) = ZU(3,5);
   pars(133) = ZU(4,0);
   pars(134) = ZU(4,1);
   pars(135) = ZU(4,2);
   pars(136) = ZU(4,3);
   pars(137) = ZU(4,4);
   pars(138) = ZU(4,5);
   pars(139) = ZU(5,0);
   pars(140) = ZU(5,1);
   pars(141) = ZU(5,2);
   pars(142) = ZU(5,3);
   pars(143) = ZU(5,4);
   pars(144) = ZU(5,5);
   pars(145) = ZE(0,0);
   pars(146) = ZE(0,1);
   pars(147) = ZE(0,2);
   pars(148) = ZE(0,3);
   pars(149) = ZE(0,4);
   pars(150) = ZE(0,5);
   pars(151) = ZE(1,0);
   pars(152) = ZE(1,1);
   pars(153) = ZE(1,2);
   pars(154) = ZE(1,3);
   pars(155) = ZE(1,4);
   pars(156) = ZE(1,5);
   pars(157) = ZE(2,0);
   pars(158) = ZE(2,1);
   pars(159) = ZE(2,2);
   pars(160) = ZE(2,3);
   pars(161) = ZE(2,4);
   pars(162) = ZE(2,5);
   pars(163) = ZE(3,0);
   pars(164) = ZE(3,1);
   pars(165) = ZE(3,2);
   pars(166) = ZE(3,3);
   pars(167) = ZE(3,4);
   pars(168) = ZE(3,5);
   pars(169) = ZE(4,0);
   pars(170) = ZE(4,1);
   pars(171) = ZE(4,2);
   pars(172) = ZE(4,3);
   pars(173) = ZE(4,4);
   pars(174) = ZE(4,5);
   pars(175) = ZE(5,0);
   pars(176) = ZE(5,1);
   pars(177) = ZE(5,2);
   pars(178) = ZE(5,3);
   pars(179) = ZE(5,4);
   pars(180) = ZE(5,5);
   pars(181) = ZH(0,0);
   pars(182) = ZH(0,1);
   pars(183) = ZH(0,2);
   pars(184) = ZH(0,3);
   pars(185) = ZH(1,0);
   pars(186) = ZH(1,1);
   pars(187) = ZH(1,2);
   pars(188) = ZH(1,3);
   pars(189) = ZH(2,0);
   pars(190) = ZH(2,1);
   pars(191) = ZH(2,2);
   pars(192) = ZH(2,3);
   pars(193) = ZH(3,0);
   pars(194) = ZH(3,1);
   pars(195) = ZH(3,2);
   pars(196) = ZH(3,3);
   pars(197) = ZA(0,0);
   pars(198) = ZA(0,1);
   pars(199) = ZA(0,2);
   pars(200) = ZA(0,3);
   pars(201) = ZA(1,0);
   pars(202) = ZA(1,1);
   pars(203) = ZA(1,2);
   pars(204) = ZA(1,3);
   pars(205) = ZA(2,0);
   pars(206) = ZA(2,1);
   pars(207) = ZA(2,2);
   pars(208) = ZA(2,3);
   pars(209) = ZA(3,0);
   pars(210) = ZA(3,1);
   pars(211) = ZA(3,2);
   pars(212) = ZA(3,3);
   pars(213) = ZHR(0,0);
   pars(214) = ZHR(0,1);
   pars(215) = ZHR(1,0);
   pars(216) = ZHR(1,1);
   pars(217) = ZP(0,0);
   pars(218) = ZP(0,1);
   pars(219) = ZP(0,2);
   pars(220) = ZP(0,3);
   pars(221) = ZP(1,0);
   pars(222) = ZP(1,1);
   pars(223) = ZP(1,2);
   pars(224) = ZP(1,3);
   pars(225) = ZP(2,0);
   pars(226) = ZP(2,1);
   pars(227) = ZP(2,2);
   pars(228) = ZP(2,3);
   pars(229) = ZP(3,0);
   pars(230) = ZP(3,1);
   pars(231) = ZP(3,2);
   pars(232) = ZP(3,3);
   pars(233) = Re(ZN1(0,0));
   pars(234) = Im(ZN1(0,0));
   pars(235) = Re(ZN1(0,1));
   pars(236) = Im(ZN1(0,1));
   pars(237) = Re(ZN1(0,2));
   pars(238) = Im(ZN1(0,2));
   pars(239) = Re(ZN1(0,3));
   pars(240) = Im(ZN1(0,3));
   pars(241) = Re(ZN1(1,0));
   pars(242) = Im(ZN1(1,0));
   pars(243) = Re(ZN1(1,1));
   pars(244) = Im(ZN1(1,1));
   pars(245) = Re(ZN1(1,2));
   pars(246) = Im(ZN1(1,2));
   pars(247) = Re(ZN1(1,3));
   pars(248) = Im(ZN1(1,3));
   pars(249) = Re(ZN1(2,0));
   pars(250) = Im(ZN1(2,0));
   pars(251) = Re(ZN1(2,1));
   pars(252) = Im(ZN1(2,1));
   pars(253) = Re(ZN1(2,2));
   pars(254) = Im(ZN1(2,2));
   pars(255) = Re(ZN1(2,3));
   pars(256) = Im(ZN1(2,3));
   pars(257) = Re(ZN1(3,0));
   pars(258) = Im(ZN1(3,0));
   pars(259) = Re(ZN1(3,1));
   pars(260) = Im(ZN1(3,1));
   pars(261) = Re(ZN1(3,2));
   pars(262) = Im(ZN1(3,2));
   pars(263) = Re(ZN1(3,3));
   pars(264) = Im(ZN1(3,3));
   pars(265) = Re(ZN2(0,0));
   pars(266) = Im(ZN2(0,0));
   pars(267) = Re(ZN2(0,1));
   pars(268) = Im(ZN2(0,1));
   pars(269) = Re(ZN2(0,2));
   pars(270) = Im(ZN2(0,2));
   pars(271) = Re(ZN2(0,3));
   pars(272) = Im(ZN2(0,3));
   pars(273) = Re(ZN2(1,0));
   pars(274) = Im(ZN2(1,0));
   pars(275) = Re(ZN2(1,1));
   pars(276) = Im(ZN2(1,1));
   pars(277) = Re(ZN2(1,2));
   pars(278) = Im(ZN2(1,2));
   pars(279) = Re(ZN2(1,3));
   pars(280) = Im(ZN2(1,3));
   pars(281) = Re(ZN2(2,0));
   pars(282) = Im(ZN2(2,0));
   pars(283) = Re(ZN2(2,1));
   pars(284) = Im(ZN2(2,1));
   pars(285) = Re(ZN2(2,2));
   pars(286) = Im(ZN2(2,2));
   pars(287) = Re(ZN2(2,3));
   pars(288) = Im(ZN2(2,3));
   pars(289) = Re(ZN2(3,0));
   pars(290) = Im(ZN2(3,0));
   pars(291) = Re(ZN2(3,1));
   pars(292) = Im(ZN2(3,1));
   pars(293) = Re(ZN2(3,2));
   pars(294) = Im(ZN2(3,2));
   pars(295) = Re(ZN2(3,3));
   pars(296) = Im(ZN2(3,3));
   pars(297) = Re(UM1(0,0));
   pars(298) = Im(UM1(0,0));
   pars(299) = Re(UM1(0,1));
   pars(300) = Im(UM1(0,1));
   pars(301) = Re(UM1(1,0));
   pars(302) = Im(UM1(1,0));
   pars(303) = Re(UM1(1,1));
   pars(304) = Im(UM1(1,1));
   pars(305) = Re(UP1(0,0));
   pars(306) = Im(UP1(0,0));
   pars(307) = Re(UP1(0,1));
   pars(308) = Im(UP1(0,1));
   pars(309) = Re(UP1(1,0));
   pars(310) = Im(UP1(1,0));
   pars(311) = Re(UP1(1,1));
   pars(312) = Im(UP1(1,1));
   pars(313) = Re(UM2(0,0));
   pars(314) = Im(UM2(0,0));
   pars(315) = Re(UM2(0,1));
   pars(316) = Im(UM2(0,1));
   pars(317) = Re(UM2(1,0));
   pars(318) = Im(UM2(1,0));
   pars(319) = Re(UM2(1,1));
   pars(320) = Im(UM2(1,1));
   pars(321) = Re(UP2(0,0));
   pars(322) = Im(UP2(0,0));
   pars(323) = Re(UP2(0,1));
   pars(324) = Im(UP2(0,1));
   pars(325) = Re(UP2(1,0));
   pars(326) = Im(UP2(1,0));
   pars(327) = Re(UP2(1,1));
   pars(328) = Im(UP2(1,1));
   pars(329) = Re(ZEL(0,0));
   pars(330) = Im(ZEL(0,0));
   pars(331) = Re(ZEL(0,1));
   pars(332) = Im(ZEL(0,1));
   pars(333) = Re(ZEL(0,2));
   pars(334) = Im(ZEL(0,2));
   pars(335) = Re(ZEL(1,0));
   pars(336) = Im(ZEL(1,0));
   pars(337) = Re(ZEL(1,1));
   pars(338) = Im(ZEL(1,1));
   pars(339) = Re(ZEL(1,2));
   pars(340) = Im(ZEL(1,2));
   pars(341) = Re(ZEL(2,0));
   pars(342) = Im(ZEL(2,0));
   pars(343) = Re(ZEL(2,1));
   pars(344) = Im(ZEL(2,1));
   pars(345) = Re(ZEL(2,2));
   pars(346) = Im(ZEL(2,2));
   pars(347) = Re(ZER(0,0));
   pars(348) = Im(ZER(0,0));
   pars(349) = Re(ZER(0,1));
   pars(350) = Im(ZER(0,1));
   pars(351) = Re(ZER(0,2));
   pars(352) = Im(ZER(0,2));
   pars(353) = Re(ZER(1,0));
   pars(354) = Im(ZER(1,0));
   pars(355) = Re(ZER(1,1));
   pars(356) = Im(ZER(1,1));
   pars(357) = Re(ZER(1,2));
   pars(358) = Im(ZER(1,2));
   pars(359) = Re(ZER(2,0));
   pars(360) = Im(ZER(2,0));
   pars(361) = Re(ZER(2,1));
   pars(362) = Im(ZER(2,1));
   pars(363) = Re(ZER(2,2));
   pars(364) = Im(ZER(2,2));
   pars(365) = Re(ZDL(0,0));
   pars(366) = Im(ZDL(0,0));
   pars(367) = Re(ZDL(0,1));
   pars(368) = Im(ZDL(0,1));
   pars(369) = Re(ZDL(0,2));
   pars(370) = Im(ZDL(0,2));
   pars(371) = Re(ZDL(1,0));
   pars(372) = Im(ZDL(1,0));
   pars(373) = Re(ZDL(1,1));
   pars(374) = Im(ZDL(1,1));
   pars(375) = Re(ZDL(1,2));
   pars(376) = Im(ZDL(1,2));
   pars(377) = Re(ZDL(2,0));
   pars(378) = Im(ZDL(2,0));
   pars(379) = Re(ZDL(2,1));
   pars(380) = Im(ZDL(2,1));
   pars(381) = Re(ZDL(2,2));
   pars(382) = Im(ZDL(2,2));
   pars(383) = Re(ZDR(0,0));
   pars(384) = Im(ZDR(0,0));
   pars(385) = Re(ZDR(0,1));
   pars(386) = Im(ZDR(0,1));
   pars(387) = Re(ZDR(0,2));
   pars(388) = Im(ZDR(0,2));
   pars(389) = Re(ZDR(1,0));
   pars(390) = Im(ZDR(1,0));
   pars(391) = Re(ZDR(1,1));
   pars(392) = Im(ZDR(1,1));
   pars(393) = Re(ZDR(1,2));
   pars(394) = Im(ZDR(1,2));
   pars(395) = Re(ZDR(2,0));
   pars(396) = Im(ZDR(2,0));
   pars(397) = Re(ZDR(2,1));
   pars(398) = Im(ZDR(2,1));
   pars(399) = Re(ZDR(2,2));
   pars(400) = Im(ZDR(2,2));
   pars(401) = Re(ZUL(0,0));
   pars(402) = Im(ZUL(0,0));
   pars(403) = Re(ZUL(0,1));
   pars(404) = Im(ZUL(0,1));
   pars(405) = Re(ZUL(0,2));
   pars(406) = Im(ZUL(0,2));
   pars(407) = Re(ZUL(1,0));
   pars(408) = Im(ZUL(1,0));
   pars(409) = Re(ZUL(1,1));
   pars(410) = Im(ZUL(1,1));
   pars(411) = Re(ZUL(1,2));
   pars(412) = Im(ZUL(1,2));
   pars(413) = Re(ZUL(2,0));
   pars(414) = Im(ZUL(2,0));
   pars(415) = Re(ZUL(2,1));
   pars(416) = Im(ZUL(2,1));
   pars(417) = Re(ZUL(2,2));
   pars(418) = Im(ZUL(2,2));
   pars(419) = Re(ZUR(0,0));
   pars(420) = Im(ZUR(0,0));
   pars(421) = Re(ZUR(0,1));
   pars(422) = Im(ZUR(0,1));
   pars(423) = Re(ZUR(0,2));
   pars(424) = Im(ZUR(0,2));
   pars(425) = Re(ZUR(1,0));
   pars(426) = Im(ZUR(1,0));
   pars(427) = Re(ZUR(1,1));
   pars(428) = Im(ZUR(1,1));
   pars(429) = Re(ZUR(1,2));
   pars(430) = Im(ZUR(1,2));
   pars(431) = Re(ZUR(2,0));
   pars(432) = Im(ZUR(2,0));
   pars(433) = Re(ZUR(2,1));
   pars(434) = Im(ZUR(2,1));
   pars(435) = Re(ZUR(2,2));
   pars(436) = Im(ZUR(2,2));
   pars(437) = ZZ(0,0);
   pars(438) = ZZ(0,1);
   pars(439) = ZZ(1,0);
   pars(440) = ZZ(1,1);


   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

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

   const double mass_matrix_SRdp = Re(0.125*(8*mRd2 + 3.0983866769659336*g1*
      MDBS*vS + 4*g2*MDWBT*vT + 8*AbsSqr(MuD) - 2*vS*(-2.8284271247461903*MuD -
      2*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamSD) +
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
      problems.flag_running_tachyon(MRSSM2_info::SRdp);
   }

   MSRdp = AbsSqrt(MSRdp);
}

double CLASSNAME::get_mass_matrix_SRum() const
{

   const double mass_matrix_SRum = Re(0.125*(8*mRu2 - 3.0983866769659336*g1*
      MDBS*vS - 4*g2*MDWBT*vT + 8*AbsSqr(MuU) + 2*vS*(2.8284271247461903*MuU +
      2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamSU) -
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
      problems.flag_running_tachyon(MRSSM2_info::SRum);
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
      problems.flag_running_tachyon(MRSSM2_info::sigmaO);
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
      problems.flag_running_tachyon(MRSSM2_info::phiO);
   }

   MphiO = AbsSqrt(MphiO);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(0,0) +
      0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) + AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) + Conj(
      Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) + Conj(
      Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = -0.7071067811865475*vu*Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = -0.7071067811865475*vu*Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = -0.7071067811865475*vu*Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,1) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(1,1) +
      0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) + AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) + Conj(
      Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = -0.7071067811865475*vu*Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = -0.7071067811865475*vu*Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = -0.7071067811865475*vu*Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,2) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(2,2) +
      0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = -0.7071067811865475*vu*Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = -0.7071067811865475*vu*Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = -0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,3) = 0.2581988897471611*g1*MDBS*vS + 0.2581988897471611*g1*
      vS*Conj(MDBS) + md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)) +
      AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) + Conj(
      Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) + Conj(
      Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = 0.2581988897471611*g1*MDBS*vS + 0.2581988897471611*g1*
      vS*Conj(MDBS) + md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) + Conj(
      Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = 0.2581988897471611*g1*MDBS*vS + 0.2581988897471611*g1*
      vS*Conj(MDBS) + md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(0,0) +
      0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(1,1) +
      0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(2,2) +
      0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(0,0) -
      0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,0)) +
      AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) + Conj(
      Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) + Conj(
      Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = -0.7071067811865475*vd*Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = -0.7071067811865475*vd*Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = -0.7071067811865475*vd*Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,1) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(1,1) -
      0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,1)) +
      AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = -0.7071067811865475*vd*Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = -0.7071067811865475*vd*Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = -0.7071067811865475*vd*Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,2) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT +
      0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(2,2) -
      0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,2)) +
      AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = -0.7071067811865475*vd*Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = -0.7071067811865475*vd*Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = -0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,3) = -0.5163977794943222*g1*MDBS*vS - 0.5163977794943222*g1
      *vS*Conj(MDBS) + mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,0)) +
      AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) + Conj(
      Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) + Conj(
      Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = -0.5163977794943222*g1*MDBS*vS - 0.5163977794943222*g1
      *vS*Conj(MDBS) + mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(1,0)) +
      AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) + Conj(
      Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = -0.5163977794943222*g1*MDBS*vS - 0.5163977794943222*g1
      *vS*Conj(MDBS) + mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(2,0)) +
      AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(0,0) +
      0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) + AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) + Conj(
      Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) + Conj(
      Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = -0.7071067811865475*vu*Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = -0.7071067811865475*vu*Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = -0.7071067811865475*vu*Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,1) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(1,1) +
      0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) + AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) + Conj(
      Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = -0.7071067811865475*vu*Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = -0.7071067811865475*vu*Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = -0.7071067811865475*vu*Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,2) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT -
      0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(2,2) +
      0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) + AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075
      *Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*
      Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = -0.7071067811865475*vu*Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = -0.7071067811865475*vu*Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = -0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,3) = 0.7745966692414834*g1*MDBS*vS + 0.7745966692414834*g1*
      vS*Conj(MDBS) + me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)) +
      AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) + Conj(
      Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) + Conj(
      Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = 0.7745966692414834*g1*MDBS*vS + 0.7745966692414834*g1*
      vS*Conj(MDBS) + me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) + Conj(
      Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = 0.7745966692414834*g1*MDBS*vS + 0.7745966692414834*g1*
      vS*Conj(MDBS) + me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,4,4> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
       + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(LamSD) +
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(LamTD) +
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) + 0.5*LamTD*vT*Conj(MuD) + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) - 0.075*Sqr(g1
      )*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*Conj(BMu) - 0.15*vd*vu*Sqr(g1) - 0.25*
      vd*vu*Sqr(g2);
   mass_matrix_hh(0,2) = -0.3872983346207417*g1*MDBS*vd + vd*vS*AbsSqr(LamSD) +
      0.7071067811865475*MuD*vd*Conj(LamSD) + 0.35355339059327373*LamTD*vd*vT*
      Conj(LamSD) + 0.35355339059327373*LamSD*vd*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*Conj(MDBS) + 0.7071067811865475*LamSD*vd*Conj(
      MuD);
   mass_matrix_hh(0,3) = 0.5*g2*MDWBT*vd + 0.5*vd*vT*AbsSqr(LamTD) +
      0.35355339059327373*LamTD*vd*vS*Conj(LamSD) + 0.5*MuD*vd*Conj(LamTD) +
      0.35355339059327373*LamSD*vd*vS*Conj(LamTD) + 0.5*g2*vd*Conj(MDWBT) + 0.5
      *LamTD*vd*Conj(MuD);
   mass_matrix_hh(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
       + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(LamSU) -
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) - 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) + 0.225*Sqr(g1
      )*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = 0.3872983346207417*g1*MDBS*vu + vS*vu*AbsSqr(LamSU) +
      0.7071067811865475*MuU*vu*Conj(LamSU) - 0.35355339059327373*LamTU*vT*vu*
      Conj(LamSU) - 0.35355339059327373*LamSU*vT*vu*Conj(LamTU) +
      0.3872983346207417*g1*vu*Conj(MDBS) + 0.7071067811865475*LamSU*vu*Conj(
      MuU);
   mass_matrix_hh(1,3) = -0.5*g2*MDWBT*vu + 0.5*vT*vu*AbsSqr(LamTU) -
      0.35355339059327373*LamTU*vS*vu*Conj(LamSU) - 0.5*MuU*vu*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vu*Conj(LamTU) - 0.5*g2*vu*Conj(MDWBT) - 0.5
      *LamTU*vu*Conj(MuU);
   mass_matrix_hh(2,2) = mS2 + 4*Sqr(MDBS) + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.5*
      AbsSqr(LamSU)*Sqr(vu);
   mass_matrix_hh(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_hh(3,3) = mT2 + 4*Sqr(MDWBT) + 0.25*AbsSqr(LamTD)*Sqr(vd) + 0.25
      *AbsSqr(LamTU)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,4,4> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
       + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(LamSD) +
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(LamTD) +
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) + 0.5*LamTD*vT*Conj(MuD) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(
      ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*
      AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) - 0.075*Sqr(g1)*Sqr(vu
      ) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW())) +
      0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) - 0.3872983346207417*g1*g2*vd*
      vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW())) -
      0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,2) = 0;
   mass_matrix_Ah(0,3) = 0;
   mass_matrix_Ah(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
       + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(LamSU) -
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) - 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW
      ())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = 0;
   mass_matrix_Ah(1,3) = 0;
   mass_matrix_Ah(2,2) = mS2 + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.5*AbsSqr(LamSU)*
      Sqr(vu);
   mass_matrix_Ah(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_Ah(3,3) = mT2 + 0.25*AbsSqr(LamTD)*Sqr(vd) + 0.25*AbsSqr(LamTU)*
      Sqr(vu);

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Rh() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Rh;

   mass_matrix_Rh(0,0) = mRd2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
       + AbsSqr(MuD) + 0.7071067811865475*MuD*vS*Conj(LamSD) +
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(LamTD) +
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) + 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) + 0.5*LamTD*vT*Conj(MuD) + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.25*AbsSqr(
      LamTD)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*
      AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) + 0.075*Sqr(g1)*Sqr(vu
      ) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Rh(0,1) = -0.5*LamSU*vd*vu*Conj(LamSD) + 0.25*LamTU*vd*vu*Conj(
      LamTD);
   mass_matrix_Rh(1,1) = mRu2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
       + AbsSqr(MuU) + 0.7071067811865475*MuU*vS*Conj(LamSU) -
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
   problems.flag_bad_mass(MRSSM2_info::Rh, eigenvalue_error > precision * Abs(
      MRh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Rh, MRh, ZHR);
#endif
   normalize_to_interval(ZHR);


   if (MRh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Rh);
   }

   MRh = AbsSqrt(MRh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,4,4> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*
      vT + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(LamSD) -
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) - 0.5*MuD*vT*Conj(LamTD) -
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) - 0.5*LamTD*vT*Conj(MuD) + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) - 0.075*Sqr(g1
      )*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = Conj(BMu);
   mass_matrix_Hpm(0,2) = 0.7071067811865475*g2*MDWBT*vd - 0.35355339059327373*
      vd*vT*AbsSqr(LamTD) + 0.5*LamTD*vd*vS*Conj(LamSD) + 0.7071067811865475*
      LamTD*vd*Conj(MuD) + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(0,3) = 0.35355339059327373*vd*vT*AbsSqr(LamTD) +
      0.7071067811865475*MuD*vd*Conj(LamTD) + 0.5*LamSD*vd*vS*Conj(LamTD) +
      0.7071067811865475*g2*vd*Conj(MDWBT);
   mass_matrix_Hpm(1,0) = BMu;
   mass_matrix_Hpm(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*
      vT + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(LamSU) +
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) + 0.5*MuU*vT*Conj(LamTU) +
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) + 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) + 0.075*Sqr(g1
      )*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(1,2) = 0.7071067811865475*g2*MDWBT*vu - 0.35355339059327373*
      vT*vu*AbsSqr(LamTU) + 0.5*LamTU*vS*vu*Conj(LamSU) + 0.7071067811865475*
      LamTU*vu*Conj(MuU);
   mass_matrix_Hpm(1,3) = 0.35355339059327373*vT*vu*AbsSqr(LamTU) +
      0.7071067811865475*MuU*vu*Conj(LamTU) + 0.5*LamSU*vS*vu*Conj(LamTU) +
      0.7071067811865475*g2*vu*Conj(MDWBT) - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(2,0) = -0.35355339059327373*vd*vT*AbsSqr(LamTD) +
      0.7071067811865475*MuD*vd*Conj(LamTD) + 0.5*LamSD*vd*vS*Conj(LamTD) +
      0.7071067811865475*g2*vd*Conj(MDWBT) + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(2,1) = -0.35355339059327373*vT*vu*AbsSqr(LamTU) +
      0.7071067811865475*MuU*vu*Conj(LamTU) + 0.5*LamSU*vS*vu*Conj(LamTU) +
      0.7071067811865475*g2*vu*Conj(MDWBT);
   mass_matrix_Hpm(2,2) = mT2 + 2*Sqr(MDWBT) + 0.5*AbsSqr(LamTD)*Sqr(vd) - 0.25
      *Sqr(g2)*Sqr(vd) + Sqr(g2)*Sqr(vT) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(2,3) = 2*Sqr(MDWBT);
   mass_matrix_Hpm(3,0) = 0.7071067811865475*g2*MDWBT*vd + 0.35355339059327373*
      vd*vT*AbsSqr(LamTD) + 0.5*LamTD*vd*vS*Conj(LamSD) + 0.7071067811865475*
      LamTD*vd*Conj(MuD);
   mass_matrix_Hpm(3,1) = 0.7071067811865475*g2*MDWBT*vu + 0.35355339059327373*
      vT*vu*AbsSqr(LamTU) + 0.5*LamTU*vS*vu*Conj(LamSU) + 0.7071067811865475*
      LamTU*vu*Conj(MuU) - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(3,2) = 2*Sqr(MDWBT);
   mass_matrix_Hpm(3,3) = mT2 + 2*Sqr(MDWBT) + 0.25*Sqr(g2)*Sqr(vd) + Sqr(g2)*
      Sqr(vT) + 0.5*AbsSqr(LamTU)*Sqr(vu) - 0.25*Sqr(g2)*Sqr(vu);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Hpm, eigenvalue_error > precision * Abs(
      MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::Hpm);
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
   mass_matrix_Chi(2,2) = -MuD - 0.7071067811865475*LamSD*vS - 0.5*LamTD*vT;
   mass_matrix_Chi(2,3) = 0;
   mass_matrix_Chi(3,0) = 0.7071067811865475*LamSU*vu;
   mass_matrix_Chi(3,1) = -0.5*LamTU*vu;
   mass_matrix_Chi(3,2) = 0;
   mass_matrix_Chi(3,3) = MuU + 0.7071067811865475*LamSU*vS - 0.5*LamTU*vT;

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Chi, MChi, ZN1, ZN2, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Chi, eigenvalue_error > precision * Abs(
      MChi(0)));
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
   mass_matrix_Cha1(1,1) = MuD + 0.7071067811865475*LamSD*vS - 0.5*LamTD*vT;

   return mass_matrix_Cha1;
}

void CLASSNAME::calculate_MCha1()
{
   const auto mass_matrix_Cha1(get_mass_matrix_Cha1());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha1, MCha1, UM1, UP1, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Cha1, eigenvalue_error > precision * Abs
      (MCha1(0)));
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
   mass_matrix_Cha2(1,1) = -MuU - 0.7071067811865475*LamSU*vS - 0.5*LamTU*vT;

   return mass_matrix_Cha2;
}

void CLASSNAME::calculate_MCha2()
{
   const auto mass_matrix_Cha2(get_mass_matrix_Cha2());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha2, MCha2, UM2, UP2, eigenvalue_error);
   problems.flag_bad_mass(MRSSM2_info::Cha2, eigenvalue_error > precision * Abs
      (MCha2(0)));
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
   problems.flag_bad_mass(MRSSM2_info::Fe, eigenvalue_error > precision * Abs(
      MFe(0)));
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
   problems.flag_bad_mass(MRSSM2_info::Fd, eigenvalue_error > precision * Abs(
      MFd(0)));
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
   problems.flag_bad_mass(MRSSM2_info::Fu, eigenvalue_error > precision * Abs(
      MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + 4*Sqr(vT) + Sqr(vu
      )));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(MRSSM2_info::VWm);
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
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error);
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   
   double result = Re(mHd2*vd - 0.3872983346207417*g1*MDBS*vd*vS + 0.5*g2*MDWBT*vd
      *vT + vd*AbsSqr(MuD) + vd*AbsSqr(Mu) - 0.5*vu*BMu + 0.7071067811865475*MuD*
      vd*vS*Conj(LamSD) + 0.35355339059327373*LamTD*vd*vS*vT*Conj(LamSD) + 0.5*MuD
      *vd*vT*Conj(LamTD) + 0.35355339059327373*LamSD*vd*vS*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*vS*Conj(MDBS) + 0.5*g2*vd*vT*Conj(MDWBT) +
      0.7071067811865475*LamSD*vd*vS*Conj(MuD) + 0.5*LamTD*vd*vT*Conj(MuD) - 0.5*
      vu*Conj(BMu) + 0.075*Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.5*vd*
      AbsSqr(LamSD)*Sqr(vS) + 0.25*vd*AbsSqr(LamTD)*Sqr(vT) - 0.075*vd*Sqr(g1)*Sqr
      (vu) - 0.125*vd*Sqr(g2)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu + 0.3872983346207417*g1*MDBS*vS*vu - 0.5*g2*MDWBT*vT
      *vu + vu*AbsSqr(MuU) + vu*AbsSqr(Mu) - 0.5*vd*BMu + 0.7071067811865475*MuU*
      vS*vu*Conj(LamSU) - 0.35355339059327373*LamTU*vS*vT*vu*Conj(LamSU) - 0.5*MuU
      *vT*vu*Conj(LamTU) - 0.35355339059327373*LamSU*vS*vT*vu*Conj(LamTU) +
      0.3872983346207417*g1*vS*vu*Conj(MDBS) - 0.5*g2*vT*vu*Conj(MDWBT) +
      0.7071067811865475*LamSU*vS*vu*Conj(MuU) - 0.5*LamTU*vT*vu*Conj(MuU) - 0.5*
      vd*Conj(BMu) + 0.075*Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) - 0.075*vu*
      Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(LamSU)*Sqr(vS) +
      0.25*vu*AbsSqr(LamTU)*Sqr(vT));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(mT2*vT + 4*vT*Sqr(MDWBT) + 0.25*g2*MDWBT*Sqr(vd) + 0.25*vT*
      AbsSqr(LamTD)*Sqr(vd) + 0.17677669529663687*LamTD*vS*Conj(LamSD)*Sqr(vd) +
      0.25*MuD*Conj(LamTD)*Sqr(vd) + 0.17677669529663687*LamSD*vS*Conj(LamTD)*Sqr(
      vd) + 0.25*g2*Conj(MDWBT)*Sqr(vd) + 0.25*LamTD*Conj(MuD)*Sqr(vd) - 0.25*g2*
      MDWBT*Sqr(vu) + 0.25*vT*AbsSqr(LamTU)*Sqr(vu) - 0.17677669529663687*LamTU*vS
      *Conj(LamSU)*Sqr(vu) - 0.25*MuU*Conj(LamTU)*Sqr(vu) - 0.17677669529663687*
      LamSU*vS*Conj(LamTU)*Sqr(vu) - 0.25*g2*Conj(MDWBT)*Sqr(vu) - 0.25*LamTU*Conj
      (MuU)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_4() const
{
   
   double result = Re(mS2*vS + 4*vS*Sqr(MDBS) - 0.19364916731037085*g1*MDBS*Sqr(vd
      ) + 0.5*vS*AbsSqr(LamSD)*Sqr(vd) + 0.35355339059327373*MuD*Conj(LamSD)*Sqr(
      vd) + 0.17677669529663687*LamTD*vT*Conj(LamSD)*Sqr(vd) + 0.17677669529663687
      *LamSD*vT*Conj(LamTD)*Sqr(vd) - 0.19364916731037085*g1*Conj(MDBS)*Sqr(vd) +
      0.35355339059327373*LamSD*Conj(MuD)*Sqr(vd) + 0.19364916731037085*g1*MDBS*
      Sqr(vu) + 0.5*vS*AbsSqr(LamSU)*Sqr(vu) + 0.35355339059327373*MuU*Conj(LamSU)
      *Sqr(vu) - 0.17677669529663687*LamTU*vT*Conj(LamSU)*Sqr(vu) -
      0.17677669529663687*LamSU*vT*Conj(LamTU)*Sqr(vu) + 0.19364916731037085*g1*
      Conj(MDBS)*Sqr(vu) + 0.35355339059327373*LamSU*Conj(MuU)*Sqr(vu));

   return result;
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

double CLASSNAME::VEV() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}



std::ostream& operator<<(std::ostream& ostr, const MRSSM2_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
