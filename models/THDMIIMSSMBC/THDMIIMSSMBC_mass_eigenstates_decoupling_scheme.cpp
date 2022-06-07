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
 * @file THDMIIMSSMBC_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the THDMIIMSSMBC model class in the decoupling scheme
 *
 * Contains the definition of the THDMIIMSSMBC model class methods
 * which solve EWSB and calculate masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "THDMIIMSSMBC_mass_eigenstates_decoupling_scheme.hpp"
#include "THDMIIMSSMBC_mass_eigenstates.hpp"
#include "THDMIIMSSMBC_info.hpp"
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

#define CLASSNAME THDMIIMSSMBC_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const THDMIIMSSMBC_input_parameters& input_)
   : THDMIIMSSMBC_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const THDMIIMSSMBC_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<THDMIIMSSMBC_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<THDMIIMSSMBC_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            THDMIIMSSMBC_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            THDMIIMSSMBC_info::g2, new_g2, get_scale());
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

      MODEL->set_v1(Re((2*MZMSbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
         ))));
      MODEL->set_v2(Re((2*MZMSbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
         TanBeta)))));

   }

   // apply user-defined low-energy constraint for the Yukawa couplings
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v2 = MODELPARAMETER(v2);
      MODEL->set_Yu((-((1.4142135623730951*upQuarksDRbar)/v2).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v1 = MODELPARAMETER(v1);
      MODEL->set_Yd((((1.4142135623730951*downQuarksDRbar)/v1).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v1 = MODELPARAMETER(v1);
      MODEL->set_Ye((((1.4142135623730951*downLeptonsDRbar)/v1).transpose()).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const THDMIIMSSMBC_mass_eigenstates& model)
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
   MFv = OTHER(MFv);
   Mhh = OTHER(Mhh);
   ZH = OTHER(ZH);
   MAh = OTHER(MAh);
   ZA = OTHER(ZA);
   MHm = OTHER(MHm);
   ZP = OTHER(ZP);
   MFd = OTHER(MFd);
   Vd = OTHER(Vd);
   Ud = OTHER(Ud);
   MFu = OTHER(MFu);
   Vu = OTHER(Vu);
   Uu = OTHER(Uu);
   MFe = OTHER(MFe);
   Ve = OTHER(Ve);
   Ue = OTHER(Ue);
   MVWm = OTHER(MVWm);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const THDMIIMSSMBC_physical& CLASSNAME::get_physical() const
{
   return physical;
}

THDMIIMSSMBC_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const THDMIIMSSMBC_physical& physical_)
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

   
   const double old_M112 = M112;
   const double old_M222 = M222;

   M112 = Re((0.25*(2*M122*v2 + 2*v2*Conj(M122) - 4*Lambda1*Cube(v1) - Lambda7*
      Cube(v2) - Conj(Lambda7)*Cube(v2) - 3*Lambda6*v2*Sqr(v1) - 3*v2*Conj(Lambda6
      )*Sqr(v1) - 2*Lambda3*v1*Sqr(v2) - 2*Lambda4*v1*Sqr(v2) - Lambda5*v1*Sqr(v2)
      - v1*Conj(Lambda5)*Sqr(v2)))/v1);
   M222 = Re((0.25*(2*M122*v1 + 2*v1*Conj(M122) - Lambda6*Cube(v1) - Conj(Lambda6)
      *Cube(v1) - 4*Lambda2*Cube(v2) - 2*Lambda3*v2*Sqr(v1) - 2*Lambda4*v2*Sqr(v1)
      - Lambda5*v2*Sqr(v1) - v2*Conj(Lambda5)*Sqr(v1) - 3*Lambda7*v1*Sqr(v2) - 3*
      v1*Conj(Lambda7)*Sqr(v2)))/v2);

   const bool is_finite = IsFinite(M112) && IsFinite(M222);

   if (!is_finite) {
      M112 = old_M112;
      M222 = old_M222;
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
           "THDMIIMSSMBC\n"
           "========================================\n";
   THDMIIMSSMBC_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_tree_level_mass_spectrum()
{
   const auto save_M112_raii = make_raii_save(M112);
   const auto save_M222_raii = make_raii_save(M222);

   solve_ewsb_equations_tree_level();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MHm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MFv();
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
   PHYSICAL(MFv) = MFv;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHm) = MHm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWm) = MVWm;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(ZZ) = ZZ;

}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_tree_level_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHm, ZP);

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
   move_goldstone_to(0, MVWm, PHYSICAL(MHm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(THDMIIMSSMBC_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(THDMIIMSSMBC_info::Ah); }
   if (PHYSICAL(MHm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(THDMIIMSSMBC_info::Hm); }
}

/**
 * calculates spectrum for model once the MSbar parameters at
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
   MFv = Eigen::Matrix<double,3,1>::Zero();
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
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
   THDMIIMSSMBC_soft_parameters::clear();
   clear_tree_level_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MFv(0) = pars(1);
   MFv(1) = pars(2);
   MFv(2) = pars(3);
   Mhh(0) = pars(4);
   Mhh(1) = pars(5);
   MAh(0) = pars(6);
   MAh(1) = pars(7);
   MHm(0) = pars(8);
   MHm(1) = pars(9);
   MFd(0) = pars(10);
   MFd(1) = pars(11);
   MFd(2) = pars(12);
   MFu(0) = pars(13);
   MFu(1) = pars(14);
   MFu(2) = pars(15);
   MFe(0) = pars(16);
   MFe(1) = pars(17);
   MFe(2) = pars(18);
   MVWm = pars(19);
   MVP = pars(20);
   MVZ = pars(21);

}

const THDMIIMSSMBC_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

THDMIIMSSMBC_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(22);

   pars(0) = MVG;
   pars(1) = MFv(0);
   pars(2) = MFv(1);
   pars(3) = MFv(2);
   pars(4) = Mhh(0);
   pars(5) = Mhh(1);
   pars(6) = MAh(0);
   pars(7) = MAh(1);
   pars(8) = MHm(0);
   pars(9) = MHm(1);
   pars(10) = MFd(0);
   pars(11) = MFd(1);
   pars(12) = MFd(2);
   pars(13) = MFu(0);
   pars(14) = MFu(1);
   pars(15) = MFu(2);
   pars(16) = MFe(0);
   pars(17) = MFe(1);
   pars(18) = MFe(2);
   pars(19) = MVWm;
   pars(20) = MVP;
   pars(21) = MVZ;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZH(0,0) = pars(22);
   ZH(0,1) = pars(23);
   ZH(1,0) = pars(24);
   ZH(1,1) = pars(25);
   ZA(0,0) = pars(26);
   ZA(0,1) = pars(27);
   ZA(1,0) = pars(28);
   ZA(1,1) = pars(29);
   ZP(0,0) = pars(30);
   ZP(0,1) = pars(31);
   ZP(1,0) = pars(32);
   ZP(1,1) = pars(33);
   Vd(0,0) = std::complex<double>(pars(34), pars(35));
   Vd(0,1) = std::complex<double>(pars(36), pars(37));
   Vd(0,2) = std::complex<double>(pars(38), pars(39));
   Vd(1,0) = std::complex<double>(pars(40), pars(41));
   Vd(1,1) = std::complex<double>(pars(42), pars(43));
   Vd(1,2) = std::complex<double>(pars(44), pars(45));
   Vd(2,0) = std::complex<double>(pars(46), pars(47));
   Vd(2,1) = std::complex<double>(pars(48), pars(49));
   Vd(2,2) = std::complex<double>(pars(50), pars(51));
   Ud(0,0) = std::complex<double>(pars(52), pars(53));
   Ud(0,1) = std::complex<double>(pars(54), pars(55));
   Ud(0,2) = std::complex<double>(pars(56), pars(57));
   Ud(1,0) = std::complex<double>(pars(58), pars(59));
   Ud(1,1) = std::complex<double>(pars(60), pars(61));
   Ud(1,2) = std::complex<double>(pars(62), pars(63));
   Ud(2,0) = std::complex<double>(pars(64), pars(65));
   Ud(2,1) = std::complex<double>(pars(66), pars(67));
   Ud(2,2) = std::complex<double>(pars(68), pars(69));
   Vu(0,0) = std::complex<double>(pars(70), pars(71));
   Vu(0,1) = std::complex<double>(pars(72), pars(73));
   Vu(0,2) = std::complex<double>(pars(74), pars(75));
   Vu(1,0) = std::complex<double>(pars(76), pars(77));
   Vu(1,1) = std::complex<double>(pars(78), pars(79));
   Vu(1,2) = std::complex<double>(pars(80), pars(81));
   Vu(2,0) = std::complex<double>(pars(82), pars(83));
   Vu(2,1) = std::complex<double>(pars(84), pars(85));
   Vu(2,2) = std::complex<double>(pars(86), pars(87));
   Uu(0,0) = std::complex<double>(pars(88), pars(89));
   Uu(0,1) = std::complex<double>(pars(90), pars(91));
   Uu(0,2) = std::complex<double>(pars(92), pars(93));
   Uu(1,0) = std::complex<double>(pars(94), pars(95));
   Uu(1,1) = std::complex<double>(pars(96), pars(97));
   Uu(1,2) = std::complex<double>(pars(98), pars(99));
   Uu(2,0) = std::complex<double>(pars(100), pars(101));
   Uu(2,1) = std::complex<double>(pars(102), pars(103));
   Uu(2,2) = std::complex<double>(pars(104), pars(105));
   Ve(0,0) = std::complex<double>(pars(106), pars(107));
   Ve(0,1) = std::complex<double>(pars(108), pars(109));
   Ve(0,2) = std::complex<double>(pars(110), pars(111));
   Ve(1,0) = std::complex<double>(pars(112), pars(113));
   Ve(1,1) = std::complex<double>(pars(114), pars(115));
   Ve(1,2) = std::complex<double>(pars(116), pars(117));
   Ve(2,0) = std::complex<double>(pars(118), pars(119));
   Ve(2,1) = std::complex<double>(pars(120), pars(121));
   Ve(2,2) = std::complex<double>(pars(122), pars(123));
   Ue(0,0) = std::complex<double>(pars(124), pars(125));
   Ue(0,1) = std::complex<double>(pars(126), pars(127));
   Ue(0,2) = std::complex<double>(pars(128), pars(129));
   Ue(1,0) = std::complex<double>(pars(130), pars(131));
   Ue(1,1) = std::complex<double>(pars(132), pars(133));
   Ue(1,2) = std::complex<double>(pars(134), pars(135));
   Ue(2,0) = std::complex<double>(pars(136), pars(137));
   Ue(2,1) = std::complex<double>(pars(138), pars(139));
   Ue(2,2) = std::complex<double>(pars(140), pars(141));
   ZZ(0,0) = pars(142);
   ZZ(0,1) = pars(143);
   ZZ(1,0) = pars(144);
   ZZ(1,1) = pars(145);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(146);

   pars(22) = ZH(0,0);
   pars(23) = ZH(0,1);
   pars(24) = ZH(1,0);
   pars(25) = ZH(1,1);
   pars(26) = ZA(0,0);
   pars(27) = ZA(0,1);
   pars(28) = ZA(1,0);
   pars(29) = ZA(1,1);
   pars(30) = ZP(0,0);
   pars(31) = ZP(0,1);
   pars(32) = ZP(1,0);
   pars(33) = ZP(1,1);
   pars(34) = Re(Vd(0,0));
   pars(35) = Im(Vd(0,0));
   pars(36) = Re(Vd(0,1));
   pars(37) = Im(Vd(0,1));
   pars(38) = Re(Vd(0,2));
   pars(39) = Im(Vd(0,2));
   pars(40) = Re(Vd(1,0));
   pars(41) = Im(Vd(1,0));
   pars(42) = Re(Vd(1,1));
   pars(43) = Im(Vd(1,1));
   pars(44) = Re(Vd(1,2));
   pars(45) = Im(Vd(1,2));
   pars(46) = Re(Vd(2,0));
   pars(47) = Im(Vd(2,0));
   pars(48) = Re(Vd(2,1));
   pars(49) = Im(Vd(2,1));
   pars(50) = Re(Vd(2,2));
   pars(51) = Im(Vd(2,2));
   pars(52) = Re(Ud(0,0));
   pars(53) = Im(Ud(0,0));
   pars(54) = Re(Ud(0,1));
   pars(55) = Im(Ud(0,1));
   pars(56) = Re(Ud(0,2));
   pars(57) = Im(Ud(0,2));
   pars(58) = Re(Ud(1,0));
   pars(59) = Im(Ud(1,0));
   pars(60) = Re(Ud(1,1));
   pars(61) = Im(Ud(1,1));
   pars(62) = Re(Ud(1,2));
   pars(63) = Im(Ud(1,2));
   pars(64) = Re(Ud(2,0));
   pars(65) = Im(Ud(2,0));
   pars(66) = Re(Ud(2,1));
   pars(67) = Im(Ud(2,1));
   pars(68) = Re(Ud(2,2));
   pars(69) = Im(Ud(2,2));
   pars(70) = Re(Vu(0,0));
   pars(71) = Im(Vu(0,0));
   pars(72) = Re(Vu(0,1));
   pars(73) = Im(Vu(0,1));
   pars(74) = Re(Vu(0,2));
   pars(75) = Im(Vu(0,2));
   pars(76) = Re(Vu(1,0));
   pars(77) = Im(Vu(1,0));
   pars(78) = Re(Vu(1,1));
   pars(79) = Im(Vu(1,1));
   pars(80) = Re(Vu(1,2));
   pars(81) = Im(Vu(1,2));
   pars(82) = Re(Vu(2,0));
   pars(83) = Im(Vu(2,0));
   pars(84) = Re(Vu(2,1));
   pars(85) = Im(Vu(2,1));
   pars(86) = Re(Vu(2,2));
   pars(87) = Im(Vu(2,2));
   pars(88) = Re(Uu(0,0));
   pars(89) = Im(Uu(0,0));
   pars(90) = Re(Uu(0,1));
   pars(91) = Im(Uu(0,1));
   pars(92) = Re(Uu(0,2));
   pars(93) = Im(Uu(0,2));
   pars(94) = Re(Uu(1,0));
   pars(95) = Im(Uu(1,0));
   pars(96) = Re(Uu(1,1));
   pars(97) = Im(Uu(1,1));
   pars(98) = Re(Uu(1,2));
   pars(99) = Im(Uu(1,2));
   pars(100) = Re(Uu(2,0));
   pars(101) = Im(Uu(2,0));
   pars(102) = Re(Uu(2,1));
   pars(103) = Im(Uu(2,1));
   pars(104) = Re(Uu(2,2));
   pars(105) = Im(Uu(2,2));
   pars(106) = Re(Ve(0,0));
   pars(107) = Im(Ve(0,0));
   pars(108) = Re(Ve(0,1));
   pars(109) = Im(Ve(0,1));
   pars(110) = Re(Ve(0,2));
   pars(111) = Im(Ve(0,2));
   pars(112) = Re(Ve(1,0));
   pars(113) = Im(Ve(1,0));
   pars(114) = Re(Ve(1,1));
   pars(115) = Im(Ve(1,1));
   pars(116) = Re(Ve(1,2));
   pars(117) = Im(Ve(1,2));
   pars(118) = Re(Ve(2,0));
   pars(119) = Im(Ve(2,0));
   pars(120) = Re(Ve(2,1));
   pars(121) = Im(Ve(2,1));
   pars(122) = Re(Ve(2,2));
   pars(123) = Im(Ve(2,2));
   pars(124) = Re(Ue(0,0));
   pars(125) = Im(Ue(0,0));
   pars(126) = Re(Ue(0,1));
   pars(127) = Im(Ue(0,1));
   pars(128) = Re(Ue(0,2));
   pars(129) = Im(Ue(0,2));
   pars(130) = Re(Ue(1,0));
   pars(131) = Im(Ue(1,0));
   pars(132) = Re(Ue(1,1));
   pars(133) = Im(Ue(1,1));
   pars(134) = Re(Ue(1,2));
   pars(135) = Im(Ue(1,2));
   pars(136) = Re(Ue(2,0));
   pars(137) = Im(Ue(2,0));
   pars(138) = Re(Ue(2,1));
   pars(139) = Im(Ue(2,1));
   pars(140) = Re(Ue(2,2));
   pars(141) = Im(Ue(2,2));
   pars(142) = ZZ(0,0);
   pars(143) = ZZ(0,1);
   pars(144) = ZZ(1,0);
   pars(145) = ZZ(1,1);


   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

}


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHm_goldstone;
   MHm_goldstone(0) = MVWm;

   return remove_if_equal(MHm, MHm_goldstone);
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
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

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = M112 + 1.5*Lambda6*v1*v2 + 1.5*v1*v2*Conj(Lambda6) + 3
      *Lambda1*Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5*Lambda4*Sqr(v2) + 0.25*
      Lambda5*Sqr(v2) + 0.25*Conj(Lambda5)*Sqr(v2);
   mass_matrix_hh(0,1) = -0.5*M122 + Lambda3*v1*v2 + Lambda4*v1*v2 + 0.5*
      Lambda5*v1*v2 + 0.5*v1*v2*Conj(Lambda5) - 0.5*Conj(M122) + 0.75*Lambda6*
      Sqr(v1) + 0.75*Conj(Lambda6)*Sqr(v1) + 0.75*Lambda7*Sqr(v2) + 0.75*Conj(
      Lambda7)*Sqr(v2);
   mass_matrix_hh(1,1) = M222 + 1.5*Lambda7*v1*v2 + 1.5*v1*v2*Conj(Lambda7) +
      0.5*Lambda3*Sqr(v1) + 0.5*Lambda4*Sqr(v1) + 0.25*Lambda5*Sqr(v1) + 0.25*
      Conj(Lambda5)*Sqr(v1) + 3*Lambda2*Sqr(v2);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::hh, eigenvalue_error > precision *
      Abs(Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIIMSSMBC_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Conj(Lambda6) +
      Lambda1*Sqr(v1) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*
      Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5*Lambda4*Sqr(v2) - 0.25*Lambda5*Sqr(v2
      ) - 0.25*Conj(Lambda5)*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v1)*Sqr(Cos(ThetaW()))
      + 0.15*Sqr(g1)*Sqr(v1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = -0.5*M122 + 0.5*Lambda5*v1*v2 + 0.5*v1*v2*Conj(Lambda5
      ) - 0.5*Conj(M122) + 0.3872983346207417*g1*g2*v1*v2*Cos(ThetaW())*Sin(
      ThetaW()) + 0.25*Lambda6*Sqr(v1) + 0.25*Conj(Lambda6)*Sqr(v1) + 0.25*
      Lambda7*Sqr(v2) + 0.25*Conj(Lambda7)*Sqr(v2) + 0.25*v1*v2*Sqr(g2)*Sqr(Cos
      (ThetaW())) + 0.15*v1*v2*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Conj(Lambda7) +
      0.5*Lambda3*Sqr(v1) + 0.5*Lambda4*Sqr(v1) - 0.25*Lambda5*Sqr(v1) - 0.25*
      Conj(Lambda5)*Sqr(v1) + Lambda2*Sqr(v2) + 0.3872983346207417*g1*g2*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v2)*Sqr(Cos(ThetaW()))
      + 0.15*Sqr(g1)*Sqr(v2)*Sqr(Sin(ThetaW()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::Ah, eigenvalue_error > precision *
      Abs(MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIIMSSMBC_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hm;

   mass_matrix_Hm(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Conj(Lambda6) +
      Lambda1*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v1) + 0.5*Lambda3*Sqr(v2);
   mass_matrix_Hm(0,1) = 0.5*Lambda4*v1*v2 + 0.5*Lambda5*v1*v2 - Conj(M122) +
      0.25*v1*v2*Sqr(g2) + 0.5*Conj(Lambda6)*Sqr(v1) + 0.5*Conj(Lambda7)*Sqr(v2
      );
   mass_matrix_Hm(1,0) = -M122 + 0.5*Lambda4*v1*v2 + 0.5*v1*v2*Conj(Lambda5) +
      0.25*v1*v2*Sqr(g2) + 0.5*Lambda6*Sqr(v1) + 0.5*Lambda7*Sqr(v2);
   mass_matrix_Hm(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Conj(Lambda7) +
      0.5*Lambda3*Sqr(v1) + Lambda2*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v2);

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{
   const auto mass_matrix_Hm(get_mass_matrix_Hm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::Hm, eigenvalue_error > precision *
      Abs(MHm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHm.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIIMSSMBC_info::Hm);
   }

   MHm = AbsSqrt(MHm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v1*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v1*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v1*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v1*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v1*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v1*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v1*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v1*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v1*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::Fd, eigenvalue_error > precision *
      Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v2*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v2*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v2*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v2*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v2*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v2*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v2*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v2*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v2*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v1*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v1*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v1*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v1*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v1*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v1*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v1*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v1*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v1*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(THDMIIMSSMBC_info::Fe, eigenvalue_error > precision *
      Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(v1) + Sqr(v2)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(THDMIIMSSMBC_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v1) + 0.15*Sqr(g1)*Sqr(v2);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v1) -
      0.19364916731037085*g1*g2*Sqr(v2);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v2);

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
   
   double result = Re(M112*v1 - 0.5*M122*v2 - 0.5*v2*Conj(M122) + Lambda1*Cube(v1)
      + 0.25*Lambda7*Cube(v2) + 0.25*Conj(Lambda7)*Cube(v2) + 0.75*Lambda6*v2*Sqr(
      v1) + 0.75*v2*Conj(Lambda6)*Sqr(v1) + 0.5*Lambda3*v1*Sqr(v2) + 0.5*Lambda4*
      v1*Sqr(v2) + 0.25*Lambda5*v1*Sqr(v2) + 0.25*v1*Conj(Lambda5)*Sqr(v2));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(-0.5*M122*v1 + M222*v2 - 0.5*v1*Conj(M122) + 0.25*Lambda6*
      Cube(v1) + 0.25*Conj(Lambda6)*Cube(v1) + Lambda2*Cube(v2) + 0.5*Lambda3*v2*
      Sqr(v1) + 0.5*Lambda4*v2*Sqr(v1) + 0.25*Lambda5*v2*Sqr(v1) + 0.25*v2*Conj(
      Lambda5)*Sqr(v1) + 0.75*Lambda7*v1*Sqr(v2) + 0.75*v1*Conj(Lambda7)*Sqr(v2));

   return result;
}



double CLASSNAME::v() const
{

   return Sqrt(Sqr(v1) + Sqr(v2));
}

double CLASSNAME::Betax() const
{

   return ArcSin(Abs(ZP(0,1)));
}

double CLASSNAME::Alpha() const
{

   return ArcTan(ZH(1,1)/ZH(0,1));
}

double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}

double CLASSNAME::VEV() const
{

   return Sqrt(Sqr(v1) + Sqr(v2));
}



std::ostream& operator<<(std::ostream& ostr, const THDMIIMSSMBC_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
