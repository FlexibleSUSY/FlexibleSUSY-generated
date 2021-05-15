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
 * @file SplitMSSM_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the SplitMSSM model class in the decoupling scheme
 *
 * Contains the definition of the SplitMSSM model class methods
 * which solve EWSB and calculate masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#include "SplitMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "SplitMSSM_mass_eigenstates.hpp"
#include "SplitMSSM_info.hpp"
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

#define CLASSNAME SplitMSSM_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const SplitMSSM_input_parameters& input_)
   : SplitMSSM_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const SplitMSSM_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<SplitMSSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<SplitMSSM_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            SplitMSSM_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            SplitMSSM_info::g2, new_g2, get_scale());
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
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);

      MODEL->set_v(Re((1.4142135623730951*MZMSbar)/Sqrt(0.6*Sqr(g1) + Sqr(g2))));

   }

   // apply user-defined low-energy constraint for the Yukawa couplings
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Yu(((upQuarksDRbar/v).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Yd(((downQuarksDRbar/v).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Ye(((downLeptonsDRbar/v).transpose()).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const SplitMSSM_mass_eigenstates& model)
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
   MHp = OTHER(MHp);
   MFv = OTHER(MFv);
   MGlu = OTHER(MGlu);
   MAh = OTHER(MAh);
   Mhh = OTHER(Mhh);
   MFd = OTHER(MFd);
   Vd = OTHER(Vd);
   Ud = OTHER(Ud);
   MFu = OTHER(MFu);
   Vu = OTHER(Vu);
   Uu = OTHER(Uu);
   MFe = OTHER(MFe);
   Ve = OTHER(Ve);
   Ue = OTHER(Ue);
   MChi = OTHER(MChi);
   ZN = OTHER(ZN);
   MCha = OTHER(MCha);
   UM = OTHER(UM);
   UP = OTHER(UP);
   MVWp = OTHER(MVWp);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   ZZ = OTHER(ZZ);

#undef OTHER
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

   
   const double old_mu2 = mu2;

   mu2 = Re(Lambdax*Sqr(v));

   const bool is_finite = IsFinite(mu2);

   if (!is_finite) {
      mu2 = old_mu2;
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
           "SplitMSSM\n"
           "========================================\n";
   SplitMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
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
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
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
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_tree_level_mass_spectrum()
{
   const auto save_mu2_raii = make_raii_save(mu2);

   solve_ewsb_equations_tree_level();

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
   if (PHYSICAL(Mhh) < 0.) { problems.flag_pole_tachyon(SplitMSSM_info::hh); }
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
   problems.clear();
}

void CLASSNAME::clear()
{
   SplitMSSM_soft_parameters::clear();
   clear_tree_level_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& pars)
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

const SplitMSSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

SplitMSSM_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
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

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   Vd(0,0) = std::complex<double>(pars(26), pars(27));
   Vd(0,1) = std::complex<double>(pars(28), pars(29));
   Vd(0,2) = std::complex<double>(pars(30), pars(31));
   Vd(1,0) = std::complex<double>(pars(32), pars(33));
   Vd(1,1) = std::complex<double>(pars(34), pars(35));
   Vd(1,2) = std::complex<double>(pars(36), pars(37));
   Vd(2,0) = std::complex<double>(pars(38), pars(39));
   Vd(2,1) = std::complex<double>(pars(40), pars(41));
   Vd(2,2) = std::complex<double>(pars(42), pars(43));
   Ud(0,0) = std::complex<double>(pars(44), pars(45));
   Ud(0,1) = std::complex<double>(pars(46), pars(47));
   Ud(0,2) = std::complex<double>(pars(48), pars(49));
   Ud(1,0) = std::complex<double>(pars(50), pars(51));
   Ud(1,1) = std::complex<double>(pars(52), pars(53));
   Ud(1,2) = std::complex<double>(pars(54), pars(55));
   Ud(2,0) = std::complex<double>(pars(56), pars(57));
   Ud(2,1) = std::complex<double>(pars(58), pars(59));
   Ud(2,2) = std::complex<double>(pars(60), pars(61));
   Vu(0,0) = std::complex<double>(pars(62), pars(63));
   Vu(0,1) = std::complex<double>(pars(64), pars(65));
   Vu(0,2) = std::complex<double>(pars(66), pars(67));
   Vu(1,0) = std::complex<double>(pars(68), pars(69));
   Vu(1,1) = std::complex<double>(pars(70), pars(71));
   Vu(1,2) = std::complex<double>(pars(72), pars(73));
   Vu(2,0) = std::complex<double>(pars(74), pars(75));
   Vu(2,1) = std::complex<double>(pars(76), pars(77));
   Vu(2,2) = std::complex<double>(pars(78), pars(79));
   Uu(0,0) = std::complex<double>(pars(80), pars(81));
   Uu(0,1) = std::complex<double>(pars(82), pars(83));
   Uu(0,2) = std::complex<double>(pars(84), pars(85));
   Uu(1,0) = std::complex<double>(pars(86), pars(87));
   Uu(1,1) = std::complex<double>(pars(88), pars(89));
   Uu(1,2) = std::complex<double>(pars(90), pars(91));
   Uu(2,0) = std::complex<double>(pars(92), pars(93));
   Uu(2,1) = std::complex<double>(pars(94), pars(95));
   Uu(2,2) = std::complex<double>(pars(96), pars(97));
   Ve(0,0) = std::complex<double>(pars(98), pars(99));
   Ve(0,1) = std::complex<double>(pars(100), pars(101));
   Ve(0,2) = std::complex<double>(pars(102), pars(103));
   Ve(1,0) = std::complex<double>(pars(104), pars(105));
   Ve(1,1) = std::complex<double>(pars(106), pars(107));
   Ve(1,2) = std::complex<double>(pars(108), pars(109));
   Ve(2,0) = std::complex<double>(pars(110), pars(111));
   Ve(2,1) = std::complex<double>(pars(112), pars(113));
   Ve(2,2) = std::complex<double>(pars(114), pars(115));
   Ue(0,0) = std::complex<double>(pars(116), pars(117));
   Ue(0,1) = std::complex<double>(pars(118), pars(119));
   Ue(0,2) = std::complex<double>(pars(120), pars(121));
   Ue(1,0) = std::complex<double>(pars(122), pars(123));
   Ue(1,1) = std::complex<double>(pars(124), pars(125));
   Ue(1,2) = std::complex<double>(pars(126), pars(127));
   Ue(2,0) = std::complex<double>(pars(128), pars(129));
   Ue(2,1) = std::complex<double>(pars(130), pars(131));
   Ue(2,2) = std::complex<double>(pars(132), pars(133));
   ZN(0,0) = std::complex<double>(pars(134), pars(135));
   ZN(0,1) = std::complex<double>(pars(136), pars(137));
   ZN(0,2) = std::complex<double>(pars(138), pars(139));
   ZN(0,3) = std::complex<double>(pars(140), pars(141));
   ZN(1,0) = std::complex<double>(pars(142), pars(143));
   ZN(1,1) = std::complex<double>(pars(144), pars(145));
   ZN(1,2) = std::complex<double>(pars(146), pars(147));
   ZN(1,3) = std::complex<double>(pars(148), pars(149));
   ZN(2,0) = std::complex<double>(pars(150), pars(151));
   ZN(2,1) = std::complex<double>(pars(152), pars(153));
   ZN(2,2) = std::complex<double>(pars(154), pars(155));
   ZN(2,3) = std::complex<double>(pars(156), pars(157));
   ZN(3,0) = std::complex<double>(pars(158), pars(159));
   ZN(3,1) = std::complex<double>(pars(160), pars(161));
   ZN(3,2) = std::complex<double>(pars(162), pars(163));
   ZN(3,3) = std::complex<double>(pars(164), pars(165));
   UM(0,0) = std::complex<double>(pars(166), pars(167));
   UM(0,1) = std::complex<double>(pars(168), pars(169));
   UM(1,0) = std::complex<double>(pars(170), pars(171));
   UM(1,1) = std::complex<double>(pars(172), pars(173));
   UP(0,0) = std::complex<double>(pars(174), pars(175));
   UP(0,1) = std::complex<double>(pars(176), pars(177));
   UP(1,0) = std::complex<double>(pars(178), pars(179));
   UP(1,1) = std::complex<double>(pars(180), pars(181));
   ZZ(0,0) = pars(182);
   ZZ(0,1) = pars(183);
   ZZ(1,0) = pars(184);
   ZZ(1,1) = pars(185);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(186);

   pars(26) = Re(Vd(0,0));
   pars(27) = Im(Vd(0,0));
   pars(28) = Re(Vd(0,1));
   pars(29) = Im(Vd(0,1));
   pars(30) = Re(Vd(0,2));
   pars(31) = Im(Vd(0,2));
   pars(32) = Re(Vd(1,0));
   pars(33) = Im(Vd(1,0));
   pars(34) = Re(Vd(1,1));
   pars(35) = Im(Vd(1,1));
   pars(36) = Re(Vd(1,2));
   pars(37) = Im(Vd(1,2));
   pars(38) = Re(Vd(2,0));
   pars(39) = Im(Vd(2,0));
   pars(40) = Re(Vd(2,1));
   pars(41) = Im(Vd(2,1));
   pars(42) = Re(Vd(2,2));
   pars(43) = Im(Vd(2,2));
   pars(44) = Re(Ud(0,0));
   pars(45) = Im(Ud(0,0));
   pars(46) = Re(Ud(0,1));
   pars(47) = Im(Ud(0,1));
   pars(48) = Re(Ud(0,2));
   pars(49) = Im(Ud(0,2));
   pars(50) = Re(Ud(1,0));
   pars(51) = Im(Ud(1,0));
   pars(52) = Re(Ud(1,1));
   pars(53) = Im(Ud(1,1));
   pars(54) = Re(Ud(1,2));
   pars(55) = Im(Ud(1,2));
   pars(56) = Re(Ud(2,0));
   pars(57) = Im(Ud(2,0));
   pars(58) = Re(Ud(2,1));
   pars(59) = Im(Ud(2,1));
   pars(60) = Re(Ud(2,2));
   pars(61) = Im(Ud(2,2));
   pars(62) = Re(Vu(0,0));
   pars(63) = Im(Vu(0,0));
   pars(64) = Re(Vu(0,1));
   pars(65) = Im(Vu(0,1));
   pars(66) = Re(Vu(0,2));
   pars(67) = Im(Vu(0,2));
   pars(68) = Re(Vu(1,0));
   pars(69) = Im(Vu(1,0));
   pars(70) = Re(Vu(1,1));
   pars(71) = Im(Vu(1,1));
   pars(72) = Re(Vu(1,2));
   pars(73) = Im(Vu(1,2));
   pars(74) = Re(Vu(2,0));
   pars(75) = Im(Vu(2,0));
   pars(76) = Re(Vu(2,1));
   pars(77) = Im(Vu(2,1));
   pars(78) = Re(Vu(2,2));
   pars(79) = Im(Vu(2,2));
   pars(80) = Re(Uu(0,0));
   pars(81) = Im(Uu(0,0));
   pars(82) = Re(Uu(0,1));
   pars(83) = Im(Uu(0,1));
   pars(84) = Re(Uu(0,2));
   pars(85) = Im(Uu(0,2));
   pars(86) = Re(Uu(1,0));
   pars(87) = Im(Uu(1,0));
   pars(88) = Re(Uu(1,1));
   pars(89) = Im(Uu(1,1));
   pars(90) = Re(Uu(1,2));
   pars(91) = Im(Uu(1,2));
   pars(92) = Re(Uu(2,0));
   pars(93) = Im(Uu(2,0));
   pars(94) = Re(Uu(2,1));
   pars(95) = Im(Uu(2,1));
   pars(96) = Re(Uu(2,2));
   pars(97) = Im(Uu(2,2));
   pars(98) = Re(Ve(0,0));
   pars(99) = Im(Ve(0,0));
   pars(100) = Re(Ve(0,1));
   pars(101) = Im(Ve(0,1));
   pars(102) = Re(Ve(0,2));
   pars(103) = Im(Ve(0,2));
   pars(104) = Re(Ve(1,0));
   pars(105) = Im(Ve(1,0));
   pars(106) = Re(Ve(1,1));
   pars(107) = Im(Ve(1,1));
   pars(108) = Re(Ve(1,2));
   pars(109) = Im(Ve(1,2));
   pars(110) = Re(Ve(2,0));
   pars(111) = Im(Ve(2,0));
   pars(112) = Re(Ve(2,1));
   pars(113) = Im(Ve(2,1));
   pars(114) = Re(Ve(2,2));
   pars(115) = Im(Ve(2,2));
   pars(116) = Re(Ue(0,0));
   pars(117) = Im(Ue(0,0));
   pars(118) = Re(Ue(0,1));
   pars(119) = Im(Ue(0,1));
   pars(120) = Re(Ue(0,2));
   pars(121) = Im(Ue(0,2));
   pars(122) = Re(Ue(1,0));
   pars(123) = Im(Ue(1,0));
   pars(124) = Re(Ue(1,1));
   pars(125) = Im(Ue(1,1));
   pars(126) = Re(Ue(1,2));
   pars(127) = Im(Ue(1,2));
   pars(128) = Re(Ue(2,0));
   pars(129) = Im(Ue(2,0));
   pars(130) = Re(Ue(2,1));
   pars(131) = Im(Ue(2,1));
   pars(132) = Re(Ue(2,2));
   pars(133) = Im(Ue(2,2));
   pars(134) = Re(ZN(0,0));
   pars(135) = Im(ZN(0,0));
   pars(136) = Re(ZN(0,1));
   pars(137) = Im(ZN(0,1));
   pars(138) = Re(ZN(0,2));
   pars(139) = Im(ZN(0,2));
   pars(140) = Re(ZN(0,3));
   pars(141) = Im(ZN(0,3));
   pars(142) = Re(ZN(1,0));
   pars(143) = Im(ZN(1,0));
   pars(144) = Re(ZN(1,1));
   pars(145) = Im(ZN(1,1));
   pars(146) = Re(ZN(1,2));
   pars(147) = Im(ZN(1,2));
   pars(148) = Re(ZN(1,3));
   pars(149) = Im(ZN(1,3));
   pars(150) = Re(ZN(2,0));
   pars(151) = Im(ZN(2,0));
   pars(152) = Re(ZN(2,1));
   pars(153) = Im(ZN(2,1));
   pars(154) = Re(ZN(2,2));
   pars(155) = Im(ZN(2,2));
   pars(156) = Re(ZN(2,3));
   pars(157) = Im(ZN(2,3));
   pars(158) = Re(ZN(3,0));
   pars(159) = Im(ZN(3,0));
   pars(160) = Re(ZN(3,1));
   pars(161) = Im(ZN(3,1));
   pars(162) = Re(ZN(3,2));
   pars(163) = Im(ZN(3,2));
   pars(164) = Re(ZN(3,3));
   pars(165) = Im(ZN(3,3));
   pars(166) = Re(UM(0,0));
   pars(167) = Im(UM(0,0));
   pars(168) = Re(UM(0,1));
   pars(169) = Im(UM(0,1));
   pars(170) = Re(UM(1,0));
   pars(171) = Im(UM(1,0));
   pars(172) = Re(UM(1,1));
   pars(173) = Im(UM(1,1));
   pars(174) = Re(UP(0,0));
   pars(175) = Im(UP(0,0));
   pars(176) = Re(UP(0,1));
   pars(177) = Im(UP(0,1));
   pars(178) = Re(UP(1,0));
   pars(179) = Im(UP(1,0));
   pars(180) = Re(UP(1,1));
   pars(181) = Im(UP(1,1));
   pars(182) = ZZ(0,0);
   pars(183) = ZZ(0,1);
   pars(184) = ZZ(1,0);
   pars(185) = ZZ(1,1);


   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

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

double CLASSNAME::get_mass_matrix_Hp() const
{

   const double mass_matrix_Hp = Re(-mu2 + 0.5*(2*Lambdax + Sqr(g2))*Sqr(v));

   return mass_matrix_Hp;
}

void CLASSNAME::calculate_MHp()
{

   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = mass_matrix_Hp;

   if (MHp < 0.) {
      problems.flag_running_tachyon(SplitMSSM_info::Hp);
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

   const double mass_matrix_Ah = Re(-mu2 + Lambdax*Sqr(v) + 0.7745966692414834*
      g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(v) + 0.5*Sqr(g2)*Sqr(v)*Sqr(Cos(
      ThetaW())) + 0.3*Sqr(g1)*Sqr(v)*Sqr(Sin(ThetaW())));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{

   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_running_tachyon(SplitMSSM_info::Ah);
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
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_running_tachyon(SplitMSSM_info::hh);
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
   problems.flag_bad_mass(SplitMSSM_info::Fd, eigenvalue_error > precision *
      Abs(MFd(0)));
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
   problems.flag_bad_mass(SplitMSSM_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
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
   problems.flag_bad_mass(SplitMSSM_info::Fe, eigenvalue_error > precision *
      Abs(MFe(0)));
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
   problems.flag_bad_mass(SplitMSSM_info::Chi, eigenvalue_error > precision *
      Abs(MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

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
   problems.flag_bad_mass(SplitMSSM_info::Cha, eigenvalue_error > precision *
      Abs(MCha(0)));
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
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_running_tachyon(SplitMSSM_info::VWp);
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
   
   double result = Re(-1.4142135623730951*mu2*v + 1.4142135623730951*Cube(v)*
      Lambdax);

   return result;
}



double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}

double CLASSNAME::VEV() const
{

   return v;
}



std::ostream& operator<<(std::ostream& ostr, const SplitMSSM_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
