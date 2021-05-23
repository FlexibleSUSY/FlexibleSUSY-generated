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
 * @file SM_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the SM model class in the decoupling scheme
 *
 * Contains the definition of the SM model class methods
 * which solve EWSB and calculate masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.4 .
 */

#include "SM_mass_eigenstates_decoupling_scheme.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_info.hpp"
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

#define CLASSNAME SM_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const SM_input_parameters& input_)
   : SM_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const SM_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<SM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<SM_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(SM_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            SM_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(SM_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            SM_info::g2, new_g2, get_scale());
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

      MODEL->set_v(Re((2*MZMSbar)/Sqrt(0.6*Sqr(g1) + Sqr(g2))));

   }

   // apply user-defined low-energy constraint for the Yukawa couplings
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Yu((((1.4142135623730951*upQuarksDRbar)/v).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Yd((((1.4142135623730951*downQuarksDRbar)/v).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto v = MODELPARAMETER(v);
      MODEL->set_Ye((((1.4142135623730951*downLeptonsDRbar)/v).transpose()).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const SM_mass_eigenstates& model)
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
   MVWp = OTHER(MVWp);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const SM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

SM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const SM_physical& physical_)
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

   mu2 = Re(-0.5*Lambdax*Sqr(v));

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
           "SM\n"
           "========================================\n";
   SM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
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
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_Mhh();
   calculate_MAh();
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
   if (PHYSICAL(Mhh) < 0.) { problems.flag_pole_tachyon(SM_info::hh); }
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
   MVWp = 0.;
   MVP = 0.;
   MVZ = 0.;



}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   SM_soft_parameters::clear();
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
   MAh = pars(5);
   Mhh = pars(6);
   MFd(0) = pars(7);
   MFd(1) = pars(8);
   MFd(2) = pars(9);
   MFu(0) = pars(10);
   MFu(1) = pars(11);
   MFu(2) = pars(12);
   MFe(0) = pars(13);
   MFe(1) = pars(14);
   MFe(2) = pars(15);
   MVWp = pars(16);
   MVP = pars(17);
   MVZ = pars(18);

}

const SM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

SM_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(19);

   pars(0) = MVG;
   pars(1) = MHp;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MAh;
   pars(6) = Mhh;
   pars(7) = MFd(0);
   pars(8) = MFd(1);
   pars(9) = MFd(2);
   pars(10) = MFu(0);
   pars(11) = MFu(1);
   pars(12) = MFu(2);
   pars(13) = MFe(0);
   pars(14) = MFe(1);
   pars(15) = MFe(2);
   pars(16) = MVWp;
   pars(17) = MVP;
   pars(18) = MVZ;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   Vd(0,0) = std::complex<double>(pars(19), pars(20));
   Vd(0,1) = std::complex<double>(pars(21), pars(22));
   Vd(0,2) = std::complex<double>(pars(23), pars(24));
   Vd(1,0) = std::complex<double>(pars(25), pars(26));
   Vd(1,1) = std::complex<double>(pars(27), pars(28));
   Vd(1,2) = std::complex<double>(pars(29), pars(30));
   Vd(2,0) = std::complex<double>(pars(31), pars(32));
   Vd(2,1) = std::complex<double>(pars(33), pars(34));
   Vd(2,2) = std::complex<double>(pars(35), pars(36));
   Ud(0,0) = std::complex<double>(pars(37), pars(38));
   Ud(0,1) = std::complex<double>(pars(39), pars(40));
   Ud(0,2) = std::complex<double>(pars(41), pars(42));
   Ud(1,0) = std::complex<double>(pars(43), pars(44));
   Ud(1,1) = std::complex<double>(pars(45), pars(46));
   Ud(1,2) = std::complex<double>(pars(47), pars(48));
   Ud(2,0) = std::complex<double>(pars(49), pars(50));
   Ud(2,1) = std::complex<double>(pars(51), pars(52));
   Ud(2,2) = std::complex<double>(pars(53), pars(54));
   Vu(0,0) = std::complex<double>(pars(55), pars(56));
   Vu(0,1) = std::complex<double>(pars(57), pars(58));
   Vu(0,2) = std::complex<double>(pars(59), pars(60));
   Vu(1,0) = std::complex<double>(pars(61), pars(62));
   Vu(1,1) = std::complex<double>(pars(63), pars(64));
   Vu(1,2) = std::complex<double>(pars(65), pars(66));
   Vu(2,0) = std::complex<double>(pars(67), pars(68));
   Vu(2,1) = std::complex<double>(pars(69), pars(70));
   Vu(2,2) = std::complex<double>(pars(71), pars(72));
   Uu(0,0) = std::complex<double>(pars(73), pars(74));
   Uu(0,1) = std::complex<double>(pars(75), pars(76));
   Uu(0,2) = std::complex<double>(pars(77), pars(78));
   Uu(1,0) = std::complex<double>(pars(79), pars(80));
   Uu(1,1) = std::complex<double>(pars(81), pars(82));
   Uu(1,2) = std::complex<double>(pars(83), pars(84));
   Uu(2,0) = std::complex<double>(pars(85), pars(86));
   Uu(2,1) = std::complex<double>(pars(87), pars(88));
   Uu(2,2) = std::complex<double>(pars(89), pars(90));
   Ve(0,0) = std::complex<double>(pars(91), pars(92));
   Ve(0,1) = std::complex<double>(pars(93), pars(94));
   Ve(0,2) = std::complex<double>(pars(95), pars(96));
   Ve(1,0) = std::complex<double>(pars(97), pars(98));
   Ve(1,1) = std::complex<double>(pars(99), pars(100));
   Ve(1,2) = std::complex<double>(pars(101), pars(102));
   Ve(2,0) = std::complex<double>(pars(103), pars(104));
   Ve(2,1) = std::complex<double>(pars(105), pars(106));
   Ve(2,2) = std::complex<double>(pars(107), pars(108));
   Ue(0,0) = std::complex<double>(pars(109), pars(110));
   Ue(0,1) = std::complex<double>(pars(111), pars(112));
   Ue(0,2) = std::complex<double>(pars(113), pars(114));
   Ue(1,0) = std::complex<double>(pars(115), pars(116));
   Ue(1,1) = std::complex<double>(pars(117), pars(118));
   Ue(1,2) = std::complex<double>(pars(119), pars(120));
   Ue(2,0) = std::complex<double>(pars(121), pars(122));
   Ue(2,1) = std::complex<double>(pars(123), pars(124));
   Ue(2,2) = std::complex<double>(pars(125), pars(126));
   ZZ(0,0) = pars(127);
   ZZ(0,1) = pars(128);
   ZZ(1,0) = pars(129);
   ZZ(1,1) = pars(130);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(131);

   pars(19) = Re(Vd(0,0));
   pars(20) = Im(Vd(0,0));
   pars(21) = Re(Vd(0,1));
   pars(22) = Im(Vd(0,1));
   pars(23) = Re(Vd(0,2));
   pars(24) = Im(Vd(0,2));
   pars(25) = Re(Vd(1,0));
   pars(26) = Im(Vd(1,0));
   pars(27) = Re(Vd(1,1));
   pars(28) = Im(Vd(1,1));
   pars(29) = Re(Vd(1,2));
   pars(30) = Im(Vd(1,2));
   pars(31) = Re(Vd(2,0));
   pars(32) = Im(Vd(2,0));
   pars(33) = Re(Vd(2,1));
   pars(34) = Im(Vd(2,1));
   pars(35) = Re(Vd(2,2));
   pars(36) = Im(Vd(2,2));
   pars(37) = Re(Ud(0,0));
   pars(38) = Im(Ud(0,0));
   pars(39) = Re(Ud(0,1));
   pars(40) = Im(Ud(0,1));
   pars(41) = Re(Ud(0,2));
   pars(42) = Im(Ud(0,2));
   pars(43) = Re(Ud(1,0));
   pars(44) = Im(Ud(1,0));
   pars(45) = Re(Ud(1,1));
   pars(46) = Im(Ud(1,1));
   pars(47) = Re(Ud(1,2));
   pars(48) = Im(Ud(1,2));
   pars(49) = Re(Ud(2,0));
   pars(50) = Im(Ud(2,0));
   pars(51) = Re(Ud(2,1));
   pars(52) = Im(Ud(2,1));
   pars(53) = Re(Ud(2,2));
   pars(54) = Im(Ud(2,2));
   pars(55) = Re(Vu(0,0));
   pars(56) = Im(Vu(0,0));
   pars(57) = Re(Vu(0,1));
   pars(58) = Im(Vu(0,1));
   pars(59) = Re(Vu(0,2));
   pars(60) = Im(Vu(0,2));
   pars(61) = Re(Vu(1,0));
   pars(62) = Im(Vu(1,0));
   pars(63) = Re(Vu(1,1));
   pars(64) = Im(Vu(1,1));
   pars(65) = Re(Vu(1,2));
   pars(66) = Im(Vu(1,2));
   pars(67) = Re(Vu(2,0));
   pars(68) = Im(Vu(2,0));
   pars(69) = Re(Vu(2,1));
   pars(70) = Im(Vu(2,1));
   pars(71) = Re(Vu(2,2));
   pars(72) = Im(Vu(2,2));
   pars(73) = Re(Uu(0,0));
   pars(74) = Im(Uu(0,0));
   pars(75) = Re(Uu(0,1));
   pars(76) = Im(Uu(0,1));
   pars(77) = Re(Uu(0,2));
   pars(78) = Im(Uu(0,2));
   pars(79) = Re(Uu(1,0));
   pars(80) = Im(Uu(1,0));
   pars(81) = Re(Uu(1,1));
   pars(82) = Im(Uu(1,1));
   pars(83) = Re(Uu(1,2));
   pars(84) = Im(Uu(1,2));
   pars(85) = Re(Uu(2,0));
   pars(86) = Im(Uu(2,0));
   pars(87) = Re(Uu(2,1));
   pars(88) = Im(Uu(2,1));
   pars(89) = Re(Uu(2,2));
   pars(90) = Im(Uu(2,2));
   pars(91) = Re(Ve(0,0));
   pars(92) = Im(Ve(0,0));
   pars(93) = Re(Ve(0,1));
   pars(94) = Im(Ve(0,1));
   pars(95) = Re(Ve(0,2));
   pars(96) = Im(Ve(0,2));
   pars(97) = Re(Ve(1,0));
   pars(98) = Im(Ve(1,0));
   pars(99) = Re(Ve(1,1));
   pars(100) = Im(Ve(1,1));
   pars(101) = Re(Ve(1,2));
   pars(102) = Im(Ve(1,2));
   pars(103) = Re(Ve(2,0));
   pars(104) = Im(Ve(2,0));
   pars(105) = Re(Ve(2,1));
   pars(106) = Im(Ve(2,1));
   pars(107) = Re(Ve(2,2));
   pars(108) = Im(Ve(2,2));
   pars(109) = Re(Ue(0,0));
   pars(110) = Im(Ue(0,0));
   pars(111) = Re(Ue(0,1));
   pars(112) = Im(Ue(0,1));
   pars(113) = Re(Ue(0,2));
   pars(114) = Im(Ue(0,2));
   pars(115) = Re(Ue(1,0));
   pars(116) = Im(Ue(1,0));
   pars(117) = Re(Ue(1,1));
   pars(118) = Im(Ue(1,1));
   pars(119) = Re(Ue(1,2));
   pars(120) = Im(Ue(1,2));
   pars(121) = Re(Ue(2,0));
   pars(122) = Im(Ue(2,0));
   pars(123) = Re(Ue(2,1));
   pars(124) = Im(Ue(2,1));
   pars(125) = Re(Ue(2,2));
   pars(126) = Im(Ue(2,2));
   pars(127) = ZZ(0,0);
   pars(128) = ZZ(0,1);
   pars(129) = ZZ(1,0);
   pars(130) = ZZ(1,1);


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

   const double mass_matrix_Hp = Re(mu2 + 0.25*(2*Lambdax + Sqr(g2))*Sqr(v));

   return mass_matrix_Hp;
}

void CLASSNAME::calculate_MHp()
{

   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = mass_matrix_Hp;

   if (MHp < 0.) {
      problems.flag_running_tachyon(SM_info::Hp);
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

double CLASSNAME::get_mass_matrix_Ah() const
{

   const double mass_matrix_Ah = Re(0.25*(4*mu2 + 2*Lambdax*Sqr(v) +
      1.5491933384829668*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(v) + Sqr(g2)*Sqr
      (v)*Sqr(Cos(ThetaW())) + 0.6*Sqr(g1)*Sqr(v)*Sqr(Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{

   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_running_tachyon(SM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double CLASSNAME::get_mass_matrix_hh() const
{

   const double mass_matrix_hh = Re(mu2 + 1.5*Lambdax*Sqr(v));

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{

   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_running_tachyon(SM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(SM_info::Fd, eigenvalue_error > precision * Abs(MFd(0
      )));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(SM_info::Fu, eigenvalue_error > precision * Abs(MFu(0
      )));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(SM_info::Fe, eigenvalue_error > precision * Abs(MFe(0
      )));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double CLASSNAME::get_mass_matrix_VWp() const
{

   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void CLASSNAME::calculate_MVWp()
{

   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_running_tachyon(SM_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

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
   
   double result = Re(mu2*v + 0.5*Cube(v)*Lambdax);

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



std::ostream& operator<<(std::ostream& ostr, const SM_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
