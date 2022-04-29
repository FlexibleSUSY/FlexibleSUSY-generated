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
 * @file lowMSSM_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the lowMSSM model class in the decoupling scheme
 *
 * Contains the definition of the lowMSSM model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include "lowMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "lowMSSM_mass_eigenstates.hpp"
#include "lowMSSM_info.hpp"
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

#define CLASSNAME lowMSSM_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const lowMSSM_input_parameters& input_)
   : lowMSSM_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const lowMSSM_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<lowMSSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<lowMSSM_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(lowMSSM_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            lowMSSM_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(lowMSSM_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            lowMSSM_info::g2, new_g2, get_scale());
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
      MODEL->set_Yu((((1.4142135623730951*upQuarksDRbar)/vu).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Yd((((1.4142135623730951*downQuarksDRbar)/vd).transpose()).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Ye((((1.4142135623730951*downLeptonsDRbar)/vd).transpose()).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const lowMSSM_mass_eigenstates& model)
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
   MHpm = OTHER(MHpm);
   ZP = OTHER(ZP);
   MChi = OTHER(MChi);
   ZN = OTHER(ZN);
   MCha = OTHER(MCha);
   UM = OTHER(UM);
   UP = OTHER(UP);
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

const lowMSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

lowMSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const lowMSSM_physical& physical_)
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

   
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) - 3*Cube(vd)*
      Sqr(g1) - 5*Cube(vd)*Sqr(g2) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu)))
      /vd);
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) - 3*Cube(vu)*
      Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd)))
      /vu);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2);

   if (!is_finite) {
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
           "lowMSSM\n"
           "========================================\n";
   lowMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
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
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
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

   solve_ewsb_equations_tree_level();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
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
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Se); }
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(lowMSSM_info::Hpm); }
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
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
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

   PhaseGlu = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   lowMSSM_soft_parameters::clear();
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
   MSd(0) = pars(5);
   MSd(1) = pars(6);
   MSd(2) = pars(7);
   MSd(3) = pars(8);
   MSd(4) = pars(9);
   MSd(5) = pars(10);
   MSv(0) = pars(11);
   MSv(1) = pars(12);
   MSv(2) = pars(13);
   MSu(0) = pars(14);
   MSu(1) = pars(15);
   MSu(2) = pars(16);
   MSu(3) = pars(17);
   MSu(4) = pars(18);
   MSu(5) = pars(19);
   MSe(0) = pars(20);
   MSe(1) = pars(21);
   MSe(2) = pars(22);
   MSe(3) = pars(23);
   MSe(4) = pars(24);
   MSe(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   MAh(0) = pars(28);
   MAh(1) = pars(29);
   MHpm(0) = pars(30);
   MHpm(1) = pars(31);
   MChi(0) = pars(32);
   MChi(1) = pars(33);
   MChi(2) = pars(34);
   MChi(3) = pars(35);
   MCha(0) = pars(36);
   MCha(1) = pars(37);
   MFe(0) = pars(38);
   MFe(1) = pars(39);
   MFe(2) = pars(40);
   MFd(0) = pars(41);
   MFd(1) = pars(42);
   MFd(2) = pars(43);
   MFu(0) = pars(44);
   MFu(1) = pars(45);
   MFu(2) = pars(46);
   MVWm = pars(47);
   MVP = pars(48);
   MVZ = pars(49);

}

const lowMSSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

lowMSSM_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(50);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSd(0);
   pars(6) = MSd(1);
   pars(7) = MSd(2);
   pars(8) = MSd(3);
   pars(9) = MSd(4);
   pars(10) = MSd(5);
   pars(11) = MSv(0);
   pars(12) = MSv(1);
   pars(13) = MSv(2);
   pars(14) = MSu(0);
   pars(15) = MSu(1);
   pars(16) = MSu(2);
   pars(17) = MSu(3);
   pars(18) = MSu(4);
   pars(19) = MSu(5);
   pars(20) = MSe(0);
   pars(21) = MSe(1);
   pars(22) = MSe(2);
   pars(23) = MSe(3);
   pars(24) = MSe(4);
   pars(25) = MSe(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = MAh(0);
   pars(29) = MAh(1);
   pars(30) = MHpm(0);
   pars(31) = MHpm(1);
   pars(32) = MChi(0);
   pars(33) = MChi(1);
   pars(34) = MChi(2);
   pars(35) = MChi(3);
   pars(36) = MCha(0);
   pars(37) = MCha(1);
   pars(38) = MFe(0);
   pars(39) = MFe(1);
   pars(40) = MFe(2);
   pars(41) = MFd(0);
   pars(42) = MFd(1);
   pars(43) = MFd(2);
   pars(44) = MFu(0);
   pars(45) = MFu(1);
   pars(46) = MFu(2);
   pars(47) = MVWm;
   pars(48) = MVP;
   pars(49) = MVZ;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZD(0,0) = pars(50);
   ZD(0,1) = pars(51);
   ZD(0,2) = pars(52);
   ZD(0,3) = pars(53);
   ZD(0,4) = pars(54);
   ZD(0,5) = pars(55);
   ZD(1,0) = pars(56);
   ZD(1,1) = pars(57);
   ZD(1,2) = pars(58);
   ZD(1,3) = pars(59);
   ZD(1,4) = pars(60);
   ZD(1,5) = pars(61);
   ZD(2,0) = pars(62);
   ZD(2,1) = pars(63);
   ZD(2,2) = pars(64);
   ZD(2,3) = pars(65);
   ZD(2,4) = pars(66);
   ZD(2,5) = pars(67);
   ZD(3,0) = pars(68);
   ZD(3,1) = pars(69);
   ZD(3,2) = pars(70);
   ZD(3,3) = pars(71);
   ZD(3,4) = pars(72);
   ZD(3,5) = pars(73);
   ZD(4,0) = pars(74);
   ZD(4,1) = pars(75);
   ZD(4,2) = pars(76);
   ZD(4,3) = pars(77);
   ZD(4,4) = pars(78);
   ZD(4,5) = pars(79);
   ZD(5,0) = pars(80);
   ZD(5,1) = pars(81);
   ZD(5,2) = pars(82);
   ZD(5,3) = pars(83);
   ZD(5,4) = pars(84);
   ZD(5,5) = pars(85);
   ZV(0,0) = pars(86);
   ZV(0,1) = pars(87);
   ZV(0,2) = pars(88);
   ZV(1,0) = pars(89);
   ZV(1,1) = pars(90);
   ZV(1,2) = pars(91);
   ZV(2,0) = pars(92);
   ZV(2,1) = pars(93);
   ZV(2,2) = pars(94);
   ZU(0,0) = pars(95);
   ZU(0,1) = pars(96);
   ZU(0,2) = pars(97);
   ZU(0,3) = pars(98);
   ZU(0,4) = pars(99);
   ZU(0,5) = pars(100);
   ZU(1,0) = pars(101);
   ZU(1,1) = pars(102);
   ZU(1,2) = pars(103);
   ZU(1,3) = pars(104);
   ZU(1,4) = pars(105);
   ZU(1,5) = pars(106);
   ZU(2,0) = pars(107);
   ZU(2,1) = pars(108);
   ZU(2,2) = pars(109);
   ZU(2,3) = pars(110);
   ZU(2,4) = pars(111);
   ZU(2,5) = pars(112);
   ZU(3,0) = pars(113);
   ZU(3,1) = pars(114);
   ZU(3,2) = pars(115);
   ZU(3,3) = pars(116);
   ZU(3,4) = pars(117);
   ZU(3,5) = pars(118);
   ZU(4,0) = pars(119);
   ZU(4,1) = pars(120);
   ZU(4,2) = pars(121);
   ZU(4,3) = pars(122);
   ZU(4,4) = pars(123);
   ZU(4,5) = pars(124);
   ZU(5,0) = pars(125);
   ZU(5,1) = pars(126);
   ZU(5,2) = pars(127);
   ZU(5,3) = pars(128);
   ZU(5,4) = pars(129);
   ZU(5,5) = pars(130);
   ZE(0,0) = pars(131);
   ZE(0,1) = pars(132);
   ZE(0,2) = pars(133);
   ZE(0,3) = pars(134);
   ZE(0,4) = pars(135);
   ZE(0,5) = pars(136);
   ZE(1,0) = pars(137);
   ZE(1,1) = pars(138);
   ZE(1,2) = pars(139);
   ZE(1,3) = pars(140);
   ZE(1,4) = pars(141);
   ZE(1,5) = pars(142);
   ZE(2,0) = pars(143);
   ZE(2,1) = pars(144);
   ZE(2,2) = pars(145);
   ZE(2,3) = pars(146);
   ZE(2,4) = pars(147);
   ZE(2,5) = pars(148);
   ZE(3,0) = pars(149);
   ZE(3,1) = pars(150);
   ZE(3,2) = pars(151);
   ZE(3,3) = pars(152);
   ZE(3,4) = pars(153);
   ZE(3,5) = pars(154);
   ZE(4,0) = pars(155);
   ZE(4,1) = pars(156);
   ZE(4,2) = pars(157);
   ZE(4,3) = pars(158);
   ZE(4,4) = pars(159);
   ZE(4,5) = pars(160);
   ZE(5,0) = pars(161);
   ZE(5,1) = pars(162);
   ZE(5,2) = pars(163);
   ZE(5,3) = pars(164);
   ZE(5,4) = pars(165);
   ZE(5,5) = pars(166);
   ZH(0,0) = pars(167);
   ZH(0,1) = pars(168);
   ZH(1,0) = pars(169);
   ZH(1,1) = pars(170);
   ZA(0,0) = pars(171);
   ZA(0,1) = pars(172);
   ZA(1,0) = pars(173);
   ZA(1,1) = pars(174);
   ZP(0,0) = pars(175);
   ZP(0,1) = pars(176);
   ZP(1,0) = pars(177);
   ZP(1,1) = pars(178);
   ZN(0,0) = std::complex<double>(pars(179), pars(180));
   ZN(0,1) = std::complex<double>(pars(181), pars(182));
   ZN(0,2) = std::complex<double>(pars(183), pars(184));
   ZN(0,3) = std::complex<double>(pars(185), pars(186));
   ZN(1,0) = std::complex<double>(pars(187), pars(188));
   ZN(1,1) = std::complex<double>(pars(189), pars(190));
   ZN(1,2) = std::complex<double>(pars(191), pars(192));
   ZN(1,3) = std::complex<double>(pars(193), pars(194));
   ZN(2,0) = std::complex<double>(pars(195), pars(196));
   ZN(2,1) = std::complex<double>(pars(197), pars(198));
   ZN(2,2) = std::complex<double>(pars(199), pars(200));
   ZN(2,3) = std::complex<double>(pars(201), pars(202));
   ZN(3,0) = std::complex<double>(pars(203), pars(204));
   ZN(3,1) = std::complex<double>(pars(205), pars(206));
   ZN(3,2) = std::complex<double>(pars(207), pars(208));
   ZN(3,3) = std::complex<double>(pars(209), pars(210));
   UM(0,0) = std::complex<double>(pars(211), pars(212));
   UM(0,1) = std::complex<double>(pars(213), pars(214));
   UM(1,0) = std::complex<double>(pars(215), pars(216));
   UM(1,1) = std::complex<double>(pars(217), pars(218));
   UP(0,0) = std::complex<double>(pars(219), pars(220));
   UP(0,1) = std::complex<double>(pars(221), pars(222));
   UP(1,0) = std::complex<double>(pars(223), pars(224));
   UP(1,1) = std::complex<double>(pars(225), pars(226));
   ZEL(0,0) = std::complex<double>(pars(227), pars(228));
   ZEL(0,1) = std::complex<double>(pars(229), pars(230));
   ZEL(0,2) = std::complex<double>(pars(231), pars(232));
   ZEL(1,0) = std::complex<double>(pars(233), pars(234));
   ZEL(1,1) = std::complex<double>(pars(235), pars(236));
   ZEL(1,2) = std::complex<double>(pars(237), pars(238));
   ZEL(2,0) = std::complex<double>(pars(239), pars(240));
   ZEL(2,1) = std::complex<double>(pars(241), pars(242));
   ZEL(2,2) = std::complex<double>(pars(243), pars(244));
   ZER(0,0) = std::complex<double>(pars(245), pars(246));
   ZER(0,1) = std::complex<double>(pars(247), pars(248));
   ZER(0,2) = std::complex<double>(pars(249), pars(250));
   ZER(1,0) = std::complex<double>(pars(251), pars(252));
   ZER(1,1) = std::complex<double>(pars(253), pars(254));
   ZER(1,2) = std::complex<double>(pars(255), pars(256));
   ZER(2,0) = std::complex<double>(pars(257), pars(258));
   ZER(2,1) = std::complex<double>(pars(259), pars(260));
   ZER(2,2) = std::complex<double>(pars(261), pars(262));
   ZDL(0,0) = std::complex<double>(pars(263), pars(264));
   ZDL(0,1) = std::complex<double>(pars(265), pars(266));
   ZDL(0,2) = std::complex<double>(pars(267), pars(268));
   ZDL(1,0) = std::complex<double>(pars(269), pars(270));
   ZDL(1,1) = std::complex<double>(pars(271), pars(272));
   ZDL(1,2) = std::complex<double>(pars(273), pars(274));
   ZDL(2,0) = std::complex<double>(pars(275), pars(276));
   ZDL(2,1) = std::complex<double>(pars(277), pars(278));
   ZDL(2,2) = std::complex<double>(pars(279), pars(280));
   ZDR(0,0) = std::complex<double>(pars(281), pars(282));
   ZDR(0,1) = std::complex<double>(pars(283), pars(284));
   ZDR(0,2) = std::complex<double>(pars(285), pars(286));
   ZDR(1,0) = std::complex<double>(pars(287), pars(288));
   ZDR(1,1) = std::complex<double>(pars(289), pars(290));
   ZDR(1,2) = std::complex<double>(pars(291), pars(292));
   ZDR(2,0) = std::complex<double>(pars(293), pars(294));
   ZDR(2,1) = std::complex<double>(pars(295), pars(296));
   ZDR(2,2) = std::complex<double>(pars(297), pars(298));
   ZUL(0,0) = std::complex<double>(pars(299), pars(300));
   ZUL(0,1) = std::complex<double>(pars(301), pars(302));
   ZUL(0,2) = std::complex<double>(pars(303), pars(304));
   ZUL(1,0) = std::complex<double>(pars(305), pars(306));
   ZUL(1,1) = std::complex<double>(pars(307), pars(308));
   ZUL(1,2) = std::complex<double>(pars(309), pars(310));
   ZUL(2,0) = std::complex<double>(pars(311), pars(312));
   ZUL(2,1) = std::complex<double>(pars(313), pars(314));
   ZUL(2,2) = std::complex<double>(pars(315), pars(316));
   ZUR(0,0) = std::complex<double>(pars(317), pars(318));
   ZUR(0,1) = std::complex<double>(pars(319), pars(320));
   ZUR(0,2) = std::complex<double>(pars(321), pars(322));
   ZUR(1,0) = std::complex<double>(pars(323), pars(324));
   ZUR(1,1) = std::complex<double>(pars(325), pars(326));
   ZUR(1,2) = std::complex<double>(pars(327), pars(328));
   ZUR(2,0) = std::complex<double>(pars(329), pars(330));
   ZUR(2,1) = std::complex<double>(pars(331), pars(332));
   ZUR(2,2) = std::complex<double>(pars(333), pars(334));
   ZZ(0,0) = pars(335);
   ZZ(0,1) = pars(336);
   ZZ(1,0) = pars(337);
   ZZ(1,1) = pars(338);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(339);

   pars(50) = ZD(0,0);
   pars(51) = ZD(0,1);
   pars(52) = ZD(0,2);
   pars(53) = ZD(0,3);
   pars(54) = ZD(0,4);
   pars(55) = ZD(0,5);
   pars(56) = ZD(1,0);
   pars(57) = ZD(1,1);
   pars(58) = ZD(1,2);
   pars(59) = ZD(1,3);
   pars(60) = ZD(1,4);
   pars(61) = ZD(1,5);
   pars(62) = ZD(2,0);
   pars(63) = ZD(2,1);
   pars(64) = ZD(2,2);
   pars(65) = ZD(2,3);
   pars(66) = ZD(2,4);
   pars(67) = ZD(2,5);
   pars(68) = ZD(3,0);
   pars(69) = ZD(3,1);
   pars(70) = ZD(3,2);
   pars(71) = ZD(3,3);
   pars(72) = ZD(3,4);
   pars(73) = ZD(3,5);
   pars(74) = ZD(4,0);
   pars(75) = ZD(4,1);
   pars(76) = ZD(4,2);
   pars(77) = ZD(4,3);
   pars(78) = ZD(4,4);
   pars(79) = ZD(4,5);
   pars(80) = ZD(5,0);
   pars(81) = ZD(5,1);
   pars(82) = ZD(5,2);
   pars(83) = ZD(5,3);
   pars(84) = ZD(5,4);
   pars(85) = ZD(5,5);
   pars(86) = ZV(0,0);
   pars(87) = ZV(0,1);
   pars(88) = ZV(0,2);
   pars(89) = ZV(1,0);
   pars(90) = ZV(1,1);
   pars(91) = ZV(1,2);
   pars(92) = ZV(2,0);
   pars(93) = ZV(2,1);
   pars(94) = ZV(2,2);
   pars(95) = ZU(0,0);
   pars(96) = ZU(0,1);
   pars(97) = ZU(0,2);
   pars(98) = ZU(0,3);
   pars(99) = ZU(0,4);
   pars(100) = ZU(0,5);
   pars(101) = ZU(1,0);
   pars(102) = ZU(1,1);
   pars(103) = ZU(1,2);
   pars(104) = ZU(1,3);
   pars(105) = ZU(1,4);
   pars(106) = ZU(1,5);
   pars(107) = ZU(2,0);
   pars(108) = ZU(2,1);
   pars(109) = ZU(2,2);
   pars(110) = ZU(2,3);
   pars(111) = ZU(2,4);
   pars(112) = ZU(2,5);
   pars(113) = ZU(3,0);
   pars(114) = ZU(3,1);
   pars(115) = ZU(3,2);
   pars(116) = ZU(3,3);
   pars(117) = ZU(3,4);
   pars(118) = ZU(3,5);
   pars(119) = ZU(4,0);
   pars(120) = ZU(4,1);
   pars(121) = ZU(4,2);
   pars(122) = ZU(4,3);
   pars(123) = ZU(4,4);
   pars(124) = ZU(4,5);
   pars(125) = ZU(5,0);
   pars(126) = ZU(5,1);
   pars(127) = ZU(5,2);
   pars(128) = ZU(5,3);
   pars(129) = ZU(5,4);
   pars(130) = ZU(5,5);
   pars(131) = ZE(0,0);
   pars(132) = ZE(0,1);
   pars(133) = ZE(0,2);
   pars(134) = ZE(0,3);
   pars(135) = ZE(0,4);
   pars(136) = ZE(0,5);
   pars(137) = ZE(1,0);
   pars(138) = ZE(1,1);
   pars(139) = ZE(1,2);
   pars(140) = ZE(1,3);
   pars(141) = ZE(1,4);
   pars(142) = ZE(1,5);
   pars(143) = ZE(2,0);
   pars(144) = ZE(2,1);
   pars(145) = ZE(2,2);
   pars(146) = ZE(2,3);
   pars(147) = ZE(2,4);
   pars(148) = ZE(2,5);
   pars(149) = ZE(3,0);
   pars(150) = ZE(3,1);
   pars(151) = ZE(3,2);
   pars(152) = ZE(3,3);
   pars(153) = ZE(3,4);
   pars(154) = ZE(3,5);
   pars(155) = ZE(4,0);
   pars(156) = ZE(4,1);
   pars(157) = ZE(4,2);
   pars(158) = ZE(4,3);
   pars(159) = ZE(4,4);
   pars(160) = ZE(4,5);
   pars(161) = ZE(5,0);
   pars(162) = ZE(5,1);
   pars(163) = ZE(5,2);
   pars(164) = ZE(5,3);
   pars(165) = ZE(5,4);
   pars(166) = ZE(5,5);
   pars(167) = ZH(0,0);
   pars(168) = ZH(0,1);
   pars(169) = ZH(1,0);
   pars(170) = ZH(1,1);
   pars(171) = ZA(0,0);
   pars(172) = ZA(0,1);
   pars(173) = ZA(1,0);
   pars(174) = ZA(1,1);
   pars(175) = ZP(0,0);
   pars(176) = ZP(0,1);
   pars(177) = ZP(1,0);
   pars(178) = ZP(1,1);
   pars(179) = Re(ZN(0,0));
   pars(180) = Im(ZN(0,0));
   pars(181) = Re(ZN(0,1));
   pars(182) = Im(ZN(0,1));
   pars(183) = Re(ZN(0,2));
   pars(184) = Im(ZN(0,2));
   pars(185) = Re(ZN(0,3));
   pars(186) = Im(ZN(0,3));
   pars(187) = Re(ZN(1,0));
   pars(188) = Im(ZN(1,0));
   pars(189) = Re(ZN(1,1));
   pars(190) = Im(ZN(1,1));
   pars(191) = Re(ZN(1,2));
   pars(192) = Im(ZN(1,2));
   pars(193) = Re(ZN(1,3));
   pars(194) = Im(ZN(1,3));
   pars(195) = Re(ZN(2,0));
   pars(196) = Im(ZN(2,0));
   pars(197) = Re(ZN(2,1));
   pars(198) = Im(ZN(2,1));
   pars(199) = Re(ZN(2,2));
   pars(200) = Im(ZN(2,2));
   pars(201) = Re(ZN(2,3));
   pars(202) = Im(ZN(2,3));
   pars(203) = Re(ZN(3,0));
   pars(204) = Im(ZN(3,0));
   pars(205) = Re(ZN(3,1));
   pars(206) = Im(ZN(3,1));
   pars(207) = Re(ZN(3,2));
   pars(208) = Im(ZN(3,2));
   pars(209) = Re(ZN(3,3));
   pars(210) = Im(ZN(3,3));
   pars(211) = Re(UM(0,0));
   pars(212) = Im(UM(0,0));
   pars(213) = Re(UM(0,1));
   pars(214) = Im(UM(0,1));
   pars(215) = Re(UM(1,0));
   pars(216) = Im(UM(1,0));
   pars(217) = Re(UM(1,1));
   pars(218) = Im(UM(1,1));
   pars(219) = Re(UP(0,0));
   pars(220) = Im(UP(0,0));
   pars(221) = Re(UP(0,1));
   pars(222) = Im(UP(0,1));
   pars(223) = Re(UP(1,0));
   pars(224) = Im(UP(1,0));
   pars(225) = Re(UP(1,1));
   pars(226) = Im(UP(1,1));
   pars(227) = Re(ZEL(0,0));
   pars(228) = Im(ZEL(0,0));
   pars(229) = Re(ZEL(0,1));
   pars(230) = Im(ZEL(0,1));
   pars(231) = Re(ZEL(0,2));
   pars(232) = Im(ZEL(0,2));
   pars(233) = Re(ZEL(1,0));
   pars(234) = Im(ZEL(1,0));
   pars(235) = Re(ZEL(1,1));
   pars(236) = Im(ZEL(1,1));
   pars(237) = Re(ZEL(1,2));
   pars(238) = Im(ZEL(1,2));
   pars(239) = Re(ZEL(2,0));
   pars(240) = Im(ZEL(2,0));
   pars(241) = Re(ZEL(2,1));
   pars(242) = Im(ZEL(2,1));
   pars(243) = Re(ZEL(2,2));
   pars(244) = Im(ZEL(2,2));
   pars(245) = Re(ZER(0,0));
   pars(246) = Im(ZER(0,0));
   pars(247) = Re(ZER(0,1));
   pars(248) = Im(ZER(0,1));
   pars(249) = Re(ZER(0,2));
   pars(250) = Im(ZER(0,2));
   pars(251) = Re(ZER(1,0));
   pars(252) = Im(ZER(1,0));
   pars(253) = Re(ZER(1,1));
   pars(254) = Im(ZER(1,1));
   pars(255) = Re(ZER(1,2));
   pars(256) = Im(ZER(1,2));
   pars(257) = Re(ZER(2,0));
   pars(258) = Im(ZER(2,0));
   pars(259) = Re(ZER(2,1));
   pars(260) = Im(ZER(2,1));
   pars(261) = Re(ZER(2,2));
   pars(262) = Im(ZER(2,2));
   pars(263) = Re(ZDL(0,0));
   pars(264) = Im(ZDL(0,0));
   pars(265) = Re(ZDL(0,1));
   pars(266) = Im(ZDL(0,1));
   pars(267) = Re(ZDL(0,2));
   pars(268) = Im(ZDL(0,2));
   pars(269) = Re(ZDL(1,0));
   pars(270) = Im(ZDL(1,0));
   pars(271) = Re(ZDL(1,1));
   pars(272) = Im(ZDL(1,1));
   pars(273) = Re(ZDL(1,2));
   pars(274) = Im(ZDL(1,2));
   pars(275) = Re(ZDL(2,0));
   pars(276) = Im(ZDL(2,0));
   pars(277) = Re(ZDL(2,1));
   pars(278) = Im(ZDL(2,1));
   pars(279) = Re(ZDL(2,2));
   pars(280) = Im(ZDL(2,2));
   pars(281) = Re(ZDR(0,0));
   pars(282) = Im(ZDR(0,0));
   pars(283) = Re(ZDR(0,1));
   pars(284) = Im(ZDR(0,1));
   pars(285) = Re(ZDR(0,2));
   pars(286) = Im(ZDR(0,2));
   pars(287) = Re(ZDR(1,0));
   pars(288) = Im(ZDR(1,0));
   pars(289) = Re(ZDR(1,1));
   pars(290) = Im(ZDR(1,1));
   pars(291) = Re(ZDR(1,2));
   pars(292) = Im(ZDR(1,2));
   pars(293) = Re(ZDR(2,0));
   pars(294) = Im(ZDR(2,0));
   pars(295) = Re(ZDR(2,1));
   pars(296) = Im(ZDR(2,1));
   pars(297) = Re(ZDR(2,2));
   pars(298) = Im(ZDR(2,2));
   pars(299) = Re(ZUL(0,0));
   pars(300) = Im(ZUL(0,0));
   pars(301) = Re(ZUL(0,1));
   pars(302) = Im(ZUL(0,1));
   pars(303) = Re(ZUL(0,2));
   pars(304) = Im(ZUL(0,2));
   pars(305) = Re(ZUL(1,0));
   pars(306) = Im(ZUL(1,0));
   pars(307) = Re(ZUL(1,1));
   pars(308) = Im(ZUL(1,1));
   pars(309) = Re(ZUL(1,2));
   pars(310) = Im(ZUL(1,2));
   pars(311) = Re(ZUL(2,0));
   pars(312) = Im(ZUL(2,0));
   pars(313) = Re(ZUL(2,1));
   pars(314) = Im(ZUL(2,1));
   pars(315) = Re(ZUL(2,2));
   pars(316) = Im(ZUL(2,2));
   pars(317) = Re(ZUR(0,0));
   pars(318) = Im(ZUR(0,0));
   pars(319) = Re(ZUR(0,1));
   pars(320) = Im(ZUR(0,1));
   pars(321) = Re(ZUR(0,2));
   pars(322) = Im(ZUR(0,2));
   pars(323) = Re(ZUR(1,0));
   pars(324) = Im(ZUR(1,0));
   pars(325) = Re(ZUR(1,1));
   pars(326) = Im(ZUR(1,1));
   pars(327) = Re(ZUR(1,2));
   pars(328) = Im(ZUR(1,2));
   pars(329) = Re(ZUR(2,0));
   pars(330) = Im(ZUR(2,0));
   pars(331) = Re(ZUR(2,1));
   pars(332) = Im(ZUR(2,1));
   pars(333) = Re(ZUR(2,2));
   pars(334) = Im(ZUR(2,2));
   pars(335) = ZZ(0,0);
   pars(336) = ZZ(0,1);
   pars(337) = ZZ(1,0);
   pars(338) = ZZ(1,1);


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
   Eigen::Array<double,1,1> MHpm_goldstone;
   MHpm_goldstone(0) = MVWm;

   return remove_if_equal(MHpm, MHpm_goldstone);
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

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) + Conj(
      Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) + Conj(
      Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) -
      0.7071067811865475*vu*Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) -
      0.7071067811865475*vu*Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) -
      0.7071067811865475*vu*Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) + Conj(
      Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) -
      0.7071067811865475*vu*Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) -
      0.7071067811865475*vu*Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) -
      0.7071067811865475*vu*Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) -
      0.7071067811865475*vu*Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) -
      0.7071067811865475*vu*Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) -
      0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)) +
      AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) + Conj(
      Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) + Conj(
      Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) + Conj(
      Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)) +
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
   problems.flag_bad_mass(lowMSSM_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) + Conj(
      Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) + Conj(
      Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) -
      0.7071067811865475*vd*Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) -
      0.7071067811865475*vd*Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) -
      0.7071067811865475*vd*Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) -
      0.7071067811865475*vd*Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) -
      0.7071067811865475*vd*Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) -
      0.7071067811865475*vd*Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) -
      0.7071067811865475*vd*Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) -
      0.7071067811865475*vd*Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) -
      0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,0))
      + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) + Conj(
      Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) + Conj(
      Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(1,0))
      + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) + Conj(
      Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(2,0))
      + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) + Conj(
      Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) + Conj(
      Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) -
      0.7071067811865475*vu*Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) -
      0.7071067811865475*vu*Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) -
      0.7071067811865475*vu*Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) + Conj(
      Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) -
      0.7071067811865475*vu*Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) -
      0.7071067811865475*vu*Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) -
      0.7071067811865475*vu*Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) -
      0.7071067811865475*vu*Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) -
      0.7071067811865475*vu*Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) -
      0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)) +
      AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) + Conj(
      Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) + Conj(
      Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) + Conj(
      Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)) +
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
   problems.flag_bad_mass(lowMSSM_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + AbsSqr(Mu) + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(
      g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*Conj(BMu) - 0.15*vd*vu*Sqr(g1) - 0.25*
      vd*vu*Sqr(g2);
   mass_matrix_hh(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(
      g2)*Sqr(vd) + 0.225*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + AbsSqr(Mu) + 0.3872983346207417*g1*g2*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*
      Sqr(vd)*Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) - 0.3872983346207417*g1*g2*vd*
      vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW())) -
      0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(
      g2)*Sqr(vd) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu
      ) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*
      Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + AbsSqr(Mu) + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr
      (g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = Conj(BMu);
   mass_matrix_Hpm(1,1) = mHu2 + AbsSqr(Mu) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr
      (g2)*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::Hpm, eigenvalue_error > precision * Abs
      (MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{

   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
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
   problems.flag_bad_mass(lowMSSM_info::Chi, eigenvalue_error > precision * Abs
      (MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(lowMSSM_info::Cha, eigenvalue_error > precision * Abs
      (MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
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
   problems.flag_bad_mass(lowMSSM_info::Fe, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(lowMSSM_info::Fd, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(lowMSSM_info::Fu, eigenvalue_error > precision * Abs(
      MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + Sqr(vu)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(lowMSSM_info::VWm);
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
   
   double result = Re(mHd2*vd + vd*AbsSqr(Mu) - 0.5*vu*BMu - 0.5*vu*Conj(BMu) +
      0.075*Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) - 0.075*vd*Sqr(g1)*Sqr(vu) -
      0.125*vd*Sqr(g2)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu + vu*AbsSqr(Mu) - 0.5*vd*BMu - 0.5*vd*Conj(BMu) +
      0.075*Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) - 0.075*vu*Sqr(g1)*Sqr(vd) -
      0.125*vu*Sqr(g2)*Sqr(vd));

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

   return Sqrt(Sqr(vd) + Sqr(vu));
}



std::ostream& operator<<(std::ostream& ostr, const lowMSSM_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
