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
 * @file NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the NMSSMEFTHiggs model class in the decoupling scheme
 *
 * Contains the definition of the NMSSMEFTHiggs model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme.hpp"
#include "NMSSMEFTHiggs_mass_eigenstates.hpp"
#include "NMSSMEFTHiggs_info.hpp"
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

#define CLASSNAME NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const NMSSMEFTHiggs_input_parameters& input_)
   : NMSSMEFTHiggs_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const NMSSMEFTHiggs_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<NMSSMEFTHiggs_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(NMSSMEFTHiggs_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            NMSSMEFTHiggs_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(NMSSMEFTHiggs_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            NMSSMEFTHiggs_info::g2, new_g2, get_scale());
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

void CLASSNAME::fill_from(const NMSSMEFTHiggs_mass_eigenstates& model)
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

const NMSSMEFTHiggs_physical& CLASSNAME::get_physical() const
{
   return physical;
}

NMSSMEFTHiggs_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const NMSSMEFTHiggs_physical& physical_)
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
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Cube(vd)*Sqr(g1)
      - 5*Cube(vd)*Sqr(g2) - 20*vd*AbsSqr(Lambdax)*Sqr(vS) + 10*vu*Conj(Lambdax)*
      Kappa*Sqr(vS) + 10*vu*Conj(Kappa)*Lambdax*Sqr(vS) - 20*vd*AbsSqr(Lambdax)*
      Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) + 14.142135623730951*
      vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Cube(vu)*Sqr(g1)
      - 5*Cube(vu)*Sqr(g2) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd)
      + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) + 10*vd*Conj(Lambdax)
      *Kappa*Sqr(vS) + 10*vd*Conj(Kappa)*Lambdax*Sqr(vS) + 14.142135623730951*vd*
      vS*TLambdax))/vu);
   ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) - 4*AbsSqr(Kappa)*Cube(
      vS) + 2*vd*vS*vu*Conj(Lambdax)*Kappa + 2*vd*vS*vu*Conj(Kappa)*Lambdax - 2*vS
      *AbsSqr(Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(TKappa)*Sqr(vS) - 2*vS*
      AbsSqr(Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*TKappa +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      ms2 = old_ms2;
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
           "NMSSMEFTHiggs\n"
           "========================================\n";
   NMSSMEFTHiggs_soft_parameters::print(ostr);
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
   const auto save_ms2_raii = make_raii_save(ms2);

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


   // move goldstone bosons to the front
   reorder_tree_level_masses();
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Se); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::hh); }
   if (PHYSICAL(MAh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(NMSSMEFTHiggs_info::Hpm); }
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
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,5,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,5,5>::Zero();
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
   NMSSMEFTHiggs_soft_parameters::clear();
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
   Mhh(2) = pars(28);
   MAh(0) = pars(29);
   MAh(1) = pars(30);
   MAh(2) = pars(31);
   MHpm(0) = pars(32);
   MHpm(1) = pars(33);
   MChi(0) = pars(34);
   MChi(1) = pars(35);
   MChi(2) = pars(36);
   MChi(3) = pars(37);
   MChi(4) = pars(38);
   MCha(0) = pars(39);
   MCha(1) = pars(40);
   MFe(0) = pars(41);
   MFe(1) = pars(42);
   MFe(2) = pars(43);
   MFd(0) = pars(44);
   MFd(1) = pars(45);
   MFd(2) = pars(46);
   MFu(0) = pars(47);
   MFu(1) = pars(48);
   MFu(2) = pars(49);
   MVWm = pars(50);
   MVP = pars(51);
   MVZ = pars(52);

}

const NMSSMEFTHiggs_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

NMSSMEFTHiggs_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(53);

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
   pars(28) = Mhh(2);
   pars(29) = MAh(0);
   pars(30) = MAh(1);
   pars(31) = MAh(2);
   pars(32) = MHpm(0);
   pars(33) = MHpm(1);
   pars(34) = MChi(0);
   pars(35) = MChi(1);
   pars(36) = MChi(2);
   pars(37) = MChi(3);
   pars(38) = MChi(4);
   pars(39) = MCha(0);
   pars(40) = MCha(1);
   pars(41) = MFe(0);
   pars(42) = MFe(1);
   pars(43) = MFe(2);
   pars(44) = MFd(0);
   pars(45) = MFd(1);
   pars(46) = MFd(2);
   pars(47) = MFu(0);
   pars(48) = MFu(1);
   pars(49) = MFu(2);
   pars(50) = MVWm;
   pars(51) = MVP;
   pars(52) = MVZ;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZD(0,0) = pars(53);
   ZD(0,1) = pars(54);
   ZD(0,2) = pars(55);
   ZD(0,3) = pars(56);
   ZD(0,4) = pars(57);
   ZD(0,5) = pars(58);
   ZD(1,0) = pars(59);
   ZD(1,1) = pars(60);
   ZD(1,2) = pars(61);
   ZD(1,3) = pars(62);
   ZD(1,4) = pars(63);
   ZD(1,5) = pars(64);
   ZD(2,0) = pars(65);
   ZD(2,1) = pars(66);
   ZD(2,2) = pars(67);
   ZD(2,3) = pars(68);
   ZD(2,4) = pars(69);
   ZD(2,5) = pars(70);
   ZD(3,0) = pars(71);
   ZD(3,1) = pars(72);
   ZD(3,2) = pars(73);
   ZD(3,3) = pars(74);
   ZD(3,4) = pars(75);
   ZD(3,5) = pars(76);
   ZD(4,0) = pars(77);
   ZD(4,1) = pars(78);
   ZD(4,2) = pars(79);
   ZD(4,3) = pars(80);
   ZD(4,4) = pars(81);
   ZD(4,5) = pars(82);
   ZD(5,0) = pars(83);
   ZD(5,1) = pars(84);
   ZD(5,2) = pars(85);
   ZD(5,3) = pars(86);
   ZD(5,4) = pars(87);
   ZD(5,5) = pars(88);
   ZV(0,0) = pars(89);
   ZV(0,1) = pars(90);
   ZV(0,2) = pars(91);
   ZV(1,0) = pars(92);
   ZV(1,1) = pars(93);
   ZV(1,2) = pars(94);
   ZV(2,0) = pars(95);
   ZV(2,1) = pars(96);
   ZV(2,2) = pars(97);
   ZU(0,0) = pars(98);
   ZU(0,1) = pars(99);
   ZU(0,2) = pars(100);
   ZU(0,3) = pars(101);
   ZU(0,4) = pars(102);
   ZU(0,5) = pars(103);
   ZU(1,0) = pars(104);
   ZU(1,1) = pars(105);
   ZU(1,2) = pars(106);
   ZU(1,3) = pars(107);
   ZU(1,4) = pars(108);
   ZU(1,5) = pars(109);
   ZU(2,0) = pars(110);
   ZU(2,1) = pars(111);
   ZU(2,2) = pars(112);
   ZU(2,3) = pars(113);
   ZU(2,4) = pars(114);
   ZU(2,5) = pars(115);
   ZU(3,0) = pars(116);
   ZU(3,1) = pars(117);
   ZU(3,2) = pars(118);
   ZU(3,3) = pars(119);
   ZU(3,4) = pars(120);
   ZU(3,5) = pars(121);
   ZU(4,0) = pars(122);
   ZU(4,1) = pars(123);
   ZU(4,2) = pars(124);
   ZU(4,3) = pars(125);
   ZU(4,4) = pars(126);
   ZU(4,5) = pars(127);
   ZU(5,0) = pars(128);
   ZU(5,1) = pars(129);
   ZU(5,2) = pars(130);
   ZU(5,3) = pars(131);
   ZU(5,4) = pars(132);
   ZU(5,5) = pars(133);
   ZE(0,0) = pars(134);
   ZE(0,1) = pars(135);
   ZE(0,2) = pars(136);
   ZE(0,3) = pars(137);
   ZE(0,4) = pars(138);
   ZE(0,5) = pars(139);
   ZE(1,0) = pars(140);
   ZE(1,1) = pars(141);
   ZE(1,2) = pars(142);
   ZE(1,3) = pars(143);
   ZE(1,4) = pars(144);
   ZE(1,5) = pars(145);
   ZE(2,0) = pars(146);
   ZE(2,1) = pars(147);
   ZE(2,2) = pars(148);
   ZE(2,3) = pars(149);
   ZE(2,4) = pars(150);
   ZE(2,5) = pars(151);
   ZE(3,0) = pars(152);
   ZE(3,1) = pars(153);
   ZE(3,2) = pars(154);
   ZE(3,3) = pars(155);
   ZE(3,4) = pars(156);
   ZE(3,5) = pars(157);
   ZE(4,0) = pars(158);
   ZE(4,1) = pars(159);
   ZE(4,2) = pars(160);
   ZE(4,3) = pars(161);
   ZE(4,4) = pars(162);
   ZE(4,5) = pars(163);
   ZE(5,0) = pars(164);
   ZE(5,1) = pars(165);
   ZE(5,2) = pars(166);
   ZE(5,3) = pars(167);
   ZE(5,4) = pars(168);
   ZE(5,5) = pars(169);
   ZH(0,0) = pars(170);
   ZH(0,1) = pars(171);
   ZH(0,2) = pars(172);
   ZH(1,0) = pars(173);
   ZH(1,1) = pars(174);
   ZH(1,2) = pars(175);
   ZH(2,0) = pars(176);
   ZH(2,1) = pars(177);
   ZH(2,2) = pars(178);
   ZA(0,0) = pars(179);
   ZA(0,1) = pars(180);
   ZA(0,2) = pars(181);
   ZA(1,0) = pars(182);
   ZA(1,1) = pars(183);
   ZA(1,2) = pars(184);
   ZA(2,0) = pars(185);
   ZA(2,1) = pars(186);
   ZA(2,2) = pars(187);
   ZP(0,0) = pars(188);
   ZP(0,1) = pars(189);
   ZP(1,0) = pars(190);
   ZP(1,1) = pars(191);
   ZN(0,0) = std::complex<double>(pars(192), pars(193));
   ZN(0,1) = std::complex<double>(pars(194), pars(195));
   ZN(0,2) = std::complex<double>(pars(196), pars(197));
   ZN(0,3) = std::complex<double>(pars(198), pars(199));
   ZN(0,4) = std::complex<double>(pars(200), pars(201));
   ZN(1,0) = std::complex<double>(pars(202), pars(203));
   ZN(1,1) = std::complex<double>(pars(204), pars(205));
   ZN(1,2) = std::complex<double>(pars(206), pars(207));
   ZN(1,3) = std::complex<double>(pars(208), pars(209));
   ZN(1,4) = std::complex<double>(pars(210), pars(211));
   ZN(2,0) = std::complex<double>(pars(212), pars(213));
   ZN(2,1) = std::complex<double>(pars(214), pars(215));
   ZN(2,2) = std::complex<double>(pars(216), pars(217));
   ZN(2,3) = std::complex<double>(pars(218), pars(219));
   ZN(2,4) = std::complex<double>(pars(220), pars(221));
   ZN(3,0) = std::complex<double>(pars(222), pars(223));
   ZN(3,1) = std::complex<double>(pars(224), pars(225));
   ZN(3,2) = std::complex<double>(pars(226), pars(227));
   ZN(3,3) = std::complex<double>(pars(228), pars(229));
   ZN(3,4) = std::complex<double>(pars(230), pars(231));
   ZN(4,0) = std::complex<double>(pars(232), pars(233));
   ZN(4,1) = std::complex<double>(pars(234), pars(235));
   ZN(4,2) = std::complex<double>(pars(236), pars(237));
   ZN(4,3) = std::complex<double>(pars(238), pars(239));
   ZN(4,4) = std::complex<double>(pars(240), pars(241));
   UM(0,0) = std::complex<double>(pars(242), pars(243));
   UM(0,1) = std::complex<double>(pars(244), pars(245));
   UM(1,0) = std::complex<double>(pars(246), pars(247));
   UM(1,1) = std::complex<double>(pars(248), pars(249));
   UP(0,0) = std::complex<double>(pars(250), pars(251));
   UP(0,1) = std::complex<double>(pars(252), pars(253));
   UP(1,0) = std::complex<double>(pars(254), pars(255));
   UP(1,1) = std::complex<double>(pars(256), pars(257));
   ZEL(0,0) = std::complex<double>(pars(258), pars(259));
   ZEL(0,1) = std::complex<double>(pars(260), pars(261));
   ZEL(0,2) = std::complex<double>(pars(262), pars(263));
   ZEL(1,0) = std::complex<double>(pars(264), pars(265));
   ZEL(1,1) = std::complex<double>(pars(266), pars(267));
   ZEL(1,2) = std::complex<double>(pars(268), pars(269));
   ZEL(2,0) = std::complex<double>(pars(270), pars(271));
   ZEL(2,1) = std::complex<double>(pars(272), pars(273));
   ZEL(2,2) = std::complex<double>(pars(274), pars(275));
   ZER(0,0) = std::complex<double>(pars(276), pars(277));
   ZER(0,1) = std::complex<double>(pars(278), pars(279));
   ZER(0,2) = std::complex<double>(pars(280), pars(281));
   ZER(1,0) = std::complex<double>(pars(282), pars(283));
   ZER(1,1) = std::complex<double>(pars(284), pars(285));
   ZER(1,2) = std::complex<double>(pars(286), pars(287));
   ZER(2,0) = std::complex<double>(pars(288), pars(289));
   ZER(2,1) = std::complex<double>(pars(290), pars(291));
   ZER(2,2) = std::complex<double>(pars(292), pars(293));
   ZDL(0,0) = std::complex<double>(pars(294), pars(295));
   ZDL(0,1) = std::complex<double>(pars(296), pars(297));
   ZDL(0,2) = std::complex<double>(pars(298), pars(299));
   ZDL(1,0) = std::complex<double>(pars(300), pars(301));
   ZDL(1,1) = std::complex<double>(pars(302), pars(303));
   ZDL(1,2) = std::complex<double>(pars(304), pars(305));
   ZDL(2,0) = std::complex<double>(pars(306), pars(307));
   ZDL(2,1) = std::complex<double>(pars(308), pars(309));
   ZDL(2,2) = std::complex<double>(pars(310), pars(311));
   ZDR(0,0) = std::complex<double>(pars(312), pars(313));
   ZDR(0,1) = std::complex<double>(pars(314), pars(315));
   ZDR(0,2) = std::complex<double>(pars(316), pars(317));
   ZDR(1,0) = std::complex<double>(pars(318), pars(319));
   ZDR(1,1) = std::complex<double>(pars(320), pars(321));
   ZDR(1,2) = std::complex<double>(pars(322), pars(323));
   ZDR(2,0) = std::complex<double>(pars(324), pars(325));
   ZDR(2,1) = std::complex<double>(pars(326), pars(327));
   ZDR(2,2) = std::complex<double>(pars(328), pars(329));
   ZUL(0,0) = std::complex<double>(pars(330), pars(331));
   ZUL(0,1) = std::complex<double>(pars(332), pars(333));
   ZUL(0,2) = std::complex<double>(pars(334), pars(335));
   ZUL(1,0) = std::complex<double>(pars(336), pars(337));
   ZUL(1,1) = std::complex<double>(pars(338), pars(339));
   ZUL(1,2) = std::complex<double>(pars(340), pars(341));
   ZUL(2,0) = std::complex<double>(pars(342), pars(343));
   ZUL(2,1) = std::complex<double>(pars(344), pars(345));
   ZUL(2,2) = std::complex<double>(pars(346), pars(347));
   ZUR(0,0) = std::complex<double>(pars(348), pars(349));
   ZUR(0,1) = std::complex<double>(pars(350), pars(351));
   ZUR(0,2) = std::complex<double>(pars(352), pars(353));
   ZUR(1,0) = std::complex<double>(pars(354), pars(355));
   ZUR(1,1) = std::complex<double>(pars(356), pars(357));
   ZUR(1,2) = std::complex<double>(pars(358), pars(359));
   ZUR(2,0) = std::complex<double>(pars(360), pars(361));
   ZUR(2,1) = std::complex<double>(pars(362), pars(363));
   ZUR(2,2) = std::complex<double>(pars(364), pars(365));
   ZZ(0,0) = pars(366);
   ZZ(0,1) = pars(367);
   ZZ(1,0) = pars(368);
   ZZ(1,1) = pars(369);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(370);

   pars(53) = ZD(0,0);
   pars(54) = ZD(0,1);
   pars(55) = ZD(0,2);
   pars(56) = ZD(0,3);
   pars(57) = ZD(0,4);
   pars(58) = ZD(0,5);
   pars(59) = ZD(1,0);
   pars(60) = ZD(1,1);
   pars(61) = ZD(1,2);
   pars(62) = ZD(1,3);
   pars(63) = ZD(1,4);
   pars(64) = ZD(1,5);
   pars(65) = ZD(2,0);
   pars(66) = ZD(2,1);
   pars(67) = ZD(2,2);
   pars(68) = ZD(2,3);
   pars(69) = ZD(2,4);
   pars(70) = ZD(2,5);
   pars(71) = ZD(3,0);
   pars(72) = ZD(3,1);
   pars(73) = ZD(3,2);
   pars(74) = ZD(3,3);
   pars(75) = ZD(3,4);
   pars(76) = ZD(3,5);
   pars(77) = ZD(4,0);
   pars(78) = ZD(4,1);
   pars(79) = ZD(4,2);
   pars(80) = ZD(4,3);
   pars(81) = ZD(4,4);
   pars(82) = ZD(4,5);
   pars(83) = ZD(5,0);
   pars(84) = ZD(5,1);
   pars(85) = ZD(5,2);
   pars(86) = ZD(5,3);
   pars(87) = ZD(5,4);
   pars(88) = ZD(5,5);
   pars(89) = ZV(0,0);
   pars(90) = ZV(0,1);
   pars(91) = ZV(0,2);
   pars(92) = ZV(1,0);
   pars(93) = ZV(1,1);
   pars(94) = ZV(1,2);
   pars(95) = ZV(2,0);
   pars(96) = ZV(2,1);
   pars(97) = ZV(2,2);
   pars(98) = ZU(0,0);
   pars(99) = ZU(0,1);
   pars(100) = ZU(0,2);
   pars(101) = ZU(0,3);
   pars(102) = ZU(0,4);
   pars(103) = ZU(0,5);
   pars(104) = ZU(1,0);
   pars(105) = ZU(1,1);
   pars(106) = ZU(1,2);
   pars(107) = ZU(1,3);
   pars(108) = ZU(1,4);
   pars(109) = ZU(1,5);
   pars(110) = ZU(2,0);
   pars(111) = ZU(2,1);
   pars(112) = ZU(2,2);
   pars(113) = ZU(2,3);
   pars(114) = ZU(2,4);
   pars(115) = ZU(2,5);
   pars(116) = ZU(3,0);
   pars(117) = ZU(3,1);
   pars(118) = ZU(3,2);
   pars(119) = ZU(3,3);
   pars(120) = ZU(3,4);
   pars(121) = ZU(3,5);
   pars(122) = ZU(4,0);
   pars(123) = ZU(4,1);
   pars(124) = ZU(4,2);
   pars(125) = ZU(4,3);
   pars(126) = ZU(4,4);
   pars(127) = ZU(4,5);
   pars(128) = ZU(5,0);
   pars(129) = ZU(5,1);
   pars(130) = ZU(5,2);
   pars(131) = ZU(5,3);
   pars(132) = ZU(5,4);
   pars(133) = ZU(5,5);
   pars(134) = ZE(0,0);
   pars(135) = ZE(0,1);
   pars(136) = ZE(0,2);
   pars(137) = ZE(0,3);
   pars(138) = ZE(0,4);
   pars(139) = ZE(0,5);
   pars(140) = ZE(1,0);
   pars(141) = ZE(1,1);
   pars(142) = ZE(1,2);
   pars(143) = ZE(1,3);
   pars(144) = ZE(1,4);
   pars(145) = ZE(1,5);
   pars(146) = ZE(2,0);
   pars(147) = ZE(2,1);
   pars(148) = ZE(2,2);
   pars(149) = ZE(2,3);
   pars(150) = ZE(2,4);
   pars(151) = ZE(2,5);
   pars(152) = ZE(3,0);
   pars(153) = ZE(3,1);
   pars(154) = ZE(3,2);
   pars(155) = ZE(3,3);
   pars(156) = ZE(3,4);
   pars(157) = ZE(3,5);
   pars(158) = ZE(4,0);
   pars(159) = ZE(4,1);
   pars(160) = ZE(4,2);
   pars(161) = ZE(4,3);
   pars(162) = ZE(4,4);
   pars(163) = ZE(4,5);
   pars(164) = ZE(5,0);
   pars(165) = ZE(5,1);
   pars(166) = ZE(5,2);
   pars(167) = ZE(5,3);
   pars(168) = ZE(5,4);
   pars(169) = ZE(5,5);
   pars(170) = ZH(0,0);
   pars(171) = ZH(0,1);
   pars(172) = ZH(0,2);
   pars(173) = ZH(1,0);
   pars(174) = ZH(1,1);
   pars(175) = ZH(1,2);
   pars(176) = ZH(2,0);
   pars(177) = ZH(2,1);
   pars(178) = ZH(2,2);
   pars(179) = ZA(0,0);
   pars(180) = ZA(0,1);
   pars(181) = ZA(0,2);
   pars(182) = ZA(1,0);
   pars(183) = ZA(1,1);
   pars(184) = ZA(1,2);
   pars(185) = ZA(2,0);
   pars(186) = ZA(2,1);
   pars(187) = ZA(2,2);
   pars(188) = ZP(0,0);
   pars(189) = ZP(0,1);
   pars(190) = ZP(1,0);
   pars(191) = ZP(1,1);
   pars(192) = Re(ZN(0,0));
   pars(193) = Im(ZN(0,0));
   pars(194) = Re(ZN(0,1));
   pars(195) = Im(ZN(0,1));
   pars(196) = Re(ZN(0,2));
   pars(197) = Im(ZN(0,2));
   pars(198) = Re(ZN(0,3));
   pars(199) = Im(ZN(0,3));
   pars(200) = Re(ZN(0,4));
   pars(201) = Im(ZN(0,4));
   pars(202) = Re(ZN(1,0));
   pars(203) = Im(ZN(1,0));
   pars(204) = Re(ZN(1,1));
   pars(205) = Im(ZN(1,1));
   pars(206) = Re(ZN(1,2));
   pars(207) = Im(ZN(1,2));
   pars(208) = Re(ZN(1,3));
   pars(209) = Im(ZN(1,3));
   pars(210) = Re(ZN(1,4));
   pars(211) = Im(ZN(1,4));
   pars(212) = Re(ZN(2,0));
   pars(213) = Im(ZN(2,0));
   pars(214) = Re(ZN(2,1));
   pars(215) = Im(ZN(2,1));
   pars(216) = Re(ZN(2,2));
   pars(217) = Im(ZN(2,2));
   pars(218) = Re(ZN(2,3));
   pars(219) = Im(ZN(2,3));
   pars(220) = Re(ZN(2,4));
   pars(221) = Im(ZN(2,4));
   pars(222) = Re(ZN(3,0));
   pars(223) = Im(ZN(3,0));
   pars(224) = Re(ZN(3,1));
   pars(225) = Im(ZN(3,1));
   pars(226) = Re(ZN(3,2));
   pars(227) = Im(ZN(3,2));
   pars(228) = Re(ZN(3,3));
   pars(229) = Im(ZN(3,3));
   pars(230) = Re(ZN(3,4));
   pars(231) = Im(ZN(3,4));
   pars(232) = Re(ZN(4,0));
   pars(233) = Im(ZN(4,0));
   pars(234) = Re(ZN(4,1));
   pars(235) = Im(ZN(4,1));
   pars(236) = Re(ZN(4,2));
   pars(237) = Im(ZN(4,2));
   pars(238) = Re(ZN(4,3));
   pars(239) = Im(ZN(4,3));
   pars(240) = Re(ZN(4,4));
   pars(241) = Im(ZN(4,4));
   pars(242) = Re(UM(0,0));
   pars(243) = Im(UM(0,0));
   pars(244) = Re(UM(0,1));
   pars(245) = Im(UM(0,1));
   pars(246) = Re(UM(1,0));
   pars(247) = Im(UM(1,0));
   pars(248) = Re(UM(1,1));
   pars(249) = Im(UM(1,1));
   pars(250) = Re(UP(0,0));
   pars(251) = Im(UP(0,0));
   pars(252) = Re(UP(0,1));
   pars(253) = Im(UP(0,1));
   pars(254) = Re(UP(1,0));
   pars(255) = Im(UP(1,0));
   pars(256) = Re(UP(1,1));
   pars(257) = Im(UP(1,1));
   pars(258) = Re(ZEL(0,0));
   pars(259) = Im(ZEL(0,0));
   pars(260) = Re(ZEL(0,1));
   pars(261) = Im(ZEL(0,1));
   pars(262) = Re(ZEL(0,2));
   pars(263) = Im(ZEL(0,2));
   pars(264) = Re(ZEL(1,0));
   pars(265) = Im(ZEL(1,0));
   pars(266) = Re(ZEL(1,1));
   pars(267) = Im(ZEL(1,1));
   pars(268) = Re(ZEL(1,2));
   pars(269) = Im(ZEL(1,2));
   pars(270) = Re(ZEL(2,0));
   pars(271) = Im(ZEL(2,0));
   pars(272) = Re(ZEL(2,1));
   pars(273) = Im(ZEL(2,1));
   pars(274) = Re(ZEL(2,2));
   pars(275) = Im(ZEL(2,2));
   pars(276) = Re(ZER(0,0));
   pars(277) = Im(ZER(0,0));
   pars(278) = Re(ZER(0,1));
   pars(279) = Im(ZER(0,1));
   pars(280) = Re(ZER(0,2));
   pars(281) = Im(ZER(0,2));
   pars(282) = Re(ZER(1,0));
   pars(283) = Im(ZER(1,0));
   pars(284) = Re(ZER(1,1));
   pars(285) = Im(ZER(1,1));
   pars(286) = Re(ZER(1,2));
   pars(287) = Im(ZER(1,2));
   pars(288) = Re(ZER(2,0));
   pars(289) = Im(ZER(2,0));
   pars(290) = Re(ZER(2,1));
   pars(291) = Im(ZER(2,1));
   pars(292) = Re(ZER(2,2));
   pars(293) = Im(ZER(2,2));
   pars(294) = Re(ZDL(0,0));
   pars(295) = Im(ZDL(0,0));
   pars(296) = Re(ZDL(0,1));
   pars(297) = Im(ZDL(0,1));
   pars(298) = Re(ZDL(0,2));
   pars(299) = Im(ZDL(0,2));
   pars(300) = Re(ZDL(1,0));
   pars(301) = Im(ZDL(1,0));
   pars(302) = Re(ZDL(1,1));
   pars(303) = Im(ZDL(1,1));
   pars(304) = Re(ZDL(1,2));
   pars(305) = Im(ZDL(1,2));
   pars(306) = Re(ZDL(2,0));
   pars(307) = Im(ZDL(2,0));
   pars(308) = Re(ZDL(2,1));
   pars(309) = Im(ZDL(2,1));
   pars(310) = Re(ZDL(2,2));
   pars(311) = Im(ZDL(2,2));
   pars(312) = Re(ZDR(0,0));
   pars(313) = Im(ZDR(0,0));
   pars(314) = Re(ZDR(0,1));
   pars(315) = Im(ZDR(0,1));
   pars(316) = Re(ZDR(0,2));
   pars(317) = Im(ZDR(0,2));
   pars(318) = Re(ZDR(1,0));
   pars(319) = Im(ZDR(1,0));
   pars(320) = Re(ZDR(1,1));
   pars(321) = Im(ZDR(1,1));
   pars(322) = Re(ZDR(1,2));
   pars(323) = Im(ZDR(1,2));
   pars(324) = Re(ZDR(2,0));
   pars(325) = Im(ZDR(2,0));
   pars(326) = Re(ZDR(2,1));
   pars(327) = Im(ZDR(2,1));
   pars(328) = Re(ZDR(2,2));
   pars(329) = Im(ZDR(2,2));
   pars(330) = Re(ZUL(0,0));
   pars(331) = Im(ZUL(0,0));
   pars(332) = Re(ZUL(0,1));
   pars(333) = Im(ZUL(0,1));
   pars(334) = Re(ZUL(0,2));
   pars(335) = Im(ZUL(0,2));
   pars(336) = Re(ZUL(1,0));
   pars(337) = Im(ZUL(1,0));
   pars(338) = Re(ZUL(1,1));
   pars(339) = Im(ZUL(1,1));
   pars(340) = Re(ZUL(1,2));
   pars(341) = Im(ZUL(1,2));
   pars(342) = Re(ZUL(2,0));
   pars(343) = Im(ZUL(2,0));
   pars(344) = Re(ZUL(2,1));
   pars(345) = Im(ZUL(2,1));
   pars(346) = Re(ZUL(2,2));
   pars(347) = Im(ZUL(2,2));
   pars(348) = Re(ZUR(0,0));
   pars(349) = Im(ZUR(0,0));
   pars(350) = Re(ZUR(0,1));
   pars(351) = Im(ZUR(0,1));
   pars(352) = Re(ZUR(0,2));
   pars(353) = Im(ZUR(0,2));
   pars(354) = Re(ZUR(1,0));
   pars(355) = Im(ZUR(1,0));
   pars(356) = Re(ZUR(1,1));
   pars(357) = Im(ZUR(1,1));
   pars(358) = Re(ZUR(1,2));
   pars(359) = Im(ZUR(1,2));
   pars(360) = Re(ZUR(2,0));
   pars(361) = Im(ZUR(2,0));
   pars(362) = Re(ZUR(2,1));
   pars(363) = Im(ZUR(2,1));
   pars(364) = Re(ZUR(2,2));
   pars(365) = Im(ZUR(2,2));
   pars(366) = ZZ(0,0);
   pars(367) = ZZ(0,1);
   pars(368) = ZZ(1,0);
   pars(369) = ZZ(1,1);


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

Eigen::Array<double,2,1> CLASSNAME::get_MPseudoscalarHiggs() const
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
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vS*vu*Conj(
      Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) - 0.5*vS*vu*Conj(
      Yd(1,0))*Lambdax;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) - 0.5*vS*vu*Conj(
      Yd(2,0))*Lambdax;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) + Conj(
      Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) - 0.5*vS*vu*Conj(
      Yd(0,1))*Lambdax;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vS*vu*Conj(
      Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) - 0.5*vS*vu*Conj(
      Yd(2,1))*Lambdax;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) - 0.5*vS*vu*Conj(
      Yd(0,2))*Lambdax;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) - 0.5*vS*vu*Conj(
      Yd(1,2))*Lambdax;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vS*vu*Conj(
      Yd(2,2))*Lambdax;
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Sd, eigenvalue_error > precision
      * Abs(MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Sd);
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Sv, eigenvalue_error > precision
      * Abs(MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Sv);
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
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vS*Conj(
      Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) - 0.5*vd*vS*Conj(
      Yu(1,0))*Lambdax;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) - 0.5*vd*vS*Conj(
      Yu(2,0))*Lambdax;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) - 0.5*vd*vS*Conj(
      Yu(0,1))*Lambdax;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vS*Conj(
      Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) - 0.5*vd*vS*Conj(
      Yu(2,1))*Lambdax;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) - 0.5*vd*vS*Conj(
      Yu(0,2))*Lambdax;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) - 0.5*vd*vS*Conj(
      Yu(1,2))*Lambdax;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vS*Conj(
      Yu(2,2))*Lambdax;
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Su, eigenvalue_error > precision
      * Abs(MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Su);
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
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vS*vu*Conj(
      Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) - 0.5*vS*vu*Conj(
      Ye(1,0))*Lambdax;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) - 0.5*vS*vu*Conj(
      Ye(2,0))*Lambdax;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) + Conj(
      Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) - 0.5*vS*vu*Conj(
      Ye(0,1))*Lambdax;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vS*vu*Conj(
      Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) - 0.5*vS*vu*Conj(
      Ye(2,1))*Lambdax;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) - 0.5*vS*vu*Conj(
      Ye(0,2))*Lambdax;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) - 0.5*vS*vu*Conj(
      Ye(1,2))*Lambdax;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vS*vu*Conj(
      Ye(2,2))*Lambdax;
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Se, eigenvalue_error > precision
      * Abs(MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)
      *Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vS*Conj(
      TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) - 0.25*Conj(Lambdax)*
      Kappa*Sqr(vS) - 0.25*Conj(Kappa)*Lambdax*Sqr(vS) - 0.35355339059327373*vS
      *TLambdax;
   mass_matrix_hh(0,2) = vd*vS*AbsSqr(Lambdax) - 0.35355339059327373*vu*Conj(
      TLambdax) - 0.5*vS*vu*Conj(Lambdax)*Kappa - 0.5*vS*vu*Conj(Kappa)*Lambdax
       - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.225*Sqr(g1
      )*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = vS*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*Conj(
      TLambdax) - 0.5*vd*vS*Conj(Lambdax)*Kappa - 0.5*vd*vS*Conj(Kappa)*Lambdax
       - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.7071067811865475*vS*Conj(TKappa) - 0.5*vd*vu*
      Conj(Lambdax)*Kappa - 0.5*vd*vu*Conj(Kappa)*Lambdax + 0.5*AbsSqr(Lambdax)
      *Sqr(vd) + 3*AbsSqr(Kappa)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) +
      0.7071067811865475*vS*TKappa;

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(NMSSMEFTHiggs_info::hh, eigenvalue_error > precision
      * Abs(Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(
      ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*
      AbsSqr(Lambdax)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr
      (vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW())) +
      0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vS*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) + 0.25*Conj(
      Lambdax)*Kappa*Sqr(vS) + 0.25*Conj(Kappa)*Lambdax*Sqr(vS) - 0.25*vd*vu*
      Sqr(g2)*Sqr(Cos(ThetaW())) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW())) +
      0.35355339059327373*vS*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) - 0.5*vS*vu*Conj
      (Lambdax)*Kappa - 0.5*vS*vu*Conj(Kappa)*Lambdax + 0.35355339059327373*vu*
      TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW
      ())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) - 0.5*vd*vS*Conj
      (Lambdax)*Kappa - 0.5*vd*vS*Conj(Kappa)*Lambdax + 0.35355339059327373*vd*
      TLambdax;
   mass_matrix_Ah(2,2) = ms2 - 0.7071067811865475*vS*Conj(TKappa) + 0.5*vd*vu*
      Conj(Lambdax)*Kappa + 0.5*vd*vu*Conj(Kappa)*Lambdax + 0.5*AbsSqr(Lambdax)
      *Sqr(vd) + AbsSqr(Kappa)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) -
      0.7071067811865475*vS*TKappa;

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Ah, eigenvalue_error > precision
      * Abs(MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd)
      + 0.5*AbsSqr(Lambdax)*Sqr(vS) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr
      (vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475*vS*
      Conj(TLambdax) + 0.5*Conj(Lambdax)*Kappa*Sqr(vS);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr
      (vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Hpm, eigenvalue_error > precision
       * Abs(MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,5,5> CLASSNAME::get_mass_matrix_Chi() const
{

   Eigen::Matrix<double,5,5> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(0,4) = 0;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.7071067811865475*vS*Lambdax;
   mass_matrix_Chi(2,4) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(4,4) = 1.4142135623730951*vS*Kappa;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Chi, eigenvalue_error > precision
       * Abs(MChi(0)));
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
   mass_matrix_Cha(1,1) = 0.7071067811865475*vS*Lambdax;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Cha, eigenvalue_error > precision
       * Abs(MCha(0)));
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Fe, eigenvalue_error > precision
      * Abs(MFe(0)));
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Fd, eigenvalue_error > precision
      * Abs(MFd(0)));
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
   problems.flag_bad_mass(NMSSMEFTHiggs_info::Fu, eigenvalue_error > precision
      * Abs(MFu(0)));
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
      problems.flag_running_tachyon(NMSSMEFTHiggs_info::VWm);
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
   
   double result = Re(mHd2*vd - 0.35355339059327373*vS*vu*Conj(TLambdax) + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vS) -
      0.25*vu*Conj(Lambdax)*Kappa*Sqr(vS) - 0.25*vu*Conj(Kappa)*Lambdax*Sqr(vS) +
      0.5*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(g2)
      *Sqr(vu) - 0.35355339059327373*vS*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu - 0.35355339059327373*vd*vS*Conj(TLambdax) + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vd) -
      0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)
      *Sqr(vS) - 0.25*vd*Conj(Lambdax)*Kappa*Sqr(vS) - 0.25*vd*Conj(Kappa)*Lambdax
      *Sqr(vS) - 0.35355339059327373*vd*vS*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(ms2*vS - 0.35355339059327373*vd*vu*Conj(TLambdax) + AbsSqr(
      Kappa)*Cube(vS) - 0.5*vd*vS*vu*Conj(Lambdax)*Kappa - 0.5*vd*vS*vu*Conj(Kappa
      )*Lambdax + 0.5*vS*AbsSqr(Lambdax)*Sqr(vd) + 0.35355339059327373*Conj(TKappa
      )*Sqr(vS) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vu) + 0.35355339059327373*Sqr(vS)*
      TKappa - 0.35355339059327373*vd*vu*TLambdax);

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



std::ostream& operator<<(std::ostream& ostr, const NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
