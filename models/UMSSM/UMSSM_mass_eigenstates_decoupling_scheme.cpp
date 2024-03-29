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
 * @file UMSSM_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the UMSSM model class in the decoupling scheme
 *
 * Contains the definition of the UMSSM model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "UMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "UMSSM_mass_eigenstates.hpp"
#include "UMSSM_info.hpp"
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

#define CLASSNAME UMSSM_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const UMSSM_input_parameters& input_)
   : UMSSM_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const UMSSM_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<UMSSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<UMSSM_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(UMSSM_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            UMSSM_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(UMSSM_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            UMSSM_info::g2, new_g2, get_scale());
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

void CLASSNAME::fill_from(const UMSSM_mass_eigenstates& model)
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
   MFv = OTHER(MFv);
   ZVL = OTHER(ZVL);
   ZVR = OTHER(ZVR);
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
   MVZp = OTHER(MVZp);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const UMSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

UMSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const UMSSM_physical& physical_)
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

   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Cube(vd)*Sqr(g1)
      - 5*Cube(vd)*Sqr(g2) - 20*Cube(vd)*Sqr(gp)*Sqr(QHd) - 20*vd*AbsSqr(Lambdax)*
      Sqr(vS) - 20*QHd*Qs*vd*Sqr(gp)*Sqr(vS) - 20*vd*AbsSqr(Lambdax)*Sqr(vu) + 3*
      vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) - 20*QHd*QHu*vd*Sqr(gp)*Sqr(vu) +
      14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Cube(vu)*Sqr(g1)
      - 5*Cube(vu)*Sqr(g2) - 20*Cube(vu)*Sqr(gp)*Sqr(QHu) - 20*vu*AbsSqr(Lambdax)*
      Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*QHd*QHu*vu*Sqr(gp
      )*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) - 20*QHu*Qs*vu*Sqr(gp)*Sqr(vS) +
      14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) - 2*Cube(vS)*Sqr(gp)*
      Sqr(Qs) - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 2*QHd*Qs*vS*Sqr(gp)*Sqr(vd) - 2*vS*
      AbsSqr(Lambdax)*Sqr(vu) - 2*QHu*Qs*vS*Sqr(gp)*Sqr(vu) + 1.4142135623730951*
      vd*vu*TLambdax))/vS);

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
           "UMSSM\n"
           "========================================\n";
   UMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MVZp = " << MVZp << '\n';

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
   ostr << "ZVL = " << ZVL << '\n';
   ostr << "ZVR = " << ZVR << '\n';
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

   calculate_MVPVZVZp();
   calculate_MVWm();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MFv();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
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
   PHYSICAL(MFv) = MFv;
   PHYSICAL(ZVL) = ZVL;
   PHYSICAL(ZVR) = ZVR;
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
   PHYSICAL(MVZp) = MVZp;
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
   move_goldstone_to(1, MVZp, MAh, ZA);
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
   move_goldstone_to(1, MVZp, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Sd); }
   if (PHYSICAL(MSv).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Se); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Hpm); }
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
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,6,1>::Zero();
   ZV = Eigen::Matrix<double,6,6>::Zero();
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
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MFv = Eigen::Matrix<double,3,1>::Zero();
   ZVL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZVR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
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
   MVZp = 0.;

   PhaseGlu = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   UMSSM_soft_parameters::clear();
   clear_tree_level_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MSd(0) = pars(2);
   MSd(1) = pars(3);
   MSd(2) = pars(4);
   MSd(3) = pars(5);
   MSd(4) = pars(6);
   MSd(5) = pars(7);
   MSv(0) = pars(8);
   MSv(1) = pars(9);
   MSv(2) = pars(10);
   MSv(3) = pars(11);
   MSv(4) = pars(12);
   MSv(5) = pars(13);
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
   MChi(5) = pars(39);
   MFv(0) = pars(40);
   MFv(1) = pars(41);
   MFv(2) = pars(42);
   MCha(0) = pars(43);
   MCha(1) = pars(44);
   MFe(0) = pars(45);
   MFe(1) = pars(46);
   MFe(2) = pars(47);
   MFd(0) = pars(48);
   MFd(1) = pars(49);
   MFd(2) = pars(50);
   MFu(0) = pars(51);
   MFu(1) = pars(52);
   MFu(2) = pars(53);
   MVWm = pars(54);
   MVP = pars(55);
   MVZ = pars(56);
   MVZp = pars(57);

}

const UMSSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

UMSSM_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(58);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MSd(0);
   pars(3) = MSd(1);
   pars(4) = MSd(2);
   pars(5) = MSd(3);
   pars(6) = MSd(4);
   pars(7) = MSd(5);
   pars(8) = MSv(0);
   pars(9) = MSv(1);
   pars(10) = MSv(2);
   pars(11) = MSv(3);
   pars(12) = MSv(4);
   pars(13) = MSv(5);
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
   pars(39) = MChi(5);
   pars(40) = MFv(0);
   pars(41) = MFv(1);
   pars(42) = MFv(2);
   pars(43) = MCha(0);
   pars(44) = MCha(1);
   pars(45) = MFe(0);
   pars(46) = MFe(1);
   pars(47) = MFe(2);
   pars(48) = MFd(0);
   pars(49) = MFd(1);
   pars(50) = MFd(2);
   pars(51) = MFu(0);
   pars(52) = MFu(1);
   pars(53) = MFu(2);
   pars(54) = MVWm;
   pars(55) = MVP;
   pars(56) = MVZ;
   pars(57) = MVZp;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZD(0,0) = pars(58);
   ZD(0,1) = pars(59);
   ZD(0,2) = pars(60);
   ZD(0,3) = pars(61);
   ZD(0,4) = pars(62);
   ZD(0,5) = pars(63);
   ZD(1,0) = pars(64);
   ZD(1,1) = pars(65);
   ZD(1,2) = pars(66);
   ZD(1,3) = pars(67);
   ZD(1,4) = pars(68);
   ZD(1,5) = pars(69);
   ZD(2,0) = pars(70);
   ZD(2,1) = pars(71);
   ZD(2,2) = pars(72);
   ZD(2,3) = pars(73);
   ZD(2,4) = pars(74);
   ZD(2,5) = pars(75);
   ZD(3,0) = pars(76);
   ZD(3,1) = pars(77);
   ZD(3,2) = pars(78);
   ZD(3,3) = pars(79);
   ZD(3,4) = pars(80);
   ZD(3,5) = pars(81);
   ZD(4,0) = pars(82);
   ZD(4,1) = pars(83);
   ZD(4,2) = pars(84);
   ZD(4,3) = pars(85);
   ZD(4,4) = pars(86);
   ZD(4,5) = pars(87);
   ZD(5,0) = pars(88);
   ZD(5,1) = pars(89);
   ZD(5,2) = pars(90);
   ZD(5,3) = pars(91);
   ZD(5,4) = pars(92);
   ZD(5,5) = pars(93);
   ZV(0,0) = pars(94);
   ZV(0,1) = pars(95);
   ZV(0,2) = pars(96);
   ZV(0,3) = pars(97);
   ZV(0,4) = pars(98);
   ZV(0,5) = pars(99);
   ZV(1,0) = pars(100);
   ZV(1,1) = pars(101);
   ZV(1,2) = pars(102);
   ZV(1,3) = pars(103);
   ZV(1,4) = pars(104);
   ZV(1,5) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZV(2,3) = pars(109);
   ZV(2,4) = pars(110);
   ZV(2,5) = pars(111);
   ZV(3,0) = pars(112);
   ZV(3,1) = pars(113);
   ZV(3,2) = pars(114);
   ZV(3,3) = pars(115);
   ZV(3,4) = pars(116);
   ZV(3,5) = pars(117);
   ZV(4,0) = pars(118);
   ZV(4,1) = pars(119);
   ZV(4,2) = pars(120);
   ZV(4,3) = pars(121);
   ZV(4,4) = pars(122);
   ZV(4,5) = pars(123);
   ZV(5,0) = pars(124);
   ZV(5,1) = pars(125);
   ZV(5,2) = pars(126);
   ZV(5,3) = pars(127);
   ZV(5,4) = pars(128);
   ZV(5,5) = pars(129);
   ZU(0,0) = pars(130);
   ZU(0,1) = pars(131);
   ZU(0,2) = pars(132);
   ZU(0,3) = pars(133);
   ZU(0,4) = pars(134);
   ZU(0,5) = pars(135);
   ZU(1,0) = pars(136);
   ZU(1,1) = pars(137);
   ZU(1,2) = pars(138);
   ZU(1,3) = pars(139);
   ZU(1,4) = pars(140);
   ZU(1,5) = pars(141);
   ZU(2,0) = pars(142);
   ZU(2,1) = pars(143);
   ZU(2,2) = pars(144);
   ZU(2,3) = pars(145);
   ZU(2,4) = pars(146);
   ZU(2,5) = pars(147);
   ZU(3,0) = pars(148);
   ZU(3,1) = pars(149);
   ZU(3,2) = pars(150);
   ZU(3,3) = pars(151);
   ZU(3,4) = pars(152);
   ZU(3,5) = pars(153);
   ZU(4,0) = pars(154);
   ZU(4,1) = pars(155);
   ZU(4,2) = pars(156);
   ZU(4,3) = pars(157);
   ZU(4,4) = pars(158);
   ZU(4,5) = pars(159);
   ZU(5,0) = pars(160);
   ZU(5,1) = pars(161);
   ZU(5,2) = pars(162);
   ZU(5,3) = pars(163);
   ZU(5,4) = pars(164);
   ZU(5,5) = pars(165);
   ZE(0,0) = pars(166);
   ZE(0,1) = pars(167);
   ZE(0,2) = pars(168);
   ZE(0,3) = pars(169);
   ZE(0,4) = pars(170);
   ZE(0,5) = pars(171);
   ZE(1,0) = pars(172);
   ZE(1,1) = pars(173);
   ZE(1,2) = pars(174);
   ZE(1,3) = pars(175);
   ZE(1,4) = pars(176);
   ZE(1,5) = pars(177);
   ZE(2,0) = pars(178);
   ZE(2,1) = pars(179);
   ZE(2,2) = pars(180);
   ZE(2,3) = pars(181);
   ZE(2,4) = pars(182);
   ZE(2,5) = pars(183);
   ZE(3,0) = pars(184);
   ZE(3,1) = pars(185);
   ZE(3,2) = pars(186);
   ZE(3,3) = pars(187);
   ZE(3,4) = pars(188);
   ZE(3,5) = pars(189);
   ZE(4,0) = pars(190);
   ZE(4,1) = pars(191);
   ZE(4,2) = pars(192);
   ZE(4,3) = pars(193);
   ZE(4,4) = pars(194);
   ZE(4,5) = pars(195);
   ZE(5,0) = pars(196);
   ZE(5,1) = pars(197);
   ZE(5,2) = pars(198);
   ZE(5,3) = pars(199);
   ZE(5,4) = pars(200);
   ZE(5,5) = pars(201);
   ZH(0,0) = pars(202);
   ZH(0,1) = pars(203);
   ZH(0,2) = pars(204);
   ZH(1,0) = pars(205);
   ZH(1,1) = pars(206);
   ZH(1,2) = pars(207);
   ZH(2,0) = pars(208);
   ZH(2,1) = pars(209);
   ZH(2,2) = pars(210);
   ZA(0,0) = pars(211);
   ZA(0,1) = pars(212);
   ZA(0,2) = pars(213);
   ZA(1,0) = pars(214);
   ZA(1,1) = pars(215);
   ZA(1,2) = pars(216);
   ZA(2,0) = pars(217);
   ZA(2,1) = pars(218);
   ZA(2,2) = pars(219);
   ZP(0,0) = pars(220);
   ZP(0,1) = pars(221);
   ZP(1,0) = pars(222);
   ZP(1,1) = pars(223);
   ZN(0,0) = std::complex<double>(pars(224), pars(225));
   ZN(0,1) = std::complex<double>(pars(226), pars(227));
   ZN(0,2) = std::complex<double>(pars(228), pars(229));
   ZN(0,3) = std::complex<double>(pars(230), pars(231));
   ZN(0,4) = std::complex<double>(pars(232), pars(233));
   ZN(0,5) = std::complex<double>(pars(234), pars(235));
   ZN(1,0) = std::complex<double>(pars(236), pars(237));
   ZN(1,1) = std::complex<double>(pars(238), pars(239));
   ZN(1,2) = std::complex<double>(pars(240), pars(241));
   ZN(1,3) = std::complex<double>(pars(242), pars(243));
   ZN(1,4) = std::complex<double>(pars(244), pars(245));
   ZN(1,5) = std::complex<double>(pars(246), pars(247));
   ZN(2,0) = std::complex<double>(pars(248), pars(249));
   ZN(2,1) = std::complex<double>(pars(250), pars(251));
   ZN(2,2) = std::complex<double>(pars(252), pars(253));
   ZN(2,3) = std::complex<double>(pars(254), pars(255));
   ZN(2,4) = std::complex<double>(pars(256), pars(257));
   ZN(2,5) = std::complex<double>(pars(258), pars(259));
   ZN(3,0) = std::complex<double>(pars(260), pars(261));
   ZN(3,1) = std::complex<double>(pars(262), pars(263));
   ZN(3,2) = std::complex<double>(pars(264), pars(265));
   ZN(3,3) = std::complex<double>(pars(266), pars(267));
   ZN(3,4) = std::complex<double>(pars(268), pars(269));
   ZN(3,5) = std::complex<double>(pars(270), pars(271));
   ZN(4,0) = std::complex<double>(pars(272), pars(273));
   ZN(4,1) = std::complex<double>(pars(274), pars(275));
   ZN(4,2) = std::complex<double>(pars(276), pars(277));
   ZN(4,3) = std::complex<double>(pars(278), pars(279));
   ZN(4,4) = std::complex<double>(pars(280), pars(281));
   ZN(4,5) = std::complex<double>(pars(282), pars(283));
   ZN(5,0) = std::complex<double>(pars(284), pars(285));
   ZN(5,1) = std::complex<double>(pars(286), pars(287));
   ZN(5,2) = std::complex<double>(pars(288), pars(289));
   ZN(5,3) = std::complex<double>(pars(290), pars(291));
   ZN(5,4) = std::complex<double>(pars(292), pars(293));
   ZN(5,5) = std::complex<double>(pars(294), pars(295));
   ZVL(0,0) = std::complex<double>(pars(296), pars(297));
   ZVL(0,1) = std::complex<double>(pars(298), pars(299));
   ZVL(0,2) = std::complex<double>(pars(300), pars(301));
   ZVL(1,0) = std::complex<double>(pars(302), pars(303));
   ZVL(1,1) = std::complex<double>(pars(304), pars(305));
   ZVL(1,2) = std::complex<double>(pars(306), pars(307));
   ZVL(2,0) = std::complex<double>(pars(308), pars(309));
   ZVL(2,1) = std::complex<double>(pars(310), pars(311));
   ZVL(2,2) = std::complex<double>(pars(312), pars(313));
   ZVR(0,0) = std::complex<double>(pars(314), pars(315));
   ZVR(0,1) = std::complex<double>(pars(316), pars(317));
   ZVR(0,2) = std::complex<double>(pars(318), pars(319));
   ZVR(1,0) = std::complex<double>(pars(320), pars(321));
   ZVR(1,1) = std::complex<double>(pars(322), pars(323));
   ZVR(1,2) = std::complex<double>(pars(324), pars(325));
   ZVR(2,0) = std::complex<double>(pars(326), pars(327));
   ZVR(2,1) = std::complex<double>(pars(328), pars(329));
   ZVR(2,2) = std::complex<double>(pars(330), pars(331));
   UM(0,0) = std::complex<double>(pars(332), pars(333));
   UM(0,1) = std::complex<double>(pars(334), pars(335));
   UM(1,0) = std::complex<double>(pars(336), pars(337));
   UM(1,1) = std::complex<double>(pars(338), pars(339));
   UP(0,0) = std::complex<double>(pars(340), pars(341));
   UP(0,1) = std::complex<double>(pars(342), pars(343));
   UP(1,0) = std::complex<double>(pars(344), pars(345));
   UP(1,1) = std::complex<double>(pars(346), pars(347));
   ZEL(0,0) = std::complex<double>(pars(348), pars(349));
   ZEL(0,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,2) = std::complex<double>(pars(352), pars(353));
   ZEL(1,0) = std::complex<double>(pars(354), pars(355));
   ZEL(1,1) = std::complex<double>(pars(356), pars(357));
   ZEL(1,2) = std::complex<double>(pars(358), pars(359));
   ZEL(2,0) = std::complex<double>(pars(360), pars(361));
   ZEL(2,1) = std::complex<double>(pars(362), pars(363));
   ZEL(2,2) = std::complex<double>(pars(364), pars(365));
   ZER(0,0) = std::complex<double>(pars(366), pars(367));
   ZER(0,1) = std::complex<double>(pars(368), pars(369));
   ZER(0,2) = std::complex<double>(pars(370), pars(371));
   ZER(1,0) = std::complex<double>(pars(372), pars(373));
   ZER(1,1) = std::complex<double>(pars(374), pars(375));
   ZER(1,2) = std::complex<double>(pars(376), pars(377));
   ZER(2,0) = std::complex<double>(pars(378), pars(379));
   ZER(2,1) = std::complex<double>(pars(380), pars(381));
   ZER(2,2) = std::complex<double>(pars(382), pars(383));
   ZDL(0,0) = std::complex<double>(pars(384), pars(385));
   ZDL(0,1) = std::complex<double>(pars(386), pars(387));
   ZDL(0,2) = std::complex<double>(pars(388), pars(389));
   ZDL(1,0) = std::complex<double>(pars(390), pars(391));
   ZDL(1,1) = std::complex<double>(pars(392), pars(393));
   ZDL(1,2) = std::complex<double>(pars(394), pars(395));
   ZDL(2,0) = std::complex<double>(pars(396), pars(397));
   ZDL(2,1) = std::complex<double>(pars(398), pars(399));
   ZDL(2,2) = std::complex<double>(pars(400), pars(401));
   ZDR(0,0) = std::complex<double>(pars(402), pars(403));
   ZDR(0,1) = std::complex<double>(pars(404), pars(405));
   ZDR(0,2) = std::complex<double>(pars(406), pars(407));
   ZDR(1,0) = std::complex<double>(pars(408), pars(409));
   ZDR(1,1) = std::complex<double>(pars(410), pars(411));
   ZDR(1,2) = std::complex<double>(pars(412), pars(413));
   ZDR(2,0) = std::complex<double>(pars(414), pars(415));
   ZDR(2,1) = std::complex<double>(pars(416), pars(417));
   ZDR(2,2) = std::complex<double>(pars(418), pars(419));
   ZUL(0,0) = std::complex<double>(pars(420), pars(421));
   ZUL(0,1) = std::complex<double>(pars(422), pars(423));
   ZUL(0,2) = std::complex<double>(pars(424), pars(425));
   ZUL(1,0) = std::complex<double>(pars(426), pars(427));
   ZUL(1,1) = std::complex<double>(pars(428), pars(429));
   ZUL(1,2) = std::complex<double>(pars(430), pars(431));
   ZUL(2,0) = std::complex<double>(pars(432), pars(433));
   ZUL(2,1) = std::complex<double>(pars(434), pars(435));
   ZUL(2,2) = std::complex<double>(pars(436), pars(437));
   ZUR(0,0) = std::complex<double>(pars(438), pars(439));
   ZUR(0,1) = std::complex<double>(pars(440), pars(441));
   ZUR(0,2) = std::complex<double>(pars(442), pars(443));
   ZUR(1,0) = std::complex<double>(pars(444), pars(445));
   ZUR(1,1) = std::complex<double>(pars(446), pars(447));
   ZUR(1,2) = std::complex<double>(pars(448), pars(449));
   ZUR(2,0) = std::complex<double>(pars(450), pars(451));
   ZUR(2,1) = std::complex<double>(pars(452), pars(453));
   ZUR(2,2) = std::complex<double>(pars(454), pars(455));
   ZZ(0,0) = pars(456);
   ZZ(0,1) = pars(457);
   ZZ(0,2) = pars(458);
   ZZ(1,0) = pars(459);
   ZZ(1,1) = pars(460);
   ZZ(1,2) = pars(461);
   ZZ(2,0) = pars(462);
   ZZ(2,1) = pars(463);
   ZZ(2,2) = pars(464);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(465);

   pars(58) = ZD(0,0);
   pars(59) = ZD(0,1);
   pars(60) = ZD(0,2);
   pars(61) = ZD(0,3);
   pars(62) = ZD(0,4);
   pars(63) = ZD(0,5);
   pars(64) = ZD(1,0);
   pars(65) = ZD(1,1);
   pars(66) = ZD(1,2);
   pars(67) = ZD(1,3);
   pars(68) = ZD(1,4);
   pars(69) = ZD(1,5);
   pars(70) = ZD(2,0);
   pars(71) = ZD(2,1);
   pars(72) = ZD(2,2);
   pars(73) = ZD(2,3);
   pars(74) = ZD(2,4);
   pars(75) = ZD(2,5);
   pars(76) = ZD(3,0);
   pars(77) = ZD(3,1);
   pars(78) = ZD(3,2);
   pars(79) = ZD(3,3);
   pars(80) = ZD(3,4);
   pars(81) = ZD(3,5);
   pars(82) = ZD(4,0);
   pars(83) = ZD(4,1);
   pars(84) = ZD(4,2);
   pars(85) = ZD(4,3);
   pars(86) = ZD(4,4);
   pars(87) = ZD(4,5);
   pars(88) = ZD(5,0);
   pars(89) = ZD(5,1);
   pars(90) = ZD(5,2);
   pars(91) = ZD(5,3);
   pars(92) = ZD(5,4);
   pars(93) = ZD(5,5);
   pars(94) = ZV(0,0);
   pars(95) = ZV(0,1);
   pars(96) = ZV(0,2);
   pars(97) = ZV(0,3);
   pars(98) = ZV(0,4);
   pars(99) = ZV(0,5);
   pars(100) = ZV(1,0);
   pars(101) = ZV(1,1);
   pars(102) = ZV(1,2);
   pars(103) = ZV(1,3);
   pars(104) = ZV(1,4);
   pars(105) = ZV(1,5);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZV(2,3);
   pars(110) = ZV(2,4);
   pars(111) = ZV(2,5);
   pars(112) = ZV(3,0);
   pars(113) = ZV(3,1);
   pars(114) = ZV(3,2);
   pars(115) = ZV(3,3);
   pars(116) = ZV(3,4);
   pars(117) = ZV(3,5);
   pars(118) = ZV(4,0);
   pars(119) = ZV(4,1);
   pars(120) = ZV(4,2);
   pars(121) = ZV(4,3);
   pars(122) = ZV(4,4);
   pars(123) = ZV(4,5);
   pars(124) = ZV(5,0);
   pars(125) = ZV(5,1);
   pars(126) = ZV(5,2);
   pars(127) = ZV(5,3);
   pars(128) = ZV(5,4);
   pars(129) = ZV(5,5);
   pars(130) = ZU(0,0);
   pars(131) = ZU(0,1);
   pars(132) = ZU(0,2);
   pars(133) = ZU(0,3);
   pars(134) = ZU(0,4);
   pars(135) = ZU(0,5);
   pars(136) = ZU(1,0);
   pars(137) = ZU(1,1);
   pars(138) = ZU(1,2);
   pars(139) = ZU(1,3);
   pars(140) = ZU(1,4);
   pars(141) = ZU(1,5);
   pars(142) = ZU(2,0);
   pars(143) = ZU(2,1);
   pars(144) = ZU(2,2);
   pars(145) = ZU(2,3);
   pars(146) = ZU(2,4);
   pars(147) = ZU(2,5);
   pars(148) = ZU(3,0);
   pars(149) = ZU(3,1);
   pars(150) = ZU(3,2);
   pars(151) = ZU(3,3);
   pars(152) = ZU(3,4);
   pars(153) = ZU(3,5);
   pars(154) = ZU(4,0);
   pars(155) = ZU(4,1);
   pars(156) = ZU(4,2);
   pars(157) = ZU(4,3);
   pars(158) = ZU(4,4);
   pars(159) = ZU(4,5);
   pars(160) = ZU(5,0);
   pars(161) = ZU(5,1);
   pars(162) = ZU(5,2);
   pars(163) = ZU(5,3);
   pars(164) = ZU(5,4);
   pars(165) = ZU(5,5);
   pars(166) = ZE(0,0);
   pars(167) = ZE(0,1);
   pars(168) = ZE(0,2);
   pars(169) = ZE(0,3);
   pars(170) = ZE(0,4);
   pars(171) = ZE(0,5);
   pars(172) = ZE(1,0);
   pars(173) = ZE(1,1);
   pars(174) = ZE(1,2);
   pars(175) = ZE(1,3);
   pars(176) = ZE(1,4);
   pars(177) = ZE(1,5);
   pars(178) = ZE(2,0);
   pars(179) = ZE(2,1);
   pars(180) = ZE(2,2);
   pars(181) = ZE(2,3);
   pars(182) = ZE(2,4);
   pars(183) = ZE(2,5);
   pars(184) = ZE(3,0);
   pars(185) = ZE(3,1);
   pars(186) = ZE(3,2);
   pars(187) = ZE(3,3);
   pars(188) = ZE(3,4);
   pars(189) = ZE(3,5);
   pars(190) = ZE(4,0);
   pars(191) = ZE(4,1);
   pars(192) = ZE(4,2);
   pars(193) = ZE(4,3);
   pars(194) = ZE(4,4);
   pars(195) = ZE(4,5);
   pars(196) = ZE(5,0);
   pars(197) = ZE(5,1);
   pars(198) = ZE(5,2);
   pars(199) = ZE(5,3);
   pars(200) = ZE(5,4);
   pars(201) = ZE(5,5);
   pars(202) = ZH(0,0);
   pars(203) = ZH(0,1);
   pars(204) = ZH(0,2);
   pars(205) = ZH(1,0);
   pars(206) = ZH(1,1);
   pars(207) = ZH(1,2);
   pars(208) = ZH(2,0);
   pars(209) = ZH(2,1);
   pars(210) = ZH(2,2);
   pars(211) = ZA(0,0);
   pars(212) = ZA(0,1);
   pars(213) = ZA(0,2);
   pars(214) = ZA(1,0);
   pars(215) = ZA(1,1);
   pars(216) = ZA(1,2);
   pars(217) = ZA(2,0);
   pars(218) = ZA(2,1);
   pars(219) = ZA(2,2);
   pars(220) = ZP(0,0);
   pars(221) = ZP(0,1);
   pars(222) = ZP(1,0);
   pars(223) = ZP(1,1);
   pars(224) = Re(ZN(0,0));
   pars(225) = Im(ZN(0,0));
   pars(226) = Re(ZN(0,1));
   pars(227) = Im(ZN(0,1));
   pars(228) = Re(ZN(0,2));
   pars(229) = Im(ZN(0,2));
   pars(230) = Re(ZN(0,3));
   pars(231) = Im(ZN(0,3));
   pars(232) = Re(ZN(0,4));
   pars(233) = Im(ZN(0,4));
   pars(234) = Re(ZN(0,5));
   pars(235) = Im(ZN(0,5));
   pars(236) = Re(ZN(1,0));
   pars(237) = Im(ZN(1,0));
   pars(238) = Re(ZN(1,1));
   pars(239) = Im(ZN(1,1));
   pars(240) = Re(ZN(1,2));
   pars(241) = Im(ZN(1,2));
   pars(242) = Re(ZN(1,3));
   pars(243) = Im(ZN(1,3));
   pars(244) = Re(ZN(1,4));
   pars(245) = Im(ZN(1,4));
   pars(246) = Re(ZN(1,5));
   pars(247) = Im(ZN(1,5));
   pars(248) = Re(ZN(2,0));
   pars(249) = Im(ZN(2,0));
   pars(250) = Re(ZN(2,1));
   pars(251) = Im(ZN(2,1));
   pars(252) = Re(ZN(2,2));
   pars(253) = Im(ZN(2,2));
   pars(254) = Re(ZN(2,3));
   pars(255) = Im(ZN(2,3));
   pars(256) = Re(ZN(2,4));
   pars(257) = Im(ZN(2,4));
   pars(258) = Re(ZN(2,5));
   pars(259) = Im(ZN(2,5));
   pars(260) = Re(ZN(3,0));
   pars(261) = Im(ZN(3,0));
   pars(262) = Re(ZN(3,1));
   pars(263) = Im(ZN(3,1));
   pars(264) = Re(ZN(3,2));
   pars(265) = Im(ZN(3,2));
   pars(266) = Re(ZN(3,3));
   pars(267) = Im(ZN(3,3));
   pars(268) = Re(ZN(3,4));
   pars(269) = Im(ZN(3,4));
   pars(270) = Re(ZN(3,5));
   pars(271) = Im(ZN(3,5));
   pars(272) = Re(ZN(4,0));
   pars(273) = Im(ZN(4,0));
   pars(274) = Re(ZN(4,1));
   pars(275) = Im(ZN(4,1));
   pars(276) = Re(ZN(4,2));
   pars(277) = Im(ZN(4,2));
   pars(278) = Re(ZN(4,3));
   pars(279) = Im(ZN(4,3));
   pars(280) = Re(ZN(4,4));
   pars(281) = Im(ZN(4,4));
   pars(282) = Re(ZN(4,5));
   pars(283) = Im(ZN(4,5));
   pars(284) = Re(ZN(5,0));
   pars(285) = Im(ZN(5,0));
   pars(286) = Re(ZN(5,1));
   pars(287) = Im(ZN(5,1));
   pars(288) = Re(ZN(5,2));
   pars(289) = Im(ZN(5,2));
   pars(290) = Re(ZN(5,3));
   pars(291) = Im(ZN(5,3));
   pars(292) = Re(ZN(5,4));
   pars(293) = Im(ZN(5,4));
   pars(294) = Re(ZN(5,5));
   pars(295) = Im(ZN(5,5));
   pars(296) = Re(ZVL(0,0));
   pars(297) = Im(ZVL(0,0));
   pars(298) = Re(ZVL(0,1));
   pars(299) = Im(ZVL(0,1));
   pars(300) = Re(ZVL(0,2));
   pars(301) = Im(ZVL(0,2));
   pars(302) = Re(ZVL(1,0));
   pars(303) = Im(ZVL(1,0));
   pars(304) = Re(ZVL(1,1));
   pars(305) = Im(ZVL(1,1));
   pars(306) = Re(ZVL(1,2));
   pars(307) = Im(ZVL(1,2));
   pars(308) = Re(ZVL(2,0));
   pars(309) = Im(ZVL(2,0));
   pars(310) = Re(ZVL(2,1));
   pars(311) = Im(ZVL(2,1));
   pars(312) = Re(ZVL(2,2));
   pars(313) = Im(ZVL(2,2));
   pars(314) = Re(ZVR(0,0));
   pars(315) = Im(ZVR(0,0));
   pars(316) = Re(ZVR(0,1));
   pars(317) = Im(ZVR(0,1));
   pars(318) = Re(ZVR(0,2));
   pars(319) = Im(ZVR(0,2));
   pars(320) = Re(ZVR(1,0));
   pars(321) = Im(ZVR(1,0));
   pars(322) = Re(ZVR(1,1));
   pars(323) = Im(ZVR(1,1));
   pars(324) = Re(ZVR(1,2));
   pars(325) = Im(ZVR(1,2));
   pars(326) = Re(ZVR(2,0));
   pars(327) = Im(ZVR(2,0));
   pars(328) = Re(ZVR(2,1));
   pars(329) = Im(ZVR(2,1));
   pars(330) = Re(ZVR(2,2));
   pars(331) = Im(ZVR(2,2));
   pars(332) = Re(UM(0,0));
   pars(333) = Im(UM(0,0));
   pars(334) = Re(UM(0,1));
   pars(335) = Im(UM(0,1));
   pars(336) = Re(UM(1,0));
   pars(337) = Im(UM(1,0));
   pars(338) = Re(UM(1,1));
   pars(339) = Im(UM(1,1));
   pars(340) = Re(UP(0,0));
   pars(341) = Im(UP(0,0));
   pars(342) = Re(UP(0,1));
   pars(343) = Im(UP(0,1));
   pars(344) = Re(UP(1,0));
   pars(345) = Im(UP(1,0));
   pars(346) = Re(UP(1,1));
   pars(347) = Im(UP(1,1));
   pars(348) = Re(ZEL(0,0));
   pars(349) = Im(ZEL(0,0));
   pars(350) = Re(ZEL(0,1));
   pars(351) = Im(ZEL(0,1));
   pars(352) = Re(ZEL(0,2));
   pars(353) = Im(ZEL(0,2));
   pars(354) = Re(ZEL(1,0));
   pars(355) = Im(ZEL(1,0));
   pars(356) = Re(ZEL(1,1));
   pars(357) = Im(ZEL(1,1));
   pars(358) = Re(ZEL(1,2));
   pars(359) = Im(ZEL(1,2));
   pars(360) = Re(ZEL(2,0));
   pars(361) = Im(ZEL(2,0));
   pars(362) = Re(ZEL(2,1));
   pars(363) = Im(ZEL(2,1));
   pars(364) = Re(ZEL(2,2));
   pars(365) = Im(ZEL(2,2));
   pars(366) = Re(ZER(0,0));
   pars(367) = Im(ZER(0,0));
   pars(368) = Re(ZER(0,1));
   pars(369) = Im(ZER(0,1));
   pars(370) = Re(ZER(0,2));
   pars(371) = Im(ZER(0,2));
   pars(372) = Re(ZER(1,0));
   pars(373) = Im(ZER(1,0));
   pars(374) = Re(ZER(1,1));
   pars(375) = Im(ZER(1,1));
   pars(376) = Re(ZER(1,2));
   pars(377) = Im(ZER(1,2));
   pars(378) = Re(ZER(2,0));
   pars(379) = Im(ZER(2,0));
   pars(380) = Re(ZER(2,1));
   pars(381) = Im(ZER(2,1));
   pars(382) = Re(ZER(2,2));
   pars(383) = Im(ZER(2,2));
   pars(384) = Re(ZDL(0,0));
   pars(385) = Im(ZDL(0,0));
   pars(386) = Re(ZDL(0,1));
   pars(387) = Im(ZDL(0,1));
   pars(388) = Re(ZDL(0,2));
   pars(389) = Im(ZDL(0,2));
   pars(390) = Re(ZDL(1,0));
   pars(391) = Im(ZDL(1,0));
   pars(392) = Re(ZDL(1,1));
   pars(393) = Im(ZDL(1,1));
   pars(394) = Re(ZDL(1,2));
   pars(395) = Im(ZDL(1,2));
   pars(396) = Re(ZDL(2,0));
   pars(397) = Im(ZDL(2,0));
   pars(398) = Re(ZDL(2,1));
   pars(399) = Im(ZDL(2,1));
   pars(400) = Re(ZDL(2,2));
   pars(401) = Im(ZDL(2,2));
   pars(402) = Re(ZDR(0,0));
   pars(403) = Im(ZDR(0,0));
   pars(404) = Re(ZDR(0,1));
   pars(405) = Im(ZDR(0,1));
   pars(406) = Re(ZDR(0,2));
   pars(407) = Im(ZDR(0,2));
   pars(408) = Re(ZDR(1,0));
   pars(409) = Im(ZDR(1,0));
   pars(410) = Re(ZDR(1,1));
   pars(411) = Im(ZDR(1,1));
   pars(412) = Re(ZDR(1,2));
   pars(413) = Im(ZDR(1,2));
   pars(414) = Re(ZDR(2,0));
   pars(415) = Im(ZDR(2,0));
   pars(416) = Re(ZDR(2,1));
   pars(417) = Im(ZDR(2,1));
   pars(418) = Re(ZDR(2,2));
   pars(419) = Im(ZDR(2,2));
   pars(420) = Re(ZUL(0,0));
   pars(421) = Im(ZUL(0,0));
   pars(422) = Re(ZUL(0,1));
   pars(423) = Im(ZUL(0,1));
   pars(424) = Re(ZUL(0,2));
   pars(425) = Im(ZUL(0,2));
   pars(426) = Re(ZUL(1,0));
   pars(427) = Im(ZUL(1,0));
   pars(428) = Re(ZUL(1,1));
   pars(429) = Im(ZUL(1,1));
   pars(430) = Re(ZUL(1,2));
   pars(431) = Im(ZUL(1,2));
   pars(432) = Re(ZUL(2,0));
   pars(433) = Im(ZUL(2,0));
   pars(434) = Re(ZUL(2,1));
   pars(435) = Im(ZUL(2,1));
   pars(436) = Re(ZUL(2,2));
   pars(437) = Im(ZUL(2,2));
   pars(438) = Re(ZUR(0,0));
   pars(439) = Im(ZUR(0,0));
   pars(440) = Re(ZUR(0,1));
   pars(441) = Im(ZUR(0,1));
   pars(442) = Re(ZUR(0,2));
   pars(443) = Im(ZUR(0,2));
   pars(444) = Re(ZUR(1,0));
   pars(445) = Im(ZUR(1,0));
   pars(446) = Re(ZUR(1,1));
   pars(447) = Im(ZUR(1,1));
   pars(448) = Re(ZUR(1,2));
   pars(449) = Im(ZUR(1,2));
   pars(450) = Re(ZUR(2,0));
   pars(451) = Im(ZUR(2,0));
   pars(452) = Re(ZUR(2,1));
   pars(453) = Im(ZUR(2,1));
   pars(454) = Re(ZUR(2,2));
   pars(455) = Im(ZUR(2,2));
   pars(456) = ZZ(0,0);
   pars(457) = ZZ(0,1);
   pars(458) = ZZ(0,2);
   pars(459) = ZZ(1,0);
   pars(460) = ZZ(1,1);
   pars(461) = ZZ(1,2);
   pars(462) = ZZ(2,0);
   pars(463) = ZZ(2,1);
   pars(464) = ZZ(2,2);


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
   Eigen::Array<double,2,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;
   MAh_goldstone(1) = MVZp;

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

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
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
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
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
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) - 0.5*vS*vu*Conj(
      Yd(0,2))*Lambdax;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) - 0.5*vS*vu*Conj(
      Yd(1,2))*Lambdax;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vS*vu*Conj(
      Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)) +
      AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) + Conj(
      Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) + Conj(
      Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) + Conj(
      Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sv() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qv = LOCALINPUT(Qv);

   Eigen::Matrix<double,6,6> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(0,0)) + AbsSqr(Yv(0,1)) + AbsSqr(Yv(0,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(1,0) + Conj(
      Yv(0,1))*Yv(1,1) + Conj(Yv(0,2))*Yv(1,2));
   mass_matrix_Sv(0,2) = ml2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(2,0) + Conj(
      Yv(0,1))*Yv(2,1) + Conj(Yv(0,2))*Yv(2,2));
   mass_matrix_Sv(0,3) = 0.7071067811865475*vu*Conj(TYv(0,0)) - 0.5*vd*vS*Conj(
      Yv(0,0))*Lambdax;
   mass_matrix_Sv(0,4) = 0.7071067811865475*vu*Conj(TYv(0,1)) - 0.5*vd*vS*Conj(
      Yv(0,1))*Lambdax;
   mass_matrix_Sv(0,5) = 0.7071067811865475*vu*Conj(TYv(0,2)) - 0.5*vd*vS*Conj(
      Yv(0,2))*Lambdax;
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(1,0)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(1,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(1,0))*Yv(2,0) + Conj(
      Yv(1,1))*Yv(2,1) + Conj(Yv(1,2))*Yv(2,2));
   mass_matrix_Sv(1,3) = 0.7071067811865475*vu*Conj(TYv(1,0)) - 0.5*vd*vS*Conj(
      Yv(1,0))*Lambdax;
   mass_matrix_Sv(1,4) = 0.7071067811865475*vu*Conj(TYv(1,1)) - 0.5*vd*vS*Conj(
      Yv(1,1))*Lambdax;
   mass_matrix_Sv(1,5) = 0.7071067811865475*vu*Conj(TYv(1,2)) - 0.5*vd*vS*Conj(
      Yv(1,2))*Lambdax;
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(2,0)) + AbsSqr(Yv(2,1)) + AbsSqr(Yv(2,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(2,3) = 0.7071067811865475*vu*Conj(TYv(2,0)) - 0.5*vd*vS*Conj(
      Yv(2,0))*Lambdax;
   mass_matrix_Sv(2,4) = 0.7071067811865475*vu*Conj(TYv(2,1)) - 0.5*vd*vS*Conj(
      Yv(2,1))*Lambdax;
   mass_matrix_Sv(2,5) = 0.7071067811865475*vu*Conj(TYv(2,2)) - 0.5*vd*vS*Conj(
      Yv(2,2))*Lambdax;
   mass_matrix_Sv(3,3) = mvR2(0,0) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,0)) + AbsSqr(Yv(1,0)) + AbsSqr(Yv(2,0)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(3,4) = mvR2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,1))*Yv(0,0) + Conj(
      Yv(1,1))*Yv(1,0) + Conj(Yv(2,1))*Yv(2,0));
   mass_matrix_Sv(3,5) = mvR2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,0) + Conj(
      Yv(1,2))*Yv(1,0) + Conj(Yv(2,2))*Yv(2,0));
   mass_matrix_Sv(4,4) = mvR2(1,1) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,1)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(2,1)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(4,5) = mvR2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,1) + Conj(
      Yv(1,2))*Yv(1,1) + Conj(Yv(2,2))*Yv(2,1));
   mass_matrix_Sv(5,5) = mvR2(2,2) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,2)) + AbsSqr(Yv(1,2)) + AbsSqr(Yv(2,2)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qu = LOCALINPUT(Qu);

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
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
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) - 0.5*vd*vS*Conj(
      Yu(0,1))*Lambdax;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vS*Conj(
      Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) - 0.5*vd*vS*Conj(
      Yu(2,1))*Lambdax;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) - 0.5*vd*vS*Conj(
      Yu(0,2))*Lambdax;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) - 0.5*vd*vS*Conj(
      Yu(1,2))*Lambdax;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vS*Conj(
      Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(0,
      1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) + Conj(
      Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) + Conj(
      Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(1,0)) + AbsSqr(Yu(1,
      1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) + Conj(
      Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(2,0)) + AbsSqr(Yu(2,
      1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qe = LOCALINPUT(Qe);

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
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
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
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
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) - 0.5*vS*vu*Conj(
      Ye(0,2))*Lambdax;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) - 0.5*vS*vu*Conj(
      Ye(1,2))*Lambdax;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vS*vu*Conj(
      Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)) +
      AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) + Conj(
      Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) + Conj(
      Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) + Conj(
      Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      1.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vS*Conj(
      TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) + QHd*QHu*vd*vu*Sqr(
      gp) - 0.35355339059327373*vS*TLambdax;
   mass_matrix_hh(0,2) = vd*vS*AbsSqr(Lambdax) - 0.35355339059327373*vu*Conj(
      TLambdax) + QHd*Qs*vd*vS*Sqr(gp) - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(
      Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.225*Sqr(g1)*Sqr(vu) +
      0.375*Sqr(g2)*Sqr(vu) + 1.5*Sqr(gp)*Sqr(QHu)*Sqr(vu);
   mass_matrix_hh(1,2) = vS*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*Conj(
      TLambdax) + QHu*Qs*vS*vu*Sqr(gp) - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*Sqr(gp)
      *Sqr(vd) + 1.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) +
      0.5*QHu*Qs*Sqr(gp)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) +
      0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu) + 0.3872983346207417*
      g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Cos(ThetaWp())) + Sqr(gp)*
      Sqr(QHd)*Sqr(vd)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Cos(ThetaWp()))
      *Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW()
      )*Sqr(vd)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHd)*Sqr(vd)*Sqr(Sin(ThetaWp(
      ))) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*
      Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vS*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(Cos(ThetaWp())) - 0.25*vd*vu*Sqr(
      g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 0.15*vd*vu*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) - 0.3872983346207417*g1*g2*vd*vu*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(
      Sin(ThetaWp())) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()
      )) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vS*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) + QHd*Qs*vd*vS*
      Sqr(gp)*Sqr(Cos(ThetaWp())) + QHd*Qs*vd*vS*Sqr(gp)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(
      Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu) + 0.5*Sqr(gp)*Sqr(QHu)*Sqr(vu) + 0.3872983346207417
      *g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Cos(ThetaWp())) + Sqr(gp)*
      Sqr(QHu)*Sqr(vu)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Cos(ThetaWp()))
      *Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW()
      )*Sqr(vu)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHu)*Sqr(vu)*Sqr(Sin(ThetaWp(
      ))) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*
      Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) + QHu*Qs*vS*vu*
      Sqr(gp)*Sqr(Cos(ThetaWp())) + QHu*Qs*vS*vu*Sqr(gp)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*Sqr(gp)
      *Sqr(vd) + 0.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) +
      0.5*QHu*Qs*Sqr(gp)*Sqr(vu) + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Cos(ThetaWp()))
      + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Sin(ThetaWp()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd)
      + 0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs
      *Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*
      QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475*vS*
      Conj(TLambdax);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHu*Qs*
      Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu) + 0.5*Sqr
      (gp)*Sqr(QHu)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Hpm, eigenvalue_error > precision * Abs(
      MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Chi() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   Eigen::Matrix<double,6,6> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassU;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = 0;
   mass_matrix_Chi(0,3) = gp*QHd*vd;
   mass_matrix_Chi(0,4) = gp*QHu*vu;
   mass_matrix_Chi(0,5) = gp*Qs*vS;
   mass_matrix_Chi(1,1) = MassB;
   mass_matrix_Chi(1,2) = 0;
   mass_matrix_Chi(1,3) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(1,4) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,5) = 0;
   mass_matrix_Chi(2,2) = MassWB;
   mass_matrix_Chi(2,3) = 0.5*g2*vd;
   mass_matrix_Chi(2,4) = -0.5*g2*vu;
   mass_matrix_Chi(2,5) = 0;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vS*Lambdax;
   mass_matrix_Chi(3,5) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(4,4) = 0;
   mass_matrix_Chi(4,5) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(5,5) = 0;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Chi, eigenvalue_error > precision * Abs(
      MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0.7071067811865475*vu*Yv(0,0);
   mass_matrix_Fv(0,1) = 0.7071067811865475*vu*Yv(0,1);
   mass_matrix_Fv(0,2) = 0.7071067811865475*vu*Yv(0,2);
   mass_matrix_Fv(1,0) = 0.7071067811865475*vu*Yv(1,0);
   mass_matrix_Fv(1,1) = 0.7071067811865475*vu*Yv(1,1);
   mass_matrix_Fv(1,2) = 0.7071067811865475*vu*Yv(1,2);
   mass_matrix_Fv(2,0) = 0.7071067811865475*vu*Yv(2,0);
   mass_matrix_Fv(2,1) = 0.7071067811865475*vu*Yv(2,1);
   mass_matrix_Fv(2,2) = 0.7071067811865475*vu*Yv(2,2);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   const auto mass_matrix_Fv(get_mass_matrix_Fv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fv, MFv, ZVL, ZVR, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Fv, eigenvalue_error > precision * Abs(
      MFv(0)));
#else
   fs_svd(mass_matrix_Fv, MFv, ZVL, ZVR);
#endif

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
   problems.flag_bad_mass(UMSSM_info::Cha, eigenvalue_error > precision * Abs(
      MCha(0)));
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
   problems.flag_bad_mass(UMSSM_info::Fe, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(UMSSM_info::Fd, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(UMSSM_info::Fu, eigenvalue_error > precision * Abs(
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
      problems.flag_running_tachyon(UMSSM_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_VPVZVZp() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   Eigen::Matrix<double,3,3> mass_matrix_VPVZVZp;

   mass_matrix_VPVZVZp(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZVZp(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZVZp(0,2) = -0.3872983346207417*g1*gp*QHd*Sqr(vd) +
      0.3872983346207417*g1*gp*QHu*Sqr(vu);
   mass_matrix_VPVZVZp(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_VPVZVZp(1,2) = 0.5*g2*gp*QHd*Sqr(vd) - 0.5*g2*gp*QHu*Sqr(vu);
   mass_matrix_VPVZVZp(2,2) = Sqr(gp)*Sqr(QHd)*Sqr(vd) + Sqr(gp)*Sqr(Qs)*Sqr(vS
      ) + Sqr(gp)*Sqr(QHu)*Sqr(vu);

   Symmetrize(mass_matrix_VPVZVZp);

   return mass_matrix_VPVZVZp;
}

void CLASSNAME::calculate_MVPVZVZp()
{
   const auto mass_matrix_VPVZVZp(get_mass_matrix_VPVZVZp());
   Eigen::Array<double,3,1> MVPVZVZp;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZVZp, MVPVZVZp, ZZ, eigenvalue_error
      );
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZVZp, MVPVZVZp, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZVZp = AbsSqrt(MVPVZVZp);

   MVP = 0.;
   MVZ = MVPVZVZp(1);
   MVZp = MVPVZVZp(2);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(mHd2*vd - 0.35355339059327373*vS*vu*Conj(TLambdax) + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.5*Cube(vd)*Sqr(gp)*Sqr(QHd) +
      0.5*vd*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*vd*Sqr(gp)*Sqr(vS) + 0.5*vd*
      AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu
      ) + 0.5*QHd*QHu*vd*Sqr(gp)*Sqr(vu) - 0.35355339059327373*vS*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   double result = Re(mHu2*vu - 0.35355339059327373*vd*vS*Conj(TLambdax) + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.5*Cube(vu)*Sqr(gp)*Sqr(QHu) +
      0.5*vu*AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*QHu*vu*Sqr(gp)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vS) +
      0.5*QHu*Qs*vu*Sqr(gp)*Sqr(vS) - 0.35355339059327373*vd*vS*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(ms2*vS - 0.35355339059327373*vd*vu*Conj(TLambdax) + 0.5*Cube
      (vS)*Sqr(gp)*Sqr(Qs) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*vS*Sqr(gp
      )*Sqr(vd) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vu) + 0.5*QHu*Qs*vS*Sqr(gp)*Sqr(vu) -
      0.35355339059327373*vd*vu*TLambdax);

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

double CLASSNAME::ThetaWp() const
{

   return ArcCos(Abs(ZZ(2,2)));
}

double CLASSNAME::VEV() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}



std::ostream& operator<<(std::ostream& ostr, const UMSSM_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
