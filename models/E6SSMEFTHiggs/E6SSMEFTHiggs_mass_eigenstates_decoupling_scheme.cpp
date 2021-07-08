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
 * @file E6SSMEFTHiggs_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the E6SSMEFTHiggs model class in the decoupling scheme
 *
 * Contains the definition of the E6SSMEFTHiggs model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#include "E6SSMEFTHiggs_mass_eigenstates_decoupling_scheme.hpp"
#include "E6SSMEFTHiggs_mass_eigenstates.hpp"
#include "E6SSMEFTHiggs_info.hpp"
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

#define CLASSNAME E6SSMEFTHiggs_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const E6SSMEFTHiggs_input_parameters& input_)
   : E6SSMEFTHiggs_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const E6SSMEFTHiggs_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<E6SSMEFTHiggs_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<E6SSMEFTHiggs_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(E6SSMEFTHiggs_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            E6SSMEFTHiggs_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(E6SSMEFTHiggs_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            E6SSMEFTHiggs_info::g2, new_g2, get_scale());
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
      MODEL->set_Yu((Diag((1.4142135623730951*upQuarksDRbar)/vu)).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Yd((Diag((1.4142135623730951*downQuarksDRbar)/vd)).real());

   }
   {
      auto& model = *this;
      auto MODEL = this;
      const auto vd = MODELPARAMETER(vd);
      MODEL->set_Ye((Diag((1.4142135623730951*downLeptonsDRbar)/vd)).real());

   }

   solve_ewsb_equations_tree_level();
   calculate_tree_level_mass_spectrum();
}

void CLASSNAME::fill_from(const E6SSMEFTHiggs_mass_eigenstates& model)
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
   MChaP = OTHER(MChaP);
   MSd = OTHER(MSd);
   ZD = OTHER(ZD);
   MSv = OTHER(MSv);
   ZV = OTHER(ZV);
   MSu = OTHER(MSu);
   ZU = OTHER(ZU);
   MSe = OTHER(MSe);
   ZE = OTHER(ZE);
   MSDX = OTHER(MSDX);
   ZDX = OTHER(ZDX);
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
   MFDX = OTHER(MFDX);
   ZDXL = OTHER(ZDXL);
   ZDXR = OTHER(ZDXR);
   MSHI0 = OTHER(MSHI0);
   UHI0 = OTHER(UHI0);
   MSHIp = OTHER(MSHIp);
   UHIp = OTHER(UHIp);
   MChaI = OTHER(MChaI);
   ZMI = OTHER(ZMI);
   ZPI = OTHER(ZPI);
   MChiI = OTHER(MChiI);
   ZNI = OTHER(ZNI);
   MSSI0 = OTHER(MSSI0);
   ZSSI = OTHER(ZSSI);
   MFSI = OTHER(MFSI);
   ZFSI = OTHER(ZFSI);
   MSHp0 = OTHER(MSHp0);
   UHp0 = OTHER(UHp0);
   MSHpp = OTHER(MSHpp);
   UHpp = OTHER(UHpp);
   MChiP = OTHER(MChiP);
   ZNp = OTHER(ZNp);
   MVWm = OTHER(MVWm);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   MVZp = OTHER(MVZp);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const E6SSMEFTHiggs_physical& CLASSNAME::get_physical() const
{
   return physical;
}

E6SSMEFTHiggs_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const E6SSMEFTHiggs_physical& physical_)
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

   mHd2 = Re((0.0125*(28.284271247461902*vs*vu*Conj(TLambdax) - 6*Cube(vd)*Sqr(g1)
      - 10*Cube(vd)*Sqr(g2) - 9*Cube(vd)*Sqr(gN) - 40*vd*AbsSqr(Lambdax)*Sqr(vs) +
      15*vd*Sqr(gN)*Sqr(vs) - 40*vd*AbsSqr(Lambdax)*Sqr(vu) + 6*vd*Sqr(g1)*Sqr(vu)
      + 10*vd*Sqr(g2)*Sqr(vu) - 6*vd*Sqr(gN)*Sqr(vu) + 28.284271247461902*vs*vu*
      TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vs*Conj(TLambdax) - 3*Cube(vu)*Sqr(g1)
      - 5*Cube(vu)*Sqr(g2) - 2*Cube(vu)*Sqr(gN) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) +
      3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 3*vu*Sqr(gN)*Sqr(vd) - 20*vu*
      AbsSqr(Lambdax)*Sqr(vs) + 5*vu*Sqr(gN)*Sqr(vs) + 14.142135623730951*vd*vs*
      TLambdax))/vu);
   ms2 = Re((0.0625*(5.656854249492381*vd*vu*Conj(TLambdax) - 5*Cube(vs)*Sqr(gN) -
      8*vs*AbsSqr(Lambdax)*Sqr(vd) + 3*vs*Sqr(gN)*Sqr(vd) - 8*vs*AbsSqr(Lambdax)*
      Sqr(vu) + 2*vs*Sqr(gN)*Sqr(vu) + 5.656854249492381*vd*vu*TLambdax))/vs);

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
           "E6SSMEFTHiggs\n"
           "========================================\n";
   E6SSMEFTHiggs_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MChaP = " << MChaP << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSDX = " << MSDX.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFDX = " << MFDX.transpose() << '\n';
   ostr << "MSHI0 = " << MSHI0.transpose() << '\n';
   ostr << "MSHIp = " << MSHIp.transpose() << '\n';
   ostr << "MChaI = " << MChaI.transpose() << '\n';
   ostr << "MChiI = " << MChiI.transpose() << '\n';
   ostr << "MSSI0 = " << MSSI0.transpose() << '\n';
   ostr << "MFSI = " << MFSI.transpose() << '\n';
   ostr << "MSHp0 = " << MSHp0.transpose() << '\n';
   ostr << "MSHpp = " << MSHpp.transpose() << '\n';
   ostr << "MChiP = " << MChiP.transpose() << '\n';
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
   ostr << "ZDX = " << ZDX << '\n';
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
   ostr << "ZDXL = " << ZDXL << '\n';
   ostr << "ZDXR = " << ZDXR << '\n';
   ostr << "UHI0 = " << UHI0 << '\n';
   ostr << "UHIp = " << UHIp << '\n';
   ostr << "ZMI = " << ZMI << '\n';
   ostr << "ZPI = " << ZPI << '\n';
   ostr << "ZNI = " << ZNI << '\n';
   ostr << "ZSSI = " << ZSSI << '\n';
   ostr << "ZFSI = " << ZFSI << '\n';
   ostr << "UHp0 = " << UHp0 << '\n';
   ostr << "UHpp = " << UHpp << '\n';
   ostr << "ZNp = " << ZNp << '\n';
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
   calculate_MChiP();
   calculate_MSHpp();
   calculate_MSHp0();
   calculate_MFSI();
   calculate_MSSI0();
   calculate_MChiI();
   calculate_MChaI();
   calculate_MSHIp();
   calculate_MSHI0();
   calculate_MFDX();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSDX();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
   calculate_MChaP();
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
   PHYSICAL(MChaP) = MChaP;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSv) = MSv;
   PHYSICAL(ZV) = ZV;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSDX) = MSDX;
   PHYSICAL(ZDX) = ZDX;
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
   PHYSICAL(MFDX) = MFDX;
   PHYSICAL(ZDXL) = ZDXL;
   PHYSICAL(ZDXR) = ZDXR;
   PHYSICAL(MSHI0) = MSHI0;
   PHYSICAL(UHI0) = UHI0;
   PHYSICAL(MSHIp) = MSHIp;
   PHYSICAL(UHIp) = UHIp;
   PHYSICAL(MChaI) = MChaI;
   PHYSICAL(ZMI) = ZMI;
   PHYSICAL(ZPI) = ZPI;
   PHYSICAL(MChiI) = MChiI;
   PHYSICAL(ZNI) = ZNI;
   PHYSICAL(MSSI0) = MSSI0;
   PHYSICAL(ZSSI) = ZSSI;
   PHYSICAL(MFSI) = MFSI;
   PHYSICAL(ZFSI) = ZFSI;
   PHYSICAL(MSHp0) = MSHp0;
   PHYSICAL(UHp0) = UHp0;
   PHYSICAL(MSHpp) = MSHpp;
   PHYSICAL(UHpp) = UHpp;
   PHYSICAL(MChiP) = MChiP;
   PHYSICAL(ZNp) = ZNp;
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Se); }
   if (PHYSICAL(MSDX).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SDX); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::Hpm); }
   if (PHYSICAL(MSHI0).tail<4>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SHI0); }
   if (PHYSICAL(MSHIp).tail<4>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SHIp); }
   if (PHYSICAL(MSSI0).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SSI0); }
   if (PHYSICAL(MSHp0).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SHp0); }
   if (PHYSICAL(MSHpp).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSMEFTHiggs_info::SHpp); }
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
   MChaP = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSDX = Eigen::Matrix<double,6,1>::Zero();
   ZDX = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
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
   MFDX = Eigen::Matrix<double,3,1>::Zero();
   ZDXL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDXR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MSHI0 = Eigen::Matrix<double,4,1>::Zero();
   UHI0 = Eigen::Matrix<double,4,4>::Zero();
   MSHIp = Eigen::Matrix<double,4,1>::Zero();
   UHIp = Eigen::Matrix<double,4,4>::Zero();
   MChaI = Eigen::Matrix<double,2,1>::Zero();
   ZMI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   ZPI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MChiI = Eigen::Matrix<double,4,1>::Zero();
   ZNI = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MSSI0 = Eigen::Matrix<double,2,1>::Zero();
   ZSSI = Eigen::Matrix<double,2,2>::Zero();
   MFSI = Eigen::Matrix<double,2,1>::Zero();
   ZFSI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MSHp0 = Eigen::Matrix<double,2,1>::Zero();
   UHp0 = Eigen::Matrix<double,2,2>::Zero();
   MSHpp = Eigen::Matrix<double,2,1>::Zero();
   UHpp = Eigen::Matrix<double,2,2>::Zero();
   MChiP = Eigen::Matrix<double,2,1>::Zero();
   ZNp = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;

   PhaseGlu = std::complex<double>(1.,0.);
   PhaseFHpup = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   E6SSMEFTHiggs_soft_parameters::clear();
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
   MChaP = pars(5);
   MSd(0) = pars(6);
   MSd(1) = pars(7);
   MSd(2) = pars(8);
   MSd(3) = pars(9);
   MSd(4) = pars(10);
   MSd(5) = pars(11);
   MSv(0) = pars(12);
   MSv(1) = pars(13);
   MSv(2) = pars(14);
   MSu(0) = pars(15);
   MSu(1) = pars(16);
   MSu(2) = pars(17);
   MSu(3) = pars(18);
   MSu(4) = pars(19);
   MSu(5) = pars(20);
   MSe(0) = pars(21);
   MSe(1) = pars(22);
   MSe(2) = pars(23);
   MSe(3) = pars(24);
   MSe(4) = pars(25);
   MSe(5) = pars(26);
   MSDX(0) = pars(27);
   MSDX(1) = pars(28);
   MSDX(2) = pars(29);
   MSDX(3) = pars(30);
   MSDX(4) = pars(31);
   MSDX(5) = pars(32);
   Mhh(0) = pars(33);
   Mhh(1) = pars(34);
   Mhh(2) = pars(35);
   MAh(0) = pars(36);
   MAh(1) = pars(37);
   MAh(2) = pars(38);
   MHpm(0) = pars(39);
   MHpm(1) = pars(40);
   MChi(0) = pars(41);
   MChi(1) = pars(42);
   MChi(2) = pars(43);
   MChi(3) = pars(44);
   MChi(4) = pars(45);
   MChi(5) = pars(46);
   MCha(0) = pars(47);
   MCha(1) = pars(48);
   MFe(0) = pars(49);
   MFe(1) = pars(50);
   MFe(2) = pars(51);
   MFd(0) = pars(52);
   MFd(1) = pars(53);
   MFd(2) = pars(54);
   MFu(0) = pars(55);
   MFu(1) = pars(56);
   MFu(2) = pars(57);
   MFDX(0) = pars(58);
   MFDX(1) = pars(59);
   MFDX(2) = pars(60);
   MSHI0(0) = pars(61);
   MSHI0(1) = pars(62);
   MSHI0(2) = pars(63);
   MSHI0(3) = pars(64);
   MSHIp(0) = pars(65);
   MSHIp(1) = pars(66);
   MSHIp(2) = pars(67);
   MSHIp(3) = pars(68);
   MChaI(0) = pars(69);
   MChaI(1) = pars(70);
   MChiI(0) = pars(71);
   MChiI(1) = pars(72);
   MChiI(2) = pars(73);
   MChiI(3) = pars(74);
   MSSI0(0) = pars(75);
   MSSI0(1) = pars(76);
   MFSI(0) = pars(77);
   MFSI(1) = pars(78);
   MSHp0(0) = pars(79);
   MSHp0(1) = pars(80);
   MSHpp(0) = pars(81);
   MSHpp(1) = pars(82);
   MChiP(0) = pars(83);
   MChiP(1) = pars(84);
   MVWm = pars(85);
   MVP = pars(86);
   MVZ = pars(87);
   MVZp = pars(88);

}

const E6SSMEFTHiggs_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

E6SSMEFTHiggs_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(89);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MChaP;
   pars(6) = MSd(0);
   pars(7) = MSd(1);
   pars(8) = MSd(2);
   pars(9) = MSd(3);
   pars(10) = MSd(4);
   pars(11) = MSd(5);
   pars(12) = MSv(0);
   pars(13) = MSv(1);
   pars(14) = MSv(2);
   pars(15) = MSu(0);
   pars(16) = MSu(1);
   pars(17) = MSu(2);
   pars(18) = MSu(3);
   pars(19) = MSu(4);
   pars(20) = MSu(5);
   pars(21) = MSe(0);
   pars(22) = MSe(1);
   pars(23) = MSe(2);
   pars(24) = MSe(3);
   pars(25) = MSe(4);
   pars(26) = MSe(5);
   pars(27) = MSDX(0);
   pars(28) = MSDX(1);
   pars(29) = MSDX(2);
   pars(30) = MSDX(3);
   pars(31) = MSDX(4);
   pars(32) = MSDX(5);
   pars(33) = Mhh(0);
   pars(34) = Mhh(1);
   pars(35) = Mhh(2);
   pars(36) = MAh(0);
   pars(37) = MAh(1);
   pars(38) = MAh(2);
   pars(39) = MHpm(0);
   pars(40) = MHpm(1);
   pars(41) = MChi(0);
   pars(42) = MChi(1);
   pars(43) = MChi(2);
   pars(44) = MChi(3);
   pars(45) = MChi(4);
   pars(46) = MChi(5);
   pars(47) = MCha(0);
   pars(48) = MCha(1);
   pars(49) = MFe(0);
   pars(50) = MFe(1);
   pars(51) = MFe(2);
   pars(52) = MFd(0);
   pars(53) = MFd(1);
   pars(54) = MFd(2);
   pars(55) = MFu(0);
   pars(56) = MFu(1);
   pars(57) = MFu(2);
   pars(58) = MFDX(0);
   pars(59) = MFDX(1);
   pars(60) = MFDX(2);
   pars(61) = MSHI0(0);
   pars(62) = MSHI0(1);
   pars(63) = MSHI0(2);
   pars(64) = MSHI0(3);
   pars(65) = MSHIp(0);
   pars(66) = MSHIp(1);
   pars(67) = MSHIp(2);
   pars(68) = MSHIp(3);
   pars(69) = MChaI(0);
   pars(70) = MChaI(1);
   pars(71) = MChiI(0);
   pars(72) = MChiI(1);
   pars(73) = MChiI(2);
   pars(74) = MChiI(3);
   pars(75) = MSSI0(0);
   pars(76) = MSSI0(1);
   pars(77) = MFSI(0);
   pars(78) = MFSI(1);
   pars(79) = MSHp0(0);
   pars(80) = MSHp0(1);
   pars(81) = MSHpp(0);
   pars(82) = MSHpp(1);
   pars(83) = MChiP(0);
   pars(84) = MChiP(1);
   pars(85) = MVWm;
   pars(86) = MVP;
   pars(87) = MVZ;
   pars(88) = MVZp;

   return pars;
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_tree_level_masses(pars);

   ZD(0,0) = pars(89);
   ZD(0,1) = pars(90);
   ZD(0,2) = pars(91);
   ZD(0,3) = pars(92);
   ZD(0,4) = pars(93);
   ZD(0,5) = pars(94);
   ZD(1,0) = pars(95);
   ZD(1,1) = pars(96);
   ZD(1,2) = pars(97);
   ZD(1,3) = pars(98);
   ZD(1,4) = pars(99);
   ZD(1,5) = pars(100);
   ZD(2,0) = pars(101);
   ZD(2,1) = pars(102);
   ZD(2,2) = pars(103);
   ZD(2,3) = pars(104);
   ZD(2,4) = pars(105);
   ZD(2,5) = pars(106);
   ZD(3,0) = pars(107);
   ZD(3,1) = pars(108);
   ZD(3,2) = pars(109);
   ZD(3,3) = pars(110);
   ZD(3,4) = pars(111);
   ZD(3,5) = pars(112);
   ZD(4,0) = pars(113);
   ZD(4,1) = pars(114);
   ZD(4,2) = pars(115);
   ZD(4,3) = pars(116);
   ZD(4,4) = pars(117);
   ZD(4,5) = pars(118);
   ZD(5,0) = pars(119);
   ZD(5,1) = pars(120);
   ZD(5,2) = pars(121);
   ZD(5,3) = pars(122);
   ZD(5,4) = pars(123);
   ZD(5,5) = pars(124);
   ZV(0,0) = pars(125);
   ZV(0,1) = pars(126);
   ZV(0,2) = pars(127);
   ZV(1,0) = pars(128);
   ZV(1,1) = pars(129);
   ZV(1,2) = pars(130);
   ZV(2,0) = pars(131);
   ZV(2,1) = pars(132);
   ZV(2,2) = pars(133);
   ZU(0,0) = pars(134);
   ZU(0,1) = pars(135);
   ZU(0,2) = pars(136);
   ZU(0,3) = pars(137);
   ZU(0,4) = pars(138);
   ZU(0,5) = pars(139);
   ZU(1,0) = pars(140);
   ZU(1,1) = pars(141);
   ZU(1,2) = pars(142);
   ZU(1,3) = pars(143);
   ZU(1,4) = pars(144);
   ZU(1,5) = pars(145);
   ZU(2,0) = pars(146);
   ZU(2,1) = pars(147);
   ZU(2,2) = pars(148);
   ZU(2,3) = pars(149);
   ZU(2,4) = pars(150);
   ZU(2,5) = pars(151);
   ZU(3,0) = pars(152);
   ZU(3,1) = pars(153);
   ZU(3,2) = pars(154);
   ZU(3,3) = pars(155);
   ZU(3,4) = pars(156);
   ZU(3,5) = pars(157);
   ZU(4,0) = pars(158);
   ZU(4,1) = pars(159);
   ZU(4,2) = pars(160);
   ZU(4,3) = pars(161);
   ZU(4,4) = pars(162);
   ZU(4,5) = pars(163);
   ZU(5,0) = pars(164);
   ZU(5,1) = pars(165);
   ZU(5,2) = pars(166);
   ZU(5,3) = pars(167);
   ZU(5,4) = pars(168);
   ZU(5,5) = pars(169);
   ZE(0,0) = pars(170);
   ZE(0,1) = pars(171);
   ZE(0,2) = pars(172);
   ZE(0,3) = pars(173);
   ZE(0,4) = pars(174);
   ZE(0,5) = pars(175);
   ZE(1,0) = pars(176);
   ZE(1,1) = pars(177);
   ZE(1,2) = pars(178);
   ZE(1,3) = pars(179);
   ZE(1,4) = pars(180);
   ZE(1,5) = pars(181);
   ZE(2,0) = pars(182);
   ZE(2,1) = pars(183);
   ZE(2,2) = pars(184);
   ZE(2,3) = pars(185);
   ZE(2,4) = pars(186);
   ZE(2,5) = pars(187);
   ZE(3,0) = pars(188);
   ZE(3,1) = pars(189);
   ZE(3,2) = pars(190);
   ZE(3,3) = pars(191);
   ZE(3,4) = pars(192);
   ZE(3,5) = pars(193);
   ZE(4,0) = pars(194);
   ZE(4,1) = pars(195);
   ZE(4,2) = pars(196);
   ZE(4,3) = pars(197);
   ZE(4,4) = pars(198);
   ZE(4,5) = pars(199);
   ZE(5,0) = pars(200);
   ZE(5,1) = pars(201);
   ZE(5,2) = pars(202);
   ZE(5,3) = pars(203);
   ZE(5,4) = pars(204);
   ZE(5,5) = pars(205);
   ZDX(0,0) = pars(206);
   ZDX(0,1) = pars(207);
   ZDX(0,2) = pars(208);
   ZDX(0,3) = pars(209);
   ZDX(0,4) = pars(210);
   ZDX(0,5) = pars(211);
   ZDX(1,0) = pars(212);
   ZDX(1,1) = pars(213);
   ZDX(1,2) = pars(214);
   ZDX(1,3) = pars(215);
   ZDX(1,4) = pars(216);
   ZDX(1,5) = pars(217);
   ZDX(2,0) = pars(218);
   ZDX(2,1) = pars(219);
   ZDX(2,2) = pars(220);
   ZDX(2,3) = pars(221);
   ZDX(2,4) = pars(222);
   ZDX(2,5) = pars(223);
   ZDX(3,0) = pars(224);
   ZDX(3,1) = pars(225);
   ZDX(3,2) = pars(226);
   ZDX(3,3) = pars(227);
   ZDX(3,4) = pars(228);
   ZDX(3,5) = pars(229);
   ZDX(4,0) = pars(230);
   ZDX(4,1) = pars(231);
   ZDX(4,2) = pars(232);
   ZDX(4,3) = pars(233);
   ZDX(4,4) = pars(234);
   ZDX(4,5) = pars(235);
   ZDX(5,0) = pars(236);
   ZDX(5,1) = pars(237);
   ZDX(5,2) = pars(238);
   ZDX(5,3) = pars(239);
   ZDX(5,4) = pars(240);
   ZDX(5,5) = pars(241);
   ZH(0,0) = pars(242);
   ZH(0,1) = pars(243);
   ZH(0,2) = pars(244);
   ZH(1,0) = pars(245);
   ZH(1,1) = pars(246);
   ZH(1,2) = pars(247);
   ZH(2,0) = pars(248);
   ZH(2,1) = pars(249);
   ZH(2,2) = pars(250);
   ZA(0,0) = pars(251);
   ZA(0,1) = pars(252);
   ZA(0,2) = pars(253);
   ZA(1,0) = pars(254);
   ZA(1,1) = pars(255);
   ZA(1,2) = pars(256);
   ZA(2,0) = pars(257);
   ZA(2,1) = pars(258);
   ZA(2,2) = pars(259);
   ZP(0,0) = pars(260);
   ZP(0,1) = pars(261);
   ZP(1,0) = pars(262);
   ZP(1,1) = pars(263);
   ZN(0,0) = std::complex<double>(pars(264), pars(265));
   ZN(0,1) = std::complex<double>(pars(266), pars(267));
   ZN(0,2) = std::complex<double>(pars(268), pars(269));
   ZN(0,3) = std::complex<double>(pars(270), pars(271));
   ZN(0,4) = std::complex<double>(pars(272), pars(273));
   ZN(0,5) = std::complex<double>(pars(274), pars(275));
   ZN(1,0) = std::complex<double>(pars(276), pars(277));
   ZN(1,1) = std::complex<double>(pars(278), pars(279));
   ZN(1,2) = std::complex<double>(pars(280), pars(281));
   ZN(1,3) = std::complex<double>(pars(282), pars(283));
   ZN(1,4) = std::complex<double>(pars(284), pars(285));
   ZN(1,5) = std::complex<double>(pars(286), pars(287));
   ZN(2,0) = std::complex<double>(pars(288), pars(289));
   ZN(2,1) = std::complex<double>(pars(290), pars(291));
   ZN(2,2) = std::complex<double>(pars(292), pars(293));
   ZN(2,3) = std::complex<double>(pars(294), pars(295));
   ZN(2,4) = std::complex<double>(pars(296), pars(297));
   ZN(2,5) = std::complex<double>(pars(298), pars(299));
   ZN(3,0) = std::complex<double>(pars(300), pars(301));
   ZN(3,1) = std::complex<double>(pars(302), pars(303));
   ZN(3,2) = std::complex<double>(pars(304), pars(305));
   ZN(3,3) = std::complex<double>(pars(306), pars(307));
   ZN(3,4) = std::complex<double>(pars(308), pars(309));
   ZN(3,5) = std::complex<double>(pars(310), pars(311));
   ZN(4,0) = std::complex<double>(pars(312), pars(313));
   ZN(4,1) = std::complex<double>(pars(314), pars(315));
   ZN(4,2) = std::complex<double>(pars(316), pars(317));
   ZN(4,3) = std::complex<double>(pars(318), pars(319));
   ZN(4,4) = std::complex<double>(pars(320), pars(321));
   ZN(4,5) = std::complex<double>(pars(322), pars(323));
   ZN(5,0) = std::complex<double>(pars(324), pars(325));
   ZN(5,1) = std::complex<double>(pars(326), pars(327));
   ZN(5,2) = std::complex<double>(pars(328), pars(329));
   ZN(5,3) = std::complex<double>(pars(330), pars(331));
   ZN(5,4) = std::complex<double>(pars(332), pars(333));
   ZN(5,5) = std::complex<double>(pars(334), pars(335));
   UM(0,0) = std::complex<double>(pars(336), pars(337));
   UM(0,1) = std::complex<double>(pars(338), pars(339));
   UM(1,0) = std::complex<double>(pars(340), pars(341));
   UM(1,1) = std::complex<double>(pars(342), pars(343));
   UP(0,0) = std::complex<double>(pars(344), pars(345));
   UP(0,1) = std::complex<double>(pars(346), pars(347));
   UP(1,0) = std::complex<double>(pars(348), pars(349));
   UP(1,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,0) = std::complex<double>(pars(352), pars(353));
   ZEL(0,1) = std::complex<double>(pars(354), pars(355));
   ZEL(0,2) = std::complex<double>(pars(356), pars(357));
   ZEL(1,0) = std::complex<double>(pars(358), pars(359));
   ZEL(1,1) = std::complex<double>(pars(360), pars(361));
   ZEL(1,2) = std::complex<double>(pars(362), pars(363));
   ZEL(2,0) = std::complex<double>(pars(364), pars(365));
   ZEL(2,1) = std::complex<double>(pars(366), pars(367));
   ZEL(2,2) = std::complex<double>(pars(368), pars(369));
   ZER(0,0) = std::complex<double>(pars(370), pars(371));
   ZER(0,1) = std::complex<double>(pars(372), pars(373));
   ZER(0,2) = std::complex<double>(pars(374), pars(375));
   ZER(1,0) = std::complex<double>(pars(376), pars(377));
   ZER(1,1) = std::complex<double>(pars(378), pars(379));
   ZER(1,2) = std::complex<double>(pars(380), pars(381));
   ZER(2,0) = std::complex<double>(pars(382), pars(383));
   ZER(2,1) = std::complex<double>(pars(384), pars(385));
   ZER(2,2) = std::complex<double>(pars(386), pars(387));
   ZDL(0,0) = std::complex<double>(pars(388), pars(389));
   ZDL(0,1) = std::complex<double>(pars(390), pars(391));
   ZDL(0,2) = std::complex<double>(pars(392), pars(393));
   ZDL(1,0) = std::complex<double>(pars(394), pars(395));
   ZDL(1,1) = std::complex<double>(pars(396), pars(397));
   ZDL(1,2) = std::complex<double>(pars(398), pars(399));
   ZDL(2,0) = std::complex<double>(pars(400), pars(401));
   ZDL(2,1) = std::complex<double>(pars(402), pars(403));
   ZDL(2,2) = std::complex<double>(pars(404), pars(405));
   ZDR(0,0) = std::complex<double>(pars(406), pars(407));
   ZDR(0,1) = std::complex<double>(pars(408), pars(409));
   ZDR(0,2) = std::complex<double>(pars(410), pars(411));
   ZDR(1,0) = std::complex<double>(pars(412), pars(413));
   ZDR(1,1) = std::complex<double>(pars(414), pars(415));
   ZDR(1,2) = std::complex<double>(pars(416), pars(417));
   ZDR(2,0) = std::complex<double>(pars(418), pars(419));
   ZDR(2,1) = std::complex<double>(pars(420), pars(421));
   ZDR(2,2) = std::complex<double>(pars(422), pars(423));
   ZUL(0,0) = std::complex<double>(pars(424), pars(425));
   ZUL(0,1) = std::complex<double>(pars(426), pars(427));
   ZUL(0,2) = std::complex<double>(pars(428), pars(429));
   ZUL(1,0) = std::complex<double>(pars(430), pars(431));
   ZUL(1,1) = std::complex<double>(pars(432), pars(433));
   ZUL(1,2) = std::complex<double>(pars(434), pars(435));
   ZUL(2,0) = std::complex<double>(pars(436), pars(437));
   ZUL(2,1) = std::complex<double>(pars(438), pars(439));
   ZUL(2,2) = std::complex<double>(pars(440), pars(441));
   ZUR(0,0) = std::complex<double>(pars(442), pars(443));
   ZUR(0,1) = std::complex<double>(pars(444), pars(445));
   ZUR(0,2) = std::complex<double>(pars(446), pars(447));
   ZUR(1,0) = std::complex<double>(pars(448), pars(449));
   ZUR(1,1) = std::complex<double>(pars(450), pars(451));
   ZUR(1,2) = std::complex<double>(pars(452), pars(453));
   ZUR(2,0) = std::complex<double>(pars(454), pars(455));
   ZUR(2,1) = std::complex<double>(pars(456), pars(457));
   ZUR(2,2) = std::complex<double>(pars(458), pars(459));
   ZDXL(0,0) = std::complex<double>(pars(460), pars(461));
   ZDXL(0,1) = std::complex<double>(pars(462), pars(463));
   ZDXL(0,2) = std::complex<double>(pars(464), pars(465));
   ZDXL(1,0) = std::complex<double>(pars(466), pars(467));
   ZDXL(1,1) = std::complex<double>(pars(468), pars(469));
   ZDXL(1,2) = std::complex<double>(pars(470), pars(471));
   ZDXL(2,0) = std::complex<double>(pars(472), pars(473));
   ZDXL(2,1) = std::complex<double>(pars(474), pars(475));
   ZDXL(2,2) = std::complex<double>(pars(476), pars(477));
   ZDXR(0,0) = std::complex<double>(pars(478), pars(479));
   ZDXR(0,1) = std::complex<double>(pars(480), pars(481));
   ZDXR(0,2) = std::complex<double>(pars(482), pars(483));
   ZDXR(1,0) = std::complex<double>(pars(484), pars(485));
   ZDXR(1,1) = std::complex<double>(pars(486), pars(487));
   ZDXR(1,2) = std::complex<double>(pars(488), pars(489));
   ZDXR(2,0) = std::complex<double>(pars(490), pars(491));
   ZDXR(2,1) = std::complex<double>(pars(492), pars(493));
   ZDXR(2,2) = std::complex<double>(pars(494), pars(495));
   UHI0(0,0) = pars(496);
   UHI0(0,1) = pars(497);
   UHI0(0,2) = pars(498);
   UHI0(0,3) = pars(499);
   UHI0(1,0) = pars(500);
   UHI0(1,1) = pars(501);
   UHI0(1,2) = pars(502);
   UHI0(1,3) = pars(503);
   UHI0(2,0) = pars(504);
   UHI0(2,1) = pars(505);
   UHI0(2,2) = pars(506);
   UHI0(2,3) = pars(507);
   UHI0(3,0) = pars(508);
   UHI0(3,1) = pars(509);
   UHI0(3,2) = pars(510);
   UHI0(3,3) = pars(511);
   UHIp(0,0) = pars(512);
   UHIp(0,1) = pars(513);
   UHIp(0,2) = pars(514);
   UHIp(0,3) = pars(515);
   UHIp(1,0) = pars(516);
   UHIp(1,1) = pars(517);
   UHIp(1,2) = pars(518);
   UHIp(1,3) = pars(519);
   UHIp(2,0) = pars(520);
   UHIp(2,1) = pars(521);
   UHIp(2,2) = pars(522);
   UHIp(2,3) = pars(523);
   UHIp(3,0) = pars(524);
   UHIp(3,1) = pars(525);
   UHIp(3,2) = pars(526);
   UHIp(3,3) = pars(527);
   ZMI(0,0) = std::complex<double>(pars(528), pars(529));
   ZMI(0,1) = std::complex<double>(pars(530), pars(531));
   ZMI(1,0) = std::complex<double>(pars(532), pars(533));
   ZMI(1,1) = std::complex<double>(pars(534), pars(535));
   ZPI(0,0) = std::complex<double>(pars(536), pars(537));
   ZPI(0,1) = std::complex<double>(pars(538), pars(539));
   ZPI(1,0) = std::complex<double>(pars(540), pars(541));
   ZPI(1,1) = std::complex<double>(pars(542), pars(543));
   ZNI(0,0) = std::complex<double>(pars(544), pars(545));
   ZNI(0,1) = std::complex<double>(pars(546), pars(547));
   ZNI(0,2) = std::complex<double>(pars(548), pars(549));
   ZNI(0,3) = std::complex<double>(pars(550), pars(551));
   ZNI(1,0) = std::complex<double>(pars(552), pars(553));
   ZNI(1,1) = std::complex<double>(pars(554), pars(555));
   ZNI(1,2) = std::complex<double>(pars(556), pars(557));
   ZNI(1,3) = std::complex<double>(pars(558), pars(559));
   ZNI(2,0) = std::complex<double>(pars(560), pars(561));
   ZNI(2,1) = std::complex<double>(pars(562), pars(563));
   ZNI(2,2) = std::complex<double>(pars(564), pars(565));
   ZNI(2,3) = std::complex<double>(pars(566), pars(567));
   ZNI(3,0) = std::complex<double>(pars(568), pars(569));
   ZNI(3,1) = std::complex<double>(pars(570), pars(571));
   ZNI(3,2) = std::complex<double>(pars(572), pars(573));
   ZNI(3,3) = std::complex<double>(pars(574), pars(575));
   ZSSI(0,0) = pars(576);
   ZSSI(0,1) = pars(577);
   ZSSI(1,0) = pars(578);
   ZSSI(1,1) = pars(579);
   ZFSI(0,0) = std::complex<double>(pars(580), pars(581));
   ZFSI(0,1) = std::complex<double>(pars(582), pars(583));
   ZFSI(1,0) = std::complex<double>(pars(584), pars(585));
   ZFSI(1,1) = std::complex<double>(pars(586), pars(587));
   UHp0(0,0) = pars(588);
   UHp0(0,1) = pars(589);
   UHp0(1,0) = pars(590);
   UHp0(1,1) = pars(591);
   UHpp(0,0) = pars(592);
   UHpp(0,1) = pars(593);
   UHpp(1,0) = pars(594);
   UHpp(1,1) = pars(595);
   ZNp(0,0) = std::complex<double>(pars(596), pars(597));
   ZNp(0,1) = std::complex<double>(pars(598), pars(599));
   ZNp(1,0) = std::complex<double>(pars(600), pars(601));
   ZNp(1,1) = std::complex<double>(pars(602), pars(603));
   ZZ(0,0) = pars(604);
   ZZ(0,1) = pars(605);
   ZZ(0,2) = pars(606);
   ZZ(1,0) = pars(607);
   ZZ(1,1) = pars(608);
   ZZ(1,2) = pars(609);
   ZZ(2,0) = pars(610);
   ZZ(2,1) = pars(611);
   ZZ(2,2) = pars(612);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(613);

   pars(89) = ZD(0,0);
   pars(90) = ZD(0,1);
   pars(91) = ZD(0,2);
   pars(92) = ZD(0,3);
   pars(93) = ZD(0,4);
   pars(94) = ZD(0,5);
   pars(95) = ZD(1,0);
   pars(96) = ZD(1,1);
   pars(97) = ZD(1,2);
   pars(98) = ZD(1,3);
   pars(99) = ZD(1,4);
   pars(100) = ZD(1,5);
   pars(101) = ZD(2,0);
   pars(102) = ZD(2,1);
   pars(103) = ZD(2,2);
   pars(104) = ZD(2,3);
   pars(105) = ZD(2,4);
   pars(106) = ZD(2,5);
   pars(107) = ZD(3,0);
   pars(108) = ZD(3,1);
   pars(109) = ZD(3,2);
   pars(110) = ZD(3,3);
   pars(111) = ZD(3,4);
   pars(112) = ZD(3,5);
   pars(113) = ZD(4,0);
   pars(114) = ZD(4,1);
   pars(115) = ZD(4,2);
   pars(116) = ZD(4,3);
   pars(117) = ZD(4,4);
   pars(118) = ZD(4,5);
   pars(119) = ZD(5,0);
   pars(120) = ZD(5,1);
   pars(121) = ZD(5,2);
   pars(122) = ZD(5,3);
   pars(123) = ZD(5,4);
   pars(124) = ZD(5,5);
   pars(125) = ZV(0,0);
   pars(126) = ZV(0,1);
   pars(127) = ZV(0,2);
   pars(128) = ZV(1,0);
   pars(129) = ZV(1,1);
   pars(130) = ZV(1,2);
   pars(131) = ZV(2,0);
   pars(132) = ZV(2,1);
   pars(133) = ZV(2,2);
   pars(134) = ZU(0,0);
   pars(135) = ZU(0,1);
   pars(136) = ZU(0,2);
   pars(137) = ZU(0,3);
   pars(138) = ZU(0,4);
   pars(139) = ZU(0,5);
   pars(140) = ZU(1,0);
   pars(141) = ZU(1,1);
   pars(142) = ZU(1,2);
   pars(143) = ZU(1,3);
   pars(144) = ZU(1,4);
   pars(145) = ZU(1,5);
   pars(146) = ZU(2,0);
   pars(147) = ZU(2,1);
   pars(148) = ZU(2,2);
   pars(149) = ZU(2,3);
   pars(150) = ZU(2,4);
   pars(151) = ZU(2,5);
   pars(152) = ZU(3,0);
   pars(153) = ZU(3,1);
   pars(154) = ZU(3,2);
   pars(155) = ZU(3,3);
   pars(156) = ZU(3,4);
   pars(157) = ZU(3,5);
   pars(158) = ZU(4,0);
   pars(159) = ZU(4,1);
   pars(160) = ZU(4,2);
   pars(161) = ZU(4,3);
   pars(162) = ZU(4,4);
   pars(163) = ZU(4,5);
   pars(164) = ZU(5,0);
   pars(165) = ZU(5,1);
   pars(166) = ZU(5,2);
   pars(167) = ZU(5,3);
   pars(168) = ZU(5,4);
   pars(169) = ZU(5,5);
   pars(170) = ZE(0,0);
   pars(171) = ZE(0,1);
   pars(172) = ZE(0,2);
   pars(173) = ZE(0,3);
   pars(174) = ZE(0,4);
   pars(175) = ZE(0,5);
   pars(176) = ZE(1,0);
   pars(177) = ZE(1,1);
   pars(178) = ZE(1,2);
   pars(179) = ZE(1,3);
   pars(180) = ZE(1,4);
   pars(181) = ZE(1,5);
   pars(182) = ZE(2,0);
   pars(183) = ZE(2,1);
   pars(184) = ZE(2,2);
   pars(185) = ZE(2,3);
   pars(186) = ZE(2,4);
   pars(187) = ZE(2,5);
   pars(188) = ZE(3,0);
   pars(189) = ZE(3,1);
   pars(190) = ZE(3,2);
   pars(191) = ZE(3,3);
   pars(192) = ZE(3,4);
   pars(193) = ZE(3,5);
   pars(194) = ZE(4,0);
   pars(195) = ZE(4,1);
   pars(196) = ZE(4,2);
   pars(197) = ZE(4,3);
   pars(198) = ZE(4,4);
   pars(199) = ZE(4,5);
   pars(200) = ZE(5,0);
   pars(201) = ZE(5,1);
   pars(202) = ZE(5,2);
   pars(203) = ZE(5,3);
   pars(204) = ZE(5,4);
   pars(205) = ZE(5,5);
   pars(206) = ZDX(0,0);
   pars(207) = ZDX(0,1);
   pars(208) = ZDX(0,2);
   pars(209) = ZDX(0,3);
   pars(210) = ZDX(0,4);
   pars(211) = ZDX(0,5);
   pars(212) = ZDX(1,0);
   pars(213) = ZDX(1,1);
   pars(214) = ZDX(1,2);
   pars(215) = ZDX(1,3);
   pars(216) = ZDX(1,4);
   pars(217) = ZDX(1,5);
   pars(218) = ZDX(2,0);
   pars(219) = ZDX(2,1);
   pars(220) = ZDX(2,2);
   pars(221) = ZDX(2,3);
   pars(222) = ZDX(2,4);
   pars(223) = ZDX(2,5);
   pars(224) = ZDX(3,0);
   pars(225) = ZDX(3,1);
   pars(226) = ZDX(3,2);
   pars(227) = ZDX(3,3);
   pars(228) = ZDX(3,4);
   pars(229) = ZDX(3,5);
   pars(230) = ZDX(4,0);
   pars(231) = ZDX(4,1);
   pars(232) = ZDX(4,2);
   pars(233) = ZDX(4,3);
   pars(234) = ZDX(4,4);
   pars(235) = ZDX(4,5);
   pars(236) = ZDX(5,0);
   pars(237) = ZDX(5,1);
   pars(238) = ZDX(5,2);
   pars(239) = ZDX(5,3);
   pars(240) = ZDX(5,4);
   pars(241) = ZDX(5,5);
   pars(242) = ZH(0,0);
   pars(243) = ZH(0,1);
   pars(244) = ZH(0,2);
   pars(245) = ZH(1,0);
   pars(246) = ZH(1,1);
   pars(247) = ZH(1,2);
   pars(248) = ZH(2,0);
   pars(249) = ZH(2,1);
   pars(250) = ZH(2,2);
   pars(251) = ZA(0,0);
   pars(252) = ZA(0,1);
   pars(253) = ZA(0,2);
   pars(254) = ZA(1,0);
   pars(255) = ZA(1,1);
   pars(256) = ZA(1,2);
   pars(257) = ZA(2,0);
   pars(258) = ZA(2,1);
   pars(259) = ZA(2,2);
   pars(260) = ZP(0,0);
   pars(261) = ZP(0,1);
   pars(262) = ZP(1,0);
   pars(263) = ZP(1,1);
   pars(264) = Re(ZN(0,0));
   pars(265) = Im(ZN(0,0));
   pars(266) = Re(ZN(0,1));
   pars(267) = Im(ZN(0,1));
   pars(268) = Re(ZN(0,2));
   pars(269) = Im(ZN(0,2));
   pars(270) = Re(ZN(0,3));
   pars(271) = Im(ZN(0,3));
   pars(272) = Re(ZN(0,4));
   pars(273) = Im(ZN(0,4));
   pars(274) = Re(ZN(0,5));
   pars(275) = Im(ZN(0,5));
   pars(276) = Re(ZN(1,0));
   pars(277) = Im(ZN(1,0));
   pars(278) = Re(ZN(1,1));
   pars(279) = Im(ZN(1,1));
   pars(280) = Re(ZN(1,2));
   pars(281) = Im(ZN(1,2));
   pars(282) = Re(ZN(1,3));
   pars(283) = Im(ZN(1,3));
   pars(284) = Re(ZN(1,4));
   pars(285) = Im(ZN(1,4));
   pars(286) = Re(ZN(1,5));
   pars(287) = Im(ZN(1,5));
   pars(288) = Re(ZN(2,0));
   pars(289) = Im(ZN(2,0));
   pars(290) = Re(ZN(2,1));
   pars(291) = Im(ZN(2,1));
   pars(292) = Re(ZN(2,2));
   pars(293) = Im(ZN(2,2));
   pars(294) = Re(ZN(2,3));
   pars(295) = Im(ZN(2,3));
   pars(296) = Re(ZN(2,4));
   pars(297) = Im(ZN(2,4));
   pars(298) = Re(ZN(2,5));
   pars(299) = Im(ZN(2,5));
   pars(300) = Re(ZN(3,0));
   pars(301) = Im(ZN(3,0));
   pars(302) = Re(ZN(3,1));
   pars(303) = Im(ZN(3,1));
   pars(304) = Re(ZN(3,2));
   pars(305) = Im(ZN(3,2));
   pars(306) = Re(ZN(3,3));
   pars(307) = Im(ZN(3,3));
   pars(308) = Re(ZN(3,4));
   pars(309) = Im(ZN(3,4));
   pars(310) = Re(ZN(3,5));
   pars(311) = Im(ZN(3,5));
   pars(312) = Re(ZN(4,0));
   pars(313) = Im(ZN(4,0));
   pars(314) = Re(ZN(4,1));
   pars(315) = Im(ZN(4,1));
   pars(316) = Re(ZN(4,2));
   pars(317) = Im(ZN(4,2));
   pars(318) = Re(ZN(4,3));
   pars(319) = Im(ZN(4,3));
   pars(320) = Re(ZN(4,4));
   pars(321) = Im(ZN(4,4));
   pars(322) = Re(ZN(4,5));
   pars(323) = Im(ZN(4,5));
   pars(324) = Re(ZN(5,0));
   pars(325) = Im(ZN(5,0));
   pars(326) = Re(ZN(5,1));
   pars(327) = Im(ZN(5,1));
   pars(328) = Re(ZN(5,2));
   pars(329) = Im(ZN(5,2));
   pars(330) = Re(ZN(5,3));
   pars(331) = Im(ZN(5,3));
   pars(332) = Re(ZN(5,4));
   pars(333) = Im(ZN(5,4));
   pars(334) = Re(ZN(5,5));
   pars(335) = Im(ZN(5,5));
   pars(336) = Re(UM(0,0));
   pars(337) = Im(UM(0,0));
   pars(338) = Re(UM(0,1));
   pars(339) = Im(UM(0,1));
   pars(340) = Re(UM(1,0));
   pars(341) = Im(UM(1,0));
   pars(342) = Re(UM(1,1));
   pars(343) = Im(UM(1,1));
   pars(344) = Re(UP(0,0));
   pars(345) = Im(UP(0,0));
   pars(346) = Re(UP(0,1));
   pars(347) = Im(UP(0,1));
   pars(348) = Re(UP(1,0));
   pars(349) = Im(UP(1,0));
   pars(350) = Re(UP(1,1));
   pars(351) = Im(UP(1,1));
   pars(352) = Re(ZEL(0,0));
   pars(353) = Im(ZEL(0,0));
   pars(354) = Re(ZEL(0,1));
   pars(355) = Im(ZEL(0,1));
   pars(356) = Re(ZEL(0,2));
   pars(357) = Im(ZEL(0,2));
   pars(358) = Re(ZEL(1,0));
   pars(359) = Im(ZEL(1,0));
   pars(360) = Re(ZEL(1,1));
   pars(361) = Im(ZEL(1,1));
   pars(362) = Re(ZEL(1,2));
   pars(363) = Im(ZEL(1,2));
   pars(364) = Re(ZEL(2,0));
   pars(365) = Im(ZEL(2,0));
   pars(366) = Re(ZEL(2,1));
   pars(367) = Im(ZEL(2,1));
   pars(368) = Re(ZEL(2,2));
   pars(369) = Im(ZEL(2,2));
   pars(370) = Re(ZER(0,0));
   pars(371) = Im(ZER(0,0));
   pars(372) = Re(ZER(0,1));
   pars(373) = Im(ZER(0,1));
   pars(374) = Re(ZER(0,2));
   pars(375) = Im(ZER(0,2));
   pars(376) = Re(ZER(1,0));
   pars(377) = Im(ZER(1,0));
   pars(378) = Re(ZER(1,1));
   pars(379) = Im(ZER(1,1));
   pars(380) = Re(ZER(1,2));
   pars(381) = Im(ZER(1,2));
   pars(382) = Re(ZER(2,0));
   pars(383) = Im(ZER(2,0));
   pars(384) = Re(ZER(2,1));
   pars(385) = Im(ZER(2,1));
   pars(386) = Re(ZER(2,2));
   pars(387) = Im(ZER(2,2));
   pars(388) = Re(ZDL(0,0));
   pars(389) = Im(ZDL(0,0));
   pars(390) = Re(ZDL(0,1));
   pars(391) = Im(ZDL(0,1));
   pars(392) = Re(ZDL(0,2));
   pars(393) = Im(ZDL(0,2));
   pars(394) = Re(ZDL(1,0));
   pars(395) = Im(ZDL(1,0));
   pars(396) = Re(ZDL(1,1));
   pars(397) = Im(ZDL(1,1));
   pars(398) = Re(ZDL(1,2));
   pars(399) = Im(ZDL(1,2));
   pars(400) = Re(ZDL(2,0));
   pars(401) = Im(ZDL(2,0));
   pars(402) = Re(ZDL(2,1));
   pars(403) = Im(ZDL(2,1));
   pars(404) = Re(ZDL(2,2));
   pars(405) = Im(ZDL(2,2));
   pars(406) = Re(ZDR(0,0));
   pars(407) = Im(ZDR(0,0));
   pars(408) = Re(ZDR(0,1));
   pars(409) = Im(ZDR(0,1));
   pars(410) = Re(ZDR(0,2));
   pars(411) = Im(ZDR(0,2));
   pars(412) = Re(ZDR(1,0));
   pars(413) = Im(ZDR(1,0));
   pars(414) = Re(ZDR(1,1));
   pars(415) = Im(ZDR(1,1));
   pars(416) = Re(ZDR(1,2));
   pars(417) = Im(ZDR(1,2));
   pars(418) = Re(ZDR(2,0));
   pars(419) = Im(ZDR(2,0));
   pars(420) = Re(ZDR(2,1));
   pars(421) = Im(ZDR(2,1));
   pars(422) = Re(ZDR(2,2));
   pars(423) = Im(ZDR(2,2));
   pars(424) = Re(ZUL(0,0));
   pars(425) = Im(ZUL(0,0));
   pars(426) = Re(ZUL(0,1));
   pars(427) = Im(ZUL(0,1));
   pars(428) = Re(ZUL(0,2));
   pars(429) = Im(ZUL(0,2));
   pars(430) = Re(ZUL(1,0));
   pars(431) = Im(ZUL(1,0));
   pars(432) = Re(ZUL(1,1));
   pars(433) = Im(ZUL(1,1));
   pars(434) = Re(ZUL(1,2));
   pars(435) = Im(ZUL(1,2));
   pars(436) = Re(ZUL(2,0));
   pars(437) = Im(ZUL(2,0));
   pars(438) = Re(ZUL(2,1));
   pars(439) = Im(ZUL(2,1));
   pars(440) = Re(ZUL(2,2));
   pars(441) = Im(ZUL(2,2));
   pars(442) = Re(ZUR(0,0));
   pars(443) = Im(ZUR(0,0));
   pars(444) = Re(ZUR(0,1));
   pars(445) = Im(ZUR(0,1));
   pars(446) = Re(ZUR(0,2));
   pars(447) = Im(ZUR(0,2));
   pars(448) = Re(ZUR(1,0));
   pars(449) = Im(ZUR(1,0));
   pars(450) = Re(ZUR(1,1));
   pars(451) = Im(ZUR(1,1));
   pars(452) = Re(ZUR(1,2));
   pars(453) = Im(ZUR(1,2));
   pars(454) = Re(ZUR(2,0));
   pars(455) = Im(ZUR(2,0));
   pars(456) = Re(ZUR(2,1));
   pars(457) = Im(ZUR(2,1));
   pars(458) = Re(ZUR(2,2));
   pars(459) = Im(ZUR(2,2));
   pars(460) = Re(ZDXL(0,0));
   pars(461) = Im(ZDXL(0,0));
   pars(462) = Re(ZDXL(0,1));
   pars(463) = Im(ZDXL(0,1));
   pars(464) = Re(ZDXL(0,2));
   pars(465) = Im(ZDXL(0,2));
   pars(466) = Re(ZDXL(1,0));
   pars(467) = Im(ZDXL(1,0));
   pars(468) = Re(ZDXL(1,1));
   pars(469) = Im(ZDXL(1,1));
   pars(470) = Re(ZDXL(1,2));
   pars(471) = Im(ZDXL(1,2));
   pars(472) = Re(ZDXL(2,0));
   pars(473) = Im(ZDXL(2,0));
   pars(474) = Re(ZDXL(2,1));
   pars(475) = Im(ZDXL(2,1));
   pars(476) = Re(ZDXL(2,2));
   pars(477) = Im(ZDXL(2,2));
   pars(478) = Re(ZDXR(0,0));
   pars(479) = Im(ZDXR(0,0));
   pars(480) = Re(ZDXR(0,1));
   pars(481) = Im(ZDXR(0,1));
   pars(482) = Re(ZDXR(0,2));
   pars(483) = Im(ZDXR(0,2));
   pars(484) = Re(ZDXR(1,0));
   pars(485) = Im(ZDXR(1,0));
   pars(486) = Re(ZDXR(1,1));
   pars(487) = Im(ZDXR(1,1));
   pars(488) = Re(ZDXR(1,2));
   pars(489) = Im(ZDXR(1,2));
   pars(490) = Re(ZDXR(2,0));
   pars(491) = Im(ZDXR(2,0));
   pars(492) = Re(ZDXR(2,1));
   pars(493) = Im(ZDXR(2,1));
   pars(494) = Re(ZDXR(2,2));
   pars(495) = Im(ZDXR(2,2));
   pars(496) = UHI0(0,0);
   pars(497) = UHI0(0,1);
   pars(498) = UHI0(0,2);
   pars(499) = UHI0(0,3);
   pars(500) = UHI0(1,0);
   pars(501) = UHI0(1,1);
   pars(502) = UHI0(1,2);
   pars(503) = UHI0(1,3);
   pars(504) = UHI0(2,0);
   pars(505) = UHI0(2,1);
   pars(506) = UHI0(2,2);
   pars(507) = UHI0(2,3);
   pars(508) = UHI0(3,0);
   pars(509) = UHI0(3,1);
   pars(510) = UHI0(3,2);
   pars(511) = UHI0(3,3);
   pars(512) = UHIp(0,0);
   pars(513) = UHIp(0,1);
   pars(514) = UHIp(0,2);
   pars(515) = UHIp(0,3);
   pars(516) = UHIp(1,0);
   pars(517) = UHIp(1,1);
   pars(518) = UHIp(1,2);
   pars(519) = UHIp(1,3);
   pars(520) = UHIp(2,0);
   pars(521) = UHIp(2,1);
   pars(522) = UHIp(2,2);
   pars(523) = UHIp(2,3);
   pars(524) = UHIp(3,0);
   pars(525) = UHIp(3,1);
   pars(526) = UHIp(3,2);
   pars(527) = UHIp(3,3);
   pars(528) = Re(ZMI(0,0));
   pars(529) = Im(ZMI(0,0));
   pars(530) = Re(ZMI(0,1));
   pars(531) = Im(ZMI(0,1));
   pars(532) = Re(ZMI(1,0));
   pars(533) = Im(ZMI(1,0));
   pars(534) = Re(ZMI(1,1));
   pars(535) = Im(ZMI(1,1));
   pars(536) = Re(ZPI(0,0));
   pars(537) = Im(ZPI(0,0));
   pars(538) = Re(ZPI(0,1));
   pars(539) = Im(ZPI(0,1));
   pars(540) = Re(ZPI(1,0));
   pars(541) = Im(ZPI(1,0));
   pars(542) = Re(ZPI(1,1));
   pars(543) = Im(ZPI(1,1));
   pars(544) = Re(ZNI(0,0));
   pars(545) = Im(ZNI(0,0));
   pars(546) = Re(ZNI(0,1));
   pars(547) = Im(ZNI(0,1));
   pars(548) = Re(ZNI(0,2));
   pars(549) = Im(ZNI(0,2));
   pars(550) = Re(ZNI(0,3));
   pars(551) = Im(ZNI(0,3));
   pars(552) = Re(ZNI(1,0));
   pars(553) = Im(ZNI(1,0));
   pars(554) = Re(ZNI(1,1));
   pars(555) = Im(ZNI(1,1));
   pars(556) = Re(ZNI(1,2));
   pars(557) = Im(ZNI(1,2));
   pars(558) = Re(ZNI(1,3));
   pars(559) = Im(ZNI(1,3));
   pars(560) = Re(ZNI(2,0));
   pars(561) = Im(ZNI(2,0));
   pars(562) = Re(ZNI(2,1));
   pars(563) = Im(ZNI(2,1));
   pars(564) = Re(ZNI(2,2));
   pars(565) = Im(ZNI(2,2));
   pars(566) = Re(ZNI(2,3));
   pars(567) = Im(ZNI(2,3));
   pars(568) = Re(ZNI(3,0));
   pars(569) = Im(ZNI(3,0));
   pars(570) = Re(ZNI(3,1));
   pars(571) = Im(ZNI(3,1));
   pars(572) = Re(ZNI(3,2));
   pars(573) = Im(ZNI(3,2));
   pars(574) = Re(ZNI(3,3));
   pars(575) = Im(ZNI(3,3));
   pars(576) = ZSSI(0,0);
   pars(577) = ZSSI(0,1);
   pars(578) = ZSSI(1,0);
   pars(579) = ZSSI(1,1);
   pars(580) = Re(ZFSI(0,0));
   pars(581) = Im(ZFSI(0,0));
   pars(582) = Re(ZFSI(0,1));
   pars(583) = Im(ZFSI(0,1));
   pars(584) = Re(ZFSI(1,0));
   pars(585) = Im(ZFSI(1,0));
   pars(586) = Re(ZFSI(1,1));
   pars(587) = Im(ZFSI(1,1));
   pars(588) = UHp0(0,0);
   pars(589) = UHp0(0,1);
   pars(590) = UHp0(1,0);
   pars(591) = UHp0(1,1);
   pars(592) = UHpp(0,0);
   pars(593) = UHpp(0,1);
   pars(594) = UHpp(1,0);
   pars(595) = UHpp(1,1);
   pars(596) = Re(ZNp(0,0));
   pars(597) = Im(ZNp(0,0));
   pars(598) = Re(ZNp(0,1));
   pars(599) = Im(ZNp(0,1));
   pars(600) = Re(ZNp(1,0));
   pars(601) = Im(ZNp(1,0));
   pars(602) = Re(ZNp(1,1));
   pars(603) = Im(ZNp(1,1));
   pars(604) = ZZ(0,0);
   pars(605) = ZZ(0,1);
   pars(606) = ZZ(0,2);
   pars(607) = ZZ(1,0);
   pars(608) = ZZ(1,1);
   pars(609) = ZZ(1,2);
   pars(610) = ZZ(2,0);
   pars(611) = ZZ(2,1);
   pars(612) = ZZ(2,2);


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

double CLASSNAME::get_mass_matrix_ChaP() const
{

   const double mass_matrix_ChaP = Re(MuPr);

   return mass_matrix_ChaP;
}

void CLASSNAME::calculate_MChaP()
{

   const auto mass_matrix_ChaP = get_mass_matrix_ChaP();
   MChaP = calculate_dirac_singlet_mass(mass_matrix_ChaP, PhaseFHpup);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1);
   mass_matrix_Sd(0,2) = mq2(0,2);
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vs*vu*Conj(
      Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0;
   mass_matrix_Sd(0,5) = 0;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2);
   mass_matrix_Sd(1,3) = 0;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vs*vu*Conj(
      Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(2,3) = 0;
   mass_matrix_Sd(2,4) = 0;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vs*vu*Conj(
      Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1);
   mass_matrix_Sd(3,5) = md2(0,2);
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2);
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Sd, eigenvalue_error > precision
      * Abs(MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Sv, eigenvalue_error > precision
      * Abs(MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(0,0
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1);
   mass_matrix_Su(0,2) = mq2(0,2);
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vs*Conj(
      Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0;
   mass_matrix_Su(0,5) = 0;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(1,1
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2);
   mass_matrix_Su(1,3) = 0;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vs*Conj(
      Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(2,2
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(2,3) = 0;
   mass_matrix_Su(2,4) = 0;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vs*Conj(
      Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(0,0))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1);
   mass_matrix_Su(3,5) = mu2(0,2);
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(1,1))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2);
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(2,2))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Su, eigenvalue_error > precision
      * Abs(MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1);
   mass_matrix_Se(0,2) = ml2(0,2);
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vs*vu*Conj(
      Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0;
   mass_matrix_Se(0,5) = 0;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2);
   mass_matrix_Se(1,3) = 0;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vs*vu*Conj(
      Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(2,3) = 0;
   mass_matrix_Se(2,4) = 0;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vs*vu*Conj(
      Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1);
   mass_matrix_Se(3,5) = me2(0,2);
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2);
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Se, eigenvalue_error > precision
      * Abs(MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_SDX() const
{

   Eigen::Matrix<double,6,6> mass_matrix_SDX;

   mass_matrix_SDX(0,0) = mDx2(0,0) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(0,0))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(0,1) = 0;
   mass_matrix_SDX(0,2) = 0;
   mass_matrix_SDX(0,3) = 0.7071067811865475*vs*Conj(TKappa(0,0)) - 0.5*vd*vu*
      Conj(Kappa(0,0))*Lambdax;
   mass_matrix_SDX(0,4) = 0;
   mass_matrix_SDX(0,5) = 0;
   mass_matrix_SDX(1,1) = mDx2(1,1) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(1,1))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(1,2) = 0;
   mass_matrix_SDX(1,3) = 0;
   mass_matrix_SDX(1,4) = 0.7071067811865475*vs*Conj(TKappa(1,1)) - 0.5*vd*vu*
      Conj(Kappa(1,1))*Lambdax;
   mass_matrix_SDX(1,5) = 0;
   mass_matrix_SDX(2,2) = mDx2(2,2) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(2,2))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(2,3) = 0;
   mass_matrix_SDX(2,4) = 0;
   mass_matrix_SDX(2,5) = 0.7071067811865475*vs*Conj(TKappa(2,2)) - 0.5*vd*vu*
      Conj(Kappa(2,2))*Lambdax;
   mass_matrix_SDX(3,3) = mDxbar2(0,0) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(0,0))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(3,4) = 0;
   mass_matrix_SDX(3,5) = 0;
   mass_matrix_SDX(4,4) = mDxbar2(1,1) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(1,1))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(4,5) = 0;
   mass_matrix_SDX(5,5) = mDxbar2(2,2) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(2,2))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SDX);

   return mass_matrix_SDX;
}

void CLASSNAME::calculate_MSDX()
{
   const auto mass_matrix_SDX(get_mass_matrix_SDX());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SDX, MSDX, ZDX, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SDX, eigenvalue_error > precision
       * Abs(MSDX(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SDX, MSDX, ZDX);
#endif
   normalize_to_interval(ZDX);


   if (MSDX.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SDX);
   }

   MSDX = AbsSqrt(MSDX);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      0.3375*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*Sqr
      (vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2
      )*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vs*Conj(
      TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) + 0.15*vd*vu*Sqr(gN)
      - 0.35355339059327373*vs*TLambdax;
   mass_matrix_hh(0,2) = vd*vs*AbsSqr(Lambdax) - 0.35355339059327373*vu*Conj(
      TLambdax) - 0.375*vd*vs*Sqr(gN) - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax
      )*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) + 0.225*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)
      *Sqr(vu) + 0.15*Sqr(gN)*Sqr(vu);
   mass_matrix_hh(1,2) = vs*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*Conj(
      TLambdax) - 0.25*vs*vu*Sqr(gN) - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*Sqr(gN)*Sqr
      (vd) + 0.9375*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.125*Sqr(
      gN)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::hh, eigenvalue_error > precision
      * Abs(Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) +
      0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*Sqr
      (vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2
      )*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu) + 0.3872983346207417*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(vd)*Sqr(Cos(ThetaWp())) + 0.225*Sqr(gN)*Sqr(vd)*Sqr(
      Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp
      ())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Sin(
      ThetaWp())) + 0.225*Sqr(gN)*Sqr(vd)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*
      Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr
      (Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vs*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + 0.15*vd*vu*Sqr(gN)*Sqr(Cos(ThetaWp())) - 0.25*vd*vu*Sqr(g2)
      *Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 0.15*vd*vu*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) - 0.3872983346207417*g1*g2*vd*vu*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 0.15*vd*vu*Sqr(gN)*Sqr(Sin(
      ThetaWp())) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) -
      0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vs*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) - 0.375*vd*vs*
      Sqr(gN)*Sqr(Cos(ThetaWp())) - 0.375*vd*vs*Sqr(gN)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax
      )*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)
      *Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu) + 0.3872983346207417*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())*Sqr(vu)*Sqr(Cos(ThetaWp())) + 0.1*Sqr(gN)*Sqr(vu)*Sqr(Cos(
      ThetaWp())) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Sin(
      ThetaWp())) + 0.1*Sqr(gN)*Sqr(vu)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*Sqr(
      vu)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin
      (ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) - 0.25*vs*vu*Sqr
      (gN)*Sqr(Cos(ThetaWp())) - 0.25*vs*vu*Sqr(gN)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*Sqr(gN)*Sqr
      (vd) + 0.3125*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.125*Sqr(
      gN)*Sqr(vu) + 0.625*Sqr(gN)*Sqr(vs)*Sqr(Cos(ThetaWp())) + 0.625*Sqr(gN)*
      Sqr(vs)*Sqr(Sin(ThetaWp()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Ah, eigenvalue_error > precision
      * Abs(MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd)
      + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.075*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475*vs*
      Conj(TLambdax);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.125*Sqr(gN)*Sqr
      (vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(
      vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Hpm, eigenvalue_error > precision
       * Abs(MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Chi() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(0,4) = 0;
   mass_matrix_Chi(0,5) = 0;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(1,5) = 0;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.7071067811865475*vs*Lambdax;
   mass_matrix_Chi(2,4) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(2,5) = -0.4743416490252569*gN*vd;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(3,5) = -0.31622776601683794*gN*vu;
   mass_matrix_Chi(4,4) = 0;
   mass_matrix_Chi(4,5) = 0.7905694150420949*gN*vs;
   mass_matrix_Chi(5,5) = MassBp;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Chi, eigenvalue_error > precision
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
   mass_matrix_Cha(1,1) = 0.7071067811865475*vs*Lambdax;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Cha, eigenvalue_error > precision
       * Abs(MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*vd*Ye(0,0);
   mass_matrix_Fe(0,1) = 0;
   mass_matrix_Fe(0,2) = 0;
   mass_matrix_Fe(1,1) = 0.7071067811865475*vd*Ye(1,1);
   mass_matrix_Fe(1,2) = 0;
   mass_matrix_Fe(2,2) = 0.7071067811865475*vd*Ye(2,2);

   Symmetrize(mass_matrix_Fe);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Fe, eigenvalue_error > precision
      * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*vd*Yd(0,0);
   mass_matrix_Fd(0,1) = 0;
   mass_matrix_Fd(0,2) = 0;
   mass_matrix_Fd(1,1) = 0.7071067811865475*vd*Yd(1,1);
   mass_matrix_Fd(1,2) = 0;
   mass_matrix_Fd(2,2) = 0.7071067811865475*vd*Yd(2,2);

   Symmetrize(mass_matrix_Fd);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Fd, eigenvalue_error > precision
      * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*vu*Yu(0,0);
   mass_matrix_Fu(0,1) = 0;
   mass_matrix_Fu(0,2) = 0;
   mass_matrix_Fu(1,1) = 0.7071067811865475*vu*Yu(1,1);
   mass_matrix_Fu(1,2) = 0;
   mass_matrix_Fu(2,2) = 0.7071067811865475*vu*Yu(2,2);

   Symmetrize(mass_matrix_Fu);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::Fu, eigenvalue_error > precision
      * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_FDX() const
{

   Eigen::Matrix<double,3,3> mass_matrix_FDX;

   mass_matrix_FDX(0,0) = 0.7071067811865475*vs*Kappa(0,0);
   mass_matrix_FDX(0,1) = 0;
   mass_matrix_FDX(0,2) = 0;
   mass_matrix_FDX(1,1) = 0.7071067811865475*vs*Kappa(1,1);
   mass_matrix_FDX(1,2) = 0;
   mass_matrix_FDX(2,2) = 0.7071067811865475*vs*Kappa(2,2);

   Symmetrize(mass_matrix_FDX);

   return mass_matrix_FDX;
}

void CLASSNAME::calculate_MFDX()
{
   const auto mass_matrix_FDX(get_mass_matrix_FDX());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_FDX, MFDX, ZDXL, ZDXR, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::FDX, eigenvalue_error > precision
       * Abs(MFDX(0)));
#else
   fs_svd(mass_matrix_FDX, MFDX, ZDXL, ZDXR);
#endif

}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_SHI0() const
{

   Eigen::Matrix<double,4,4> mass_matrix_SHI0;

   mass_matrix_SHI0(0,0) = mH1I2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(0,1) = 0;
   mass_matrix_SHI0(0,2) = -0.7071067811865475*vs*Conj(TLambda12(0,0)) + 0.5*vd
      *vu*Conj(Lambda12(0,0))*Lambdax;
   mass_matrix_SHI0(0,3) = 0;
   mass_matrix_SHI0(1,1) = mH1I2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(1,2) = 0;
   mass_matrix_SHI0(1,3) = -0.7071067811865475*vs*Conj(TLambda12(1,1)) + 0.5*vd
      *vu*Conj(Lambda12(1,1))*Lambdax;
   mass_matrix_SHI0(2,2) = mH2I2(0,0) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(2,3) = 0;
   mass_matrix_SHI0(3,3) = mH2I2(1,1) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHI0);

   return mass_matrix_SHI0;
}

void CLASSNAME::calculate_MSHI0()
{
   const auto mass_matrix_SHI0(get_mass_matrix_SHI0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHI0, MSHI0, UHI0, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SHI0, eigenvalue_error >
      precision * Abs(MSHI0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHI0, MSHI0, UHI0);
#endif
   normalize_to_interval(UHI0);


   if (MSHI0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SHI0);
   }

   MSHI0 = AbsSqrt(MSHI0);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_SHIp() const
{

   Eigen::Matrix<double,4,4> mass_matrix_SHIp;

   mass_matrix_SHIp(0,0) = mH1I2(0,0) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(0,1) = 0;
   mass_matrix_SHIp(0,2) = 0.7071067811865475*vs*Conj(TLambda12(0,0)) - 0.5*vd*
      vu*Conj(Lambda12(0,0))*Lambdax;
   mass_matrix_SHIp(0,3) = 0;
   mass_matrix_SHIp(1,1) = mH1I2(1,1) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(1,2) = 0;
   mass_matrix_SHIp(1,3) = 0.7071067811865475*vs*Conj(TLambda12(1,1)) - 0.5*vd*
      vu*Conj(Lambda12(1,1))*Lambdax;
   mass_matrix_SHIp(2,2) = mH2I2(0,0) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(2,3) = 0;
   mass_matrix_SHIp(3,3) = mH2I2(1,1) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHIp);

   return mass_matrix_SHIp;
}

void CLASSNAME::calculate_MSHIp()
{
   const auto mass_matrix_SHIp(get_mass_matrix_SHIp());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHIp, MSHIp, UHIp, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SHIp, eigenvalue_error >
      precision * Abs(MSHIp(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHIp, MSHIp, UHIp);
#endif
   normalize_to_interval(UHIp);


   if (MSHIp.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SHIp);
   }

   MSHIp = AbsSqrt(MSHIp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_ChaI() const
{

   Eigen::Matrix<double,2,2> mass_matrix_ChaI;

   mass_matrix_ChaI(0,0) = 0.7071067811865475*vs*Lambda12(0,0);
   mass_matrix_ChaI(0,1) = 0;
   mass_matrix_ChaI(1,1) = 0.7071067811865475*vs*Lambda12(1,1);

   Symmetrize(mass_matrix_ChaI);

   return mass_matrix_ChaI;
}

void CLASSNAME::calculate_MChaI()
{
   const auto mass_matrix_ChaI(get_mass_matrix_ChaI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_ChaI, MChaI, ZMI, ZPI, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::ChaI, eigenvalue_error >
      precision * Abs(MChaI(0)));
#else
   fs_svd(mass_matrix_ChaI, MChaI, ZMI, ZPI);
#endif

}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_ChiI() const
{

   Eigen::Matrix<double,4,4> mass_matrix_ChiI;

   mass_matrix_ChiI(0,0) = 0;
   mass_matrix_ChiI(0,1) = 0;
   mass_matrix_ChiI(0,2) = -0.7071067811865475*vs*Lambda12(0,0);
   mass_matrix_ChiI(0,3) = 0;
   mass_matrix_ChiI(1,1) = 0;
   mass_matrix_ChiI(1,2) = 0;
   mass_matrix_ChiI(1,3) = -0.7071067811865475*vs*Lambda12(1,1);
   mass_matrix_ChiI(2,2) = 0;
   mass_matrix_ChiI(2,3) = 0;
   mass_matrix_ChiI(3,3) = 0;

   Symmetrize(mass_matrix_ChiI);

   return mass_matrix_ChiI;
}

void CLASSNAME::calculate_MChiI()
{
   const auto mass_matrix_ChiI(get_mass_matrix_ChiI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_ChiI, MChiI, ZNI, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::ChiI, eigenvalue_error >
      precision * Abs(MChiI(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_ChiI, MChiI, ZNI);
#endif
   normalize_to_interval(ZNI);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SSI0() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SSI0;

   mass_matrix_SSI0(0,0) = msI2(0,0) - 0.1875*Sqr(gN)*Sqr(vd) + 0.3125*Sqr(gN)*
      Sqr(vs) - 0.125*Sqr(gN)*Sqr(vu);
   mass_matrix_SSI0(0,1) = 0;
   mass_matrix_SSI0(1,1) = msI2(1,1) - 0.1875*Sqr(gN)*Sqr(vd) + 0.3125*Sqr(gN)*
      Sqr(vs) - 0.125*Sqr(gN)*Sqr(vu);

   Symmetrize(mass_matrix_SSI0);

   return mass_matrix_SSI0;
}

void CLASSNAME::calculate_MSSI0()
{
   const auto mass_matrix_SSI0(get_mass_matrix_SSI0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SSI0, MSSI0, ZSSI, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SSI0, eigenvalue_error >
      precision * Abs(MSSI0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SSI0, MSSI0, ZSSI);
#endif
   normalize_to_interval(ZSSI);


   if (MSSI0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SSI0);
   }

   MSSI0 = AbsSqrt(MSSI0);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_FSI() const
{

   Eigen::Matrix<double,2,2> mass_matrix_FSI;

   mass_matrix_FSI(0,0) = 0;
   mass_matrix_FSI(0,1) = 0;
   mass_matrix_FSI(1,1) = 0;

   Symmetrize(mass_matrix_FSI);

   return mass_matrix_FSI;
}

void CLASSNAME::calculate_MFSI()
{
   const auto mass_matrix_FSI(get_mass_matrix_FSI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_FSI, MFSI, ZFSI, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::FSI, eigenvalue_error > precision
       * Abs(MFSI(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_FSI, MFSI, ZFSI);
#endif
   normalize_to_interval(ZFSI);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SHp0() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SHp0;

   mass_matrix_SHp0(0,0) = mHp2 + AbsSqr(MuPr) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*
      Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*
      Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHp0(0,1) = -Conj(BMuPr);
   mass_matrix_SHp0(1,1) = mHpbar2 + AbsSqr(MuPr) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) - 0.125*Sqr(gN)*Sqr(vs) +
      0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHp0);

   return mass_matrix_SHp0;
}

void CLASSNAME::calculate_MSHp0()
{
   const auto mass_matrix_SHp0(get_mass_matrix_SHp0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHp0, MSHp0, UHp0, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SHp0, eigenvalue_error >
      precision * Abs(MSHp0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHp0, MSHp0, UHp0);
#endif
   normalize_to_interval(UHp0);


   if (MSHp0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SHp0);
   }

   MSHp0 = AbsSqrt(MSHp0);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SHpp() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SHpp;

   mass_matrix_SHpp(0,0) = mHp2 + AbsSqr(MuPr) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*
      Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*
      Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHpp(0,1) = Conj(BMuPr);
   mass_matrix_SHpp(1,1) = mHpbar2 + AbsSqr(MuPr) - 0.075*Sqr(g1)*Sqr(vd) +
      0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) - 0.125*Sqr(gN)*Sqr(vs) +
      0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHpp);

   return mass_matrix_SHpp;
}

void CLASSNAME::calculate_MSHpp()
{
   const auto mass_matrix_SHpp(get_mass_matrix_SHpp());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHpp, MSHpp, UHpp, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::SHpp, eigenvalue_error >
      precision * Abs(MSHpp(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHpp, MSHpp, UHpp);
#endif
   normalize_to_interval(UHpp);


   if (MSHpp.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::SHpp);
   }

   MSHpp = AbsSqrt(MSHpp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_ChiP() const
{

   Eigen::Matrix<double,2,2> mass_matrix_ChiP;

   mass_matrix_ChiP(0,0) = 0;
   mass_matrix_ChiP(0,1) = -MuPr;
   mass_matrix_ChiP(1,1) = 0;

   Symmetrize(mass_matrix_ChiP);

   return mass_matrix_ChiP;
}

void CLASSNAME::calculate_MChiP()
{
   const auto mass_matrix_ChiP(get_mass_matrix_ChiP());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_ChiP, MChiP, ZNp, eigenvalue_error);
   problems.flag_bad_mass(E6SSMEFTHiggs_info::ChiP, eigenvalue_error >
      precision * Abs(MChiP(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_ChiP, MChiP, ZNp);
#endif
   normalize_to_interval(ZNp);

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
      problems.flag_running_tachyon(E6SSMEFTHiggs_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_VPVZVZp() const
{

   Eigen::Matrix<double,3,3> mass_matrix_VPVZVZp;

   mass_matrix_VPVZVZp(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZVZp(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZVZp(0,2) = 0.18371173070873834*g1*gN*Sqr(vd) -
      0.1224744871391589*g1*gN*Sqr(vu);
   mass_matrix_VPVZVZp(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_VPVZVZp(1,2) = -0.23717082451262844*g2*gN*Sqr(vd) +
      0.15811388300841897*g2*gN*Sqr(vu);
   mass_matrix_VPVZVZp(2,2) = 0.225*Sqr(gN)*Sqr(vd) + 0.625*Sqr(gN)*Sqr(vs) +
      0.1*Sqr(gN)*Sqr(vu);

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
   
   double result = Re(mHd2*vd - 0.35355339059327373*vs*vu*Conj(TLambdax) + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.1125*Cube(vd)*Sqr(gN) + 0.5*vd
      *AbsSqr(Lambdax)*Sqr(vs) - 0.1875*vd*Sqr(gN)*Sqr(vs) + 0.5*vd*AbsSqr(Lambdax
      )*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu) + 0.075*vd*
      Sqr(gN)*Sqr(vu) - 0.35355339059327373*vs*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu - 0.35355339059327373*vd*vs*Conj(TLambdax) + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.05*Cube(vu)*Sqr(gN) + 0.5*vu*
      AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd
      ) + 0.075*vu*Sqr(gN)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vs) - 0.125*vu*Sqr
      (gN)*Sqr(vs) - 0.35355339059327373*vd*vs*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(ms2*vs - 0.35355339059327373*vd*vu*Conj(TLambdax) + 0.3125*
      Cube(vs)*Sqr(gN) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*vs*Sqr(gN)*Sqr(vd
      ) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vu) - 0.125*vs*Sqr(gN)*Sqr(vu) -
      0.35355339059327373*vd*vu*TLambdax);

   return result;
}



double CLASSNAME::v() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}

double CLASSNAME::Betax() const
{

   return -ArcSin(ZA(0,1));
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



std::ostream& operator<<(std::ostream& ostr, const E6SSMEFTHiggs_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
