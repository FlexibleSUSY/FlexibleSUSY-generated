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
 * @file CMSSMNoFV_mass_eigenstates_decoupling_scheme.cpp
 * @brief implementation of the CMSSMNoFV model class in the decoupling scheme
 *
 * Contains the definition of the CMSSMNoFV model class methods
 * which solve EWSB and calculate masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#include "CMSSMNoFV_mass_eigenstates_decoupling_scheme.hpp"
#include "CMSSMNoFV_mass_eigenstates.hpp"
#include "CMSSMNoFV_info.hpp"
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

#define CLASSNAME CMSSMNoFV_mass_eigenstates_decoupling_scheme

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

CLASSNAME::CLASSNAME(const CMSSMNoFV_input_parameters& input_)
   : CMSSMNoFV_soft_parameters(input_)
{
}

CLASSNAME::CLASSNAME(const CMSSMNoFV_mass_eigenstates& model)
{
   fill_from(model);
}

std::unique_ptr<CMSSMNoFV_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<CMSSMNoFV_mass_eigenstates_decoupling_scheme>(*this);
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
         model->get_problems().unflag_non_perturbative_parameter(CMSSMNoFV_info::g1);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            CMSSMNoFV_info::g1, new_g1, get_scale());
         new_g1 = Electroweak_constants::g1;
      }

      if (IsFinite(new_g2)) {
         model->get_problems().unflag_non_perturbative_parameter(CMSSMNoFV_info::g2);
      } else {
         model->get_problems().flag_non_perturbative_parameter(
            CMSSMNoFV_info::g2, new_g2, get_scale());
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
      MODEL->set_Yu(((1.4142135623730951*upQuarksDRbar)/vu).real());

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

void CLASSNAME::fill_from(const CMSSMNoFV_mass_eigenstates& model)
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
   MFd = OTHER(MFd);
   MFs = OTHER(MFs);
   MFb = OTHER(MFb);
   MFu = OTHER(MFu);
   MFc = OTHER(MFc);
   MFt = OTHER(MFt);
   MFve = OTHER(MFve);
   MFvm = OTHER(MFvm);
   MFvt = OTHER(MFvt);
   MFe = OTHER(MFe);
   MFm = OTHER(MFm);
   MFtau = OTHER(MFtau);
   MSveL = OTHER(MSveL);
   MSvmL = OTHER(MSvmL);
   MSvtL = OTHER(MSvtL);
   MSd = OTHER(MSd);
   ZD = OTHER(ZD);
   MSu = OTHER(MSu);
   ZU = OTHER(ZU);
   MSe = OTHER(MSe);
   ZE = OTHER(ZE);
   MSm = OTHER(MSm);
   ZM = OTHER(ZM);
   MStau = OTHER(MStau);
   ZTau = OTHER(ZTau);
   MSs = OTHER(MSs);
   ZS = OTHER(ZS);
   MSc = OTHER(MSc);
   ZC = OTHER(ZC);
   MSb = OTHER(MSb);
   ZB = OTHER(ZB);
   MSt = OTHER(MSt);
   ZT = OTHER(ZT);
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
   MVWm = OTHER(MVWm);
   MVP = OTHER(MVP);
   MVZ = OTHER(MVZ);
   ZZ = OTHER(ZZ);

#undef OTHER
}

const CMSSMNoFV_physical& CLASSNAME::get_physical() const
{
   return physical;
}

CMSSMNoFV_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const CMSSMNoFV_physical& physical_)
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
           "CMSSMNoFV\n"
           "========================================\n";
   CMSSMNoFV_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFd = " << MFd << '\n';
   ostr << "MFs = " << MFs << '\n';
   ostr << "MFb = " << MFb << '\n';
   ostr << "MFu = " << MFu << '\n';
   ostr << "MFc = " << MFc << '\n';
   ostr << "MFt = " << MFt << '\n';
   ostr << "MFve = " << MFve << '\n';
   ostr << "MFvm = " << MFvm << '\n';
   ostr << "MFvt = " << MFvt << '\n';
   ostr << "MFe = " << MFe << '\n';
   ostr << "MFm = " << MFm << '\n';
   ostr << "MFtau = " << MFtau << '\n';
   ostr << "MSveL = " << MSveL << '\n';
   ostr << "MSvmL = " << MSvmL << '\n';
   ostr << "MSvtL = " << MSvtL << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSm = " << MSm.transpose() << '\n';
   ostr << "MStau = " << MStau.transpose() << '\n';
   ostr << "MSs = " << MSs.transpose() << '\n';
   ostr << "MSc = " << MSc.transpose() << '\n';
   ostr << "MSb = " << MSb.transpose() << '\n';
   ostr << "MSt = " << MSt.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZM = " << ZM << '\n';
   ostr << "ZTau = " << ZTau << '\n';
   ostr << "ZS = " << ZS << '\n';
   ostr << "ZC = " << ZC << '\n';
   ostr << "ZB = " << ZB << '\n';
   ostr << "ZT = " << ZT << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
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
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSt();
   calculate_MSb();
   calculate_MSc();
   calculate_MSs();
   calculate_MStau();
   calculate_MSm();
   calculate_MSe();
   calculate_MSu();
   calculate_MSd();
   calculate_MSvtL();
   calculate_MSvmL();
   calculate_MSveL();
   calculate_MFtau();
   calculate_MFm();
   calculate_MFe();
   calculate_MFvt();
   calculate_MFvm();
   calculate_MFve();
   calculate_MFt();
   calculate_MFc();
   calculate_MFu();
   calculate_MFb();
   calculate_MFs();
   calculate_MFd();
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
   PHYSICAL(MFd) = MFd;
   PHYSICAL(MFs) = MFs;
   PHYSICAL(MFb) = MFb;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(MFc) = MFc;
   PHYSICAL(MFt) = MFt;
   PHYSICAL(MFve) = MFve;
   PHYSICAL(MFvm) = MFvm;
   PHYSICAL(MFvt) = MFvt;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(MFm) = MFm;
   PHYSICAL(MFtau) = MFtau;
   PHYSICAL(MSveL) = MSveL;
   PHYSICAL(MSvmL) = MSvmL;
   PHYSICAL(MSvtL) = MSvtL;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSm) = MSm;
   PHYSICAL(ZM) = ZM;
   PHYSICAL(MStau) = MStau;
   PHYSICAL(ZTau) = ZTau;
   PHYSICAL(MSs) = MSs;
   PHYSICAL(ZS) = ZS;
   PHYSICAL(MSc) = MSc;
   PHYSICAL(ZC) = ZC;
   PHYSICAL(MSb) = MSb;
   PHYSICAL(ZB) = ZB;
   PHYSICAL(MSt) = MSt;
   PHYSICAL(ZT) = ZT;
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
   if (PHYSICAL(MSveL) < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::SveL); }
   if (PHYSICAL(MSvmL) < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::SvmL); }
   if (PHYSICAL(MSvtL) < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::SvtL); }
   if (PHYSICAL(MSd).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Sd); }
   if (PHYSICAL(MSu).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Su); }
   if (PHYSICAL(MSe).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Se); }
   if (PHYSICAL(MSm).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Sm); }
   if (PHYSICAL(MStau).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Stau); }
   if (PHYSICAL(MSs).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Ss); }
   if (PHYSICAL(MSc).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Sc); }
   if (PHYSICAL(MSb).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Sb); }
   if (PHYSICAL(MSt).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::St); }
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(CMSSMNoFV_info::Hpm); }
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
   MFd = 0.;
   MFs = 0.;
   MFb = 0.;
   MFu = 0.;
   MFc = 0.;
   MFt = 0.;
   MFve = 0.;
   MFvm = 0.;
   MFvt = 0.;
   MFe = 0.;
   MFm = 0.;
   MFtau = 0.;
   MSveL = 0.;
   MSvmL = 0.;
   MSvtL = 0.;
   MSd = Eigen::Matrix<double,2,1>::Zero();
   ZD = Eigen::Matrix<double,2,2>::Zero();
   MSu = Eigen::Matrix<double,2,1>::Zero();
   ZU = Eigen::Matrix<double,2,2>::Zero();
   MSe = Eigen::Matrix<double,2,1>::Zero();
   ZE = Eigen::Matrix<double,2,2>::Zero();
   MSm = Eigen::Matrix<double,2,1>::Zero();
   ZM = Eigen::Matrix<double,2,2>::Zero();
   MStau = Eigen::Matrix<double,2,1>::Zero();
   ZTau = Eigen::Matrix<double,2,2>::Zero();
   MSs = Eigen::Matrix<double,2,1>::Zero();
   ZS = Eigen::Matrix<double,2,2>::Zero();
   MSc = Eigen::Matrix<double,2,1>::Zero();
   ZC = Eigen::Matrix<double,2,2>::Zero();
   MSb = Eigen::Matrix<double,2,1>::Zero();
   ZB = Eigen::Matrix<double,2,2>::Zero();
   MSt = Eigen::Matrix<double,2,1>::Zero();
   ZT = Eigen::Matrix<double,2,2>::Zero();
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
   CMSSMNoFV_soft_parameters::clear();
   clear_tree_level_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFd = pars(2);
   MFs = pars(3);
   MFb = pars(4);
   MFu = pars(5);
   MFc = pars(6);
   MFt = pars(7);
   MFve = pars(8);
   MFvm = pars(9);
   MFvt = pars(10);
   MFe = pars(11);
   MFm = pars(12);
   MFtau = pars(13);
   MSveL = pars(14);
   MSvmL = pars(15);
   MSvtL = pars(16);
   MSd(0) = pars(17);
   MSd(1) = pars(18);
   MSu(0) = pars(19);
   MSu(1) = pars(20);
   MSe(0) = pars(21);
   MSe(1) = pars(22);
   MSm(0) = pars(23);
   MSm(1) = pars(24);
   MStau(0) = pars(25);
   MStau(1) = pars(26);
   MSs(0) = pars(27);
   MSs(1) = pars(28);
   MSc(0) = pars(29);
   MSc(1) = pars(30);
   MSb(0) = pars(31);
   MSb(1) = pars(32);
   MSt(0) = pars(33);
   MSt(1) = pars(34);
   Mhh(0) = pars(35);
   Mhh(1) = pars(36);
   MAh(0) = pars(37);
   MAh(1) = pars(38);
   MHpm(0) = pars(39);
   MHpm(1) = pars(40);
   MChi(0) = pars(41);
   MChi(1) = pars(42);
   MChi(2) = pars(43);
   MChi(3) = pars(44);
   MCha(0) = pars(45);
   MCha(1) = pars(46);
   MVWm = pars(47);
   MVP = pars(48);
   MVZ = pars(49);

}

const CMSSMNoFV_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

CMSSMNoFV_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   Eigen::ArrayXd pars(50);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFd;
   pars(3) = MFs;
   pars(4) = MFb;
   pars(5) = MFu;
   pars(6) = MFc;
   pars(7) = MFt;
   pars(8) = MFve;
   pars(9) = MFvm;
   pars(10) = MFvt;
   pars(11) = MFe;
   pars(12) = MFm;
   pars(13) = MFtau;
   pars(14) = MSveL;
   pars(15) = MSvmL;
   pars(16) = MSvtL;
   pars(17) = MSd(0);
   pars(18) = MSd(1);
   pars(19) = MSu(0);
   pars(20) = MSu(1);
   pars(21) = MSe(0);
   pars(22) = MSe(1);
   pars(23) = MSm(0);
   pars(24) = MSm(1);
   pars(25) = MStau(0);
   pars(26) = MStau(1);
   pars(27) = MSs(0);
   pars(28) = MSs(1);
   pars(29) = MSc(0);
   pars(30) = MSc(1);
   pars(31) = MSb(0);
   pars(32) = MSb(1);
   pars(33) = MSt(0);
   pars(34) = MSt(1);
   pars(35) = Mhh(0);
   pars(36) = Mhh(1);
   pars(37) = MAh(0);
   pars(38) = MAh(1);
   pars(39) = MHpm(0);
   pars(40) = MHpm(1);
   pars(41) = MChi(0);
   pars(42) = MChi(1);
   pars(43) = MChi(2);
   pars(44) = MChi(3);
   pars(45) = MCha(0);
   pars(46) = MCha(1);
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
   ZD(1,0) = pars(52);
   ZD(1,1) = pars(53);
   ZU(0,0) = pars(54);
   ZU(0,1) = pars(55);
   ZU(1,0) = pars(56);
   ZU(1,1) = pars(57);
   ZE(0,0) = pars(58);
   ZE(0,1) = pars(59);
   ZE(1,0) = pars(60);
   ZE(1,1) = pars(61);
   ZM(0,0) = pars(62);
   ZM(0,1) = pars(63);
   ZM(1,0) = pars(64);
   ZM(1,1) = pars(65);
   ZTau(0,0) = pars(66);
   ZTau(0,1) = pars(67);
   ZTau(1,0) = pars(68);
   ZTau(1,1) = pars(69);
   ZS(0,0) = pars(70);
   ZS(0,1) = pars(71);
   ZS(1,0) = pars(72);
   ZS(1,1) = pars(73);
   ZC(0,0) = pars(74);
   ZC(0,1) = pars(75);
   ZC(1,0) = pars(76);
   ZC(1,1) = pars(77);
   ZB(0,0) = pars(78);
   ZB(0,1) = pars(79);
   ZB(1,0) = pars(80);
   ZB(1,1) = pars(81);
   ZT(0,0) = pars(82);
   ZT(0,1) = pars(83);
   ZT(1,0) = pars(84);
   ZT(1,1) = pars(85);
   ZH(0,0) = pars(86);
   ZH(0,1) = pars(87);
   ZH(1,0) = pars(88);
   ZH(1,1) = pars(89);
   ZA(0,0) = pars(90);
   ZA(0,1) = pars(91);
   ZA(1,0) = pars(92);
   ZA(1,1) = pars(93);
   ZP(0,0) = pars(94);
   ZP(0,1) = pars(95);
   ZP(1,0) = pars(96);
   ZP(1,1) = pars(97);
   ZN(0,0) = std::complex<double>(pars(98), pars(99));
   ZN(0,1) = std::complex<double>(pars(100), pars(101));
   ZN(0,2) = std::complex<double>(pars(102), pars(103));
   ZN(0,3) = std::complex<double>(pars(104), pars(105));
   ZN(1,0) = std::complex<double>(pars(106), pars(107));
   ZN(1,1) = std::complex<double>(pars(108), pars(109));
   ZN(1,2) = std::complex<double>(pars(110), pars(111));
   ZN(1,3) = std::complex<double>(pars(112), pars(113));
   ZN(2,0) = std::complex<double>(pars(114), pars(115));
   ZN(2,1) = std::complex<double>(pars(116), pars(117));
   ZN(2,2) = std::complex<double>(pars(118), pars(119));
   ZN(2,3) = std::complex<double>(pars(120), pars(121));
   ZN(3,0) = std::complex<double>(pars(122), pars(123));
   ZN(3,1) = std::complex<double>(pars(124), pars(125));
   ZN(3,2) = std::complex<double>(pars(126), pars(127));
   ZN(3,3) = std::complex<double>(pars(128), pars(129));
   UM(0,0) = std::complex<double>(pars(130), pars(131));
   UM(0,1) = std::complex<double>(pars(132), pars(133));
   UM(1,0) = std::complex<double>(pars(134), pars(135));
   UM(1,1) = std::complex<double>(pars(136), pars(137));
   UP(0,0) = std::complex<double>(pars(138), pars(139));
   UP(0,1) = std::complex<double>(pars(140), pars(141));
   UP(1,0) = std::complex<double>(pars(142), pars(143));
   UP(1,1) = std::complex<double>(pars(144), pars(145));
   ZZ(0,0) = pars(146);
   ZZ(0,1) = pars(147);
   ZZ(1,0) = pars(148);
   ZZ(1,1) = pars(149);

}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_tree_level_masses());

   pars.conservativeResize(150);

   pars(50) = ZD(0,0);
   pars(51) = ZD(0,1);
   pars(52) = ZD(1,0);
   pars(53) = ZD(1,1);
   pars(54) = ZU(0,0);
   pars(55) = ZU(0,1);
   pars(56) = ZU(1,0);
   pars(57) = ZU(1,1);
   pars(58) = ZE(0,0);
   pars(59) = ZE(0,1);
   pars(60) = ZE(1,0);
   pars(61) = ZE(1,1);
   pars(62) = ZM(0,0);
   pars(63) = ZM(0,1);
   pars(64) = ZM(1,0);
   pars(65) = ZM(1,1);
   pars(66) = ZTau(0,0);
   pars(67) = ZTau(0,1);
   pars(68) = ZTau(1,0);
   pars(69) = ZTau(1,1);
   pars(70) = ZS(0,0);
   pars(71) = ZS(0,1);
   pars(72) = ZS(1,0);
   pars(73) = ZS(1,1);
   pars(74) = ZC(0,0);
   pars(75) = ZC(0,1);
   pars(76) = ZC(1,0);
   pars(77) = ZC(1,1);
   pars(78) = ZB(0,0);
   pars(79) = ZB(0,1);
   pars(80) = ZB(1,0);
   pars(81) = ZB(1,1);
   pars(82) = ZT(0,0);
   pars(83) = ZT(0,1);
   pars(84) = ZT(1,0);
   pars(85) = ZT(1,1);
   pars(86) = ZH(0,0);
   pars(87) = ZH(0,1);
   pars(88) = ZH(1,0);
   pars(89) = ZH(1,1);
   pars(90) = ZA(0,0);
   pars(91) = ZA(0,1);
   pars(92) = ZA(1,0);
   pars(93) = ZA(1,1);
   pars(94) = ZP(0,0);
   pars(95) = ZP(0,1);
   pars(96) = ZP(1,0);
   pars(97) = ZP(1,1);
   pars(98) = Re(ZN(0,0));
   pars(99) = Im(ZN(0,0));
   pars(100) = Re(ZN(0,1));
   pars(101) = Im(ZN(0,1));
   pars(102) = Re(ZN(0,2));
   pars(103) = Im(ZN(0,2));
   pars(104) = Re(ZN(0,3));
   pars(105) = Im(ZN(0,3));
   pars(106) = Re(ZN(1,0));
   pars(107) = Im(ZN(1,0));
   pars(108) = Re(ZN(1,1));
   pars(109) = Im(ZN(1,1));
   pars(110) = Re(ZN(1,2));
   pars(111) = Im(ZN(1,2));
   pars(112) = Re(ZN(1,3));
   pars(113) = Im(ZN(1,3));
   pars(114) = Re(ZN(2,0));
   pars(115) = Im(ZN(2,0));
   pars(116) = Re(ZN(2,1));
   pars(117) = Im(ZN(2,1));
   pars(118) = Re(ZN(2,2));
   pars(119) = Im(ZN(2,2));
   pars(120) = Re(ZN(2,3));
   pars(121) = Im(ZN(2,3));
   pars(122) = Re(ZN(3,0));
   pars(123) = Im(ZN(3,0));
   pars(124) = Re(ZN(3,1));
   pars(125) = Im(ZN(3,1));
   pars(126) = Re(ZN(3,2));
   pars(127) = Im(ZN(3,2));
   pars(128) = Re(ZN(3,3));
   pars(129) = Im(ZN(3,3));
   pars(130) = Re(UM(0,0));
   pars(131) = Im(UM(0,0));
   pars(132) = Re(UM(0,1));
   pars(133) = Im(UM(0,1));
   pars(134) = Re(UM(1,0));
   pars(135) = Im(UM(1,0));
   pars(136) = Re(UM(1,1));
   pars(137) = Im(UM(1,1));
   pars(138) = Re(UP(0,0));
   pars(139) = Im(UP(0,0));
   pars(140) = Re(UP(0,1));
   pars(141) = Im(UP(0,1));
   pars(142) = Re(UP(1,0));
   pars(143) = Im(UP(1,0));
   pars(144) = Re(UP(1,1));
   pars(145) = Im(UP(1,1));
   pars(146) = ZZ(0,0);
   pars(147) = ZZ(0,1);
   pars(148) = ZZ(1,0);
   pars(149) = ZZ(1,1);


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

double CLASSNAME::get_mass_matrix_Fd() const
{

   const double mass_matrix_Fd = Re(0.7071067811865475*vd*Yd(0,0));

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{

   const auto mass_matrix_Fd = get_mass_matrix_Fd();
   MFd = calculate_singlet_mass(mass_matrix_Fd);
}

double CLASSNAME::get_mass_matrix_Fs() const
{

   const double mass_matrix_Fs = Re(0.7071067811865475*vd*Yd(1,1));

   return mass_matrix_Fs;
}

void CLASSNAME::calculate_MFs()
{

   const auto mass_matrix_Fs = get_mass_matrix_Fs();
   MFs = calculate_singlet_mass(mass_matrix_Fs);
}

double CLASSNAME::get_mass_matrix_Fb() const
{

   const double mass_matrix_Fb = Re(0.7071067811865475*vd*Yd(2,2));

   return mass_matrix_Fb;
}

void CLASSNAME::calculate_MFb()
{

   const auto mass_matrix_Fb = get_mass_matrix_Fb();
   MFb = calculate_singlet_mass(mass_matrix_Fb);
}

double CLASSNAME::get_mass_matrix_Fu() const
{

   const double mass_matrix_Fu = Re(0.7071067811865475*vu*Yu(0,0));

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{

   const auto mass_matrix_Fu = get_mass_matrix_Fu();
   MFu = calculate_singlet_mass(mass_matrix_Fu);
}

double CLASSNAME::get_mass_matrix_Fc() const
{

   const double mass_matrix_Fc = Re(0.7071067811865475*vu*Yu(1,1));

   return mass_matrix_Fc;
}

void CLASSNAME::calculate_MFc()
{

   const auto mass_matrix_Fc = get_mass_matrix_Fc();
   MFc = calculate_singlet_mass(mass_matrix_Fc);
}

double CLASSNAME::get_mass_matrix_Ft() const
{

   const double mass_matrix_Ft = Re(0.7071067811865475*vu*Yu(2,2));

   return mass_matrix_Ft;
}

void CLASSNAME::calculate_MFt()
{

   const auto mass_matrix_Ft = get_mass_matrix_Ft();
   MFt = calculate_singlet_mass(mass_matrix_Ft);
}

double CLASSNAME::get_mass_matrix_Fve() const
{

   const double mass_matrix_Fve = Re(0);

   return mass_matrix_Fve;
}

void CLASSNAME::calculate_MFve()
{

   const auto mass_matrix_Fve = get_mass_matrix_Fve();
   MFve = calculate_singlet_mass(mass_matrix_Fve);
}

double CLASSNAME::get_mass_matrix_Fvm() const
{

   const double mass_matrix_Fvm = Re(0);

   return mass_matrix_Fvm;
}

void CLASSNAME::calculate_MFvm()
{

   const auto mass_matrix_Fvm = get_mass_matrix_Fvm();
   MFvm = calculate_singlet_mass(mass_matrix_Fvm);
}

double CLASSNAME::get_mass_matrix_Fvt() const
{

   const double mass_matrix_Fvt = Re(0);

   return mass_matrix_Fvt;
}

void CLASSNAME::calculate_MFvt()
{

   const auto mass_matrix_Fvt = get_mass_matrix_Fvt();
   MFvt = calculate_singlet_mass(mass_matrix_Fvt);
}

double CLASSNAME::get_mass_matrix_Fe() const
{

   const double mass_matrix_Fe = Re(0.7071067811865475*vd*Ye(0,0));

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{

   const auto mass_matrix_Fe = get_mass_matrix_Fe();
   MFe = calculate_singlet_mass(mass_matrix_Fe);
}

double CLASSNAME::get_mass_matrix_Fm() const
{

   const double mass_matrix_Fm = Re(0.7071067811865475*vd*Ye(1,1));

   return mass_matrix_Fm;
}

void CLASSNAME::calculate_MFm()
{

   const auto mass_matrix_Fm = get_mass_matrix_Fm();
   MFm = calculate_singlet_mass(mass_matrix_Fm);
}

double CLASSNAME::get_mass_matrix_Ftau() const
{

   const double mass_matrix_Ftau = Re(0.7071067811865475*vd*Ye(2,2));

   return mass_matrix_Ftau;
}

void CLASSNAME::calculate_MFtau()
{

   const auto mass_matrix_Ftau = get_mass_matrix_Ftau();
   MFtau = calculate_singlet_mass(mass_matrix_Ftau);
}

double CLASSNAME::get_mass_matrix_SveL() const
{

   const double mass_matrix_SveL = Re(0.125*(8*ml2(0,0) - 0.6*Sqr(g1)*(-Sqr(vd)
      + Sqr(vu)) - Sqr(g2)*(-Sqr(vd) + Sqr(vu))));

   return mass_matrix_SveL;
}

void CLASSNAME::calculate_MSveL()
{

   const auto mass_matrix_SveL = get_mass_matrix_SveL();
   MSveL = mass_matrix_SveL;

   if (MSveL < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::SveL);
   }

   MSveL = AbsSqrt(MSveL);
}

double CLASSNAME::get_mass_matrix_SvmL() const
{

   const double mass_matrix_SvmL = Re(0.125*(8*ml2(1,1) - 0.6*Sqr(g1)*(-Sqr(vd)
      + Sqr(vu)) - Sqr(g2)*(-Sqr(vd) + Sqr(vu))));

   return mass_matrix_SvmL;
}

void CLASSNAME::calculate_MSvmL()
{

   const auto mass_matrix_SvmL = get_mass_matrix_SvmL();
   MSvmL = mass_matrix_SvmL;

   if (MSvmL < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::SvmL);
   }

   MSvmL = AbsSqrt(MSvmL);
}

double CLASSNAME::get_mass_matrix_SvtL() const
{

   const double mass_matrix_SvtL = Re(0.125*(8*ml2(2,2) - 0.6*Sqr(g1)*(-Sqr(vd)
      + Sqr(vu)) - Sqr(g2)*(-Sqr(vd) + Sqr(vu))));

   return mass_matrix_SvtL;
}

void CLASSNAME::calculate_MSvtL()
{

   const auto mass_matrix_SvtL = get_mass_matrix_SvtL();
   MSvtL = mass_matrix_SvtL;

   if (MSvtL < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::SvtL);
   }

   MSvtL = AbsSqrt(MSvtL);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*
      Sqr(vu);
   mass_matrix_Sd(0,1) = 0.7071067811865475*vd*Conj(TYd(0,0)) -
      0.7071067811865475*vu*Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(1,1) = md2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Sd, eigenvalue_error > precision *
      Abs(MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*AbsSqr(Yu(0,0))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)
      *Sqr(vu);
   mass_matrix_Su(0,1) = 0.7071067811865475*vu*Conj(TYu(0,0)) -
      0.7071067811865475*vd*Conj(Yu(0,0))*Mu;
   mass_matrix_Su(1,1) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*AbsSqr(Yu(0,0))*
      Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Su, eigenvalue_error > precision *
      Abs(MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*
      Sqr(vu);
   mass_matrix_Se(0,1) = 0.7071067811865475*vd*Conj(TYe(0,0)) -
      0.7071067811865475*vu*Conj(Ye(0,0))*Mu;
   mass_matrix_Se(1,1) = me2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Se, eigenvalue_error > precision *
      Abs(MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Sm;

   mass_matrix_Sm(0,0) = ml2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*
      Sqr(vu);
   mass_matrix_Sm(0,1) = 0.7071067811865475*vd*Conj(TYe(1,1)) -
      0.7071067811865475*vu*Conj(Ye(1,1))*Mu;
   mass_matrix_Sm(1,1) = me2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sm);

   return mass_matrix_Sm;
}

void CLASSNAME::calculate_MSm()
{
   const auto mass_matrix_Sm(get_mass_matrix_Sm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sm, MSm, ZM, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Sm, eigenvalue_error > precision *
      Abs(MSm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sm, MSm, ZM);
#endif
   normalize_to_interval(ZM);


   if (MSm.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Sm);
   }

   MSm = AbsSqrt(MSm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Stau() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Stau;

   mass_matrix_Stau(0,0) = ml2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) + 0.075*Sqr(
      g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(
      g2)*Sqr(vu);
   mass_matrix_Stau(0,1) = 0.7071067811865475*vd*Conj(TYe(2,2)) -
      0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Stau(1,1) = me2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) - 0.15*Sqr(g1
      )*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Stau);

   return mass_matrix_Stau;
}

void CLASSNAME::calculate_MStau()
{
   const auto mass_matrix_Stau(get_mass_matrix_Stau());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Stau, MStau, ZTau, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Stau, eigenvalue_error > precision *
      Abs(MStau(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Stau, MStau, ZTau);
#endif
   normalize_to_interval(ZTau);


   if (MStau.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Stau);
   }

   MStau = AbsSqrt(MStau);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ss() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Ss;

   mass_matrix_Ss(0,0) = mq2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*
      Sqr(vu);
   mass_matrix_Ss(0,1) = 0.7071067811865475*vd*Conj(TYd(1,1)) -
      0.7071067811865475*vu*Conj(Yd(1,1))*Mu;
   mass_matrix_Ss(1,1) = md2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Ss);

   return mass_matrix_Ss;
}

void CLASSNAME::calculate_MSs()
{
   const auto mass_matrix_Ss(get_mass_matrix_Ss());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ss, MSs, ZS, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Ss, eigenvalue_error > precision *
      Abs(MSs(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ss, MSs, ZS);
#endif
   normalize_to_interval(ZS);


   if (MSs.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Ss);
   }

   MSs = AbsSqrt(MSs);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sc() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Sc;

   mass_matrix_Sc(0,0) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*AbsSqr(Yu(1,1))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)
      *Sqr(vu);
   mass_matrix_Sc(0,1) = 0.7071067811865475*vu*Conj(TYu(1,1)) -
      0.7071067811865475*vd*Conj(Yu(1,1))*Mu;
   mass_matrix_Sc(1,1) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*AbsSqr(Yu(1,1))*
      Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sc);

   return mass_matrix_Sc;
}

void CLASSNAME::calculate_MSc()
{
   const auto mass_matrix_Sc(get_mass_matrix_Sc());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sc, MSc, ZC, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Sc, eigenvalue_error > precision *
      Abs(MSc(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sc, MSc, ZC);
#endif
   normalize_to_interval(ZC);


   if (MSc.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Sc);
   }

   MSc = AbsSqrt(MSc);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sb() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Sb;

   mass_matrix_Sb(0,0) = mq2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*
      Sqr(vu);
   mass_matrix_Sb(0,1) = 0.7071067811865475*vd*Conj(TYd(2,2)) -
      0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sb(1,1) = md2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sb);

   return mass_matrix_Sb;
}

void CLASSNAME::calculate_MSb()
{
   const auto mass_matrix_Sb(get_mass_matrix_Sb());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sb, MSb, ZB, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::Sb, eigenvalue_error > precision *
      Abs(MSb(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sb, MSb, ZB);
#endif
   normalize_to_interval(ZB);


   if (MSb.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Sb);
   }

   MSb = AbsSqrt(MSb);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_St() const
{

   Eigen::Matrix<double,2,2> mass_matrix_St;

   mass_matrix_St(0,0) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*AbsSqr(Yu(2,2))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)
      *Sqr(vu);
   mass_matrix_St(0,1) = 0.7071067811865475*vu*Conj(TYu(2,2)) -
      0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_St(1,1) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*AbsSqr(Yu(2,2))*
      Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_St);

   return mass_matrix_St;
}

void CLASSNAME::calculate_MSt()
{
   const auto mass_matrix_St(get_mass_matrix_St());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_St, MSt, ZT, eigenvalue_error);
   problems.flag_bad_mass(CMSSMNoFV_info::St, eigenvalue_error > precision *
      Abs(MSt(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_St, MSt, ZT);
#endif
   normalize_to_interval(ZT);


   if (MSt.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::St);
   }

   MSt = AbsSqrt(MSt);
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
   problems.flag_bad_mass(CMSSMNoFV_info::hh, eigenvalue_error > precision *
      Abs(Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::hh);
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
   problems.flag_bad_mass(CMSSMNoFV_info::Ah, eigenvalue_error > precision *
      Abs(MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Ah);
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
   problems.flag_bad_mass(CMSSMNoFV_info::Hpm, eigenvalue_error > precision *
      Abs(MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(CMSSMNoFV_info::Hpm);
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
   problems.flag_bad_mass(CMSSMNoFV_info::Chi, eigenvalue_error > precision *
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
   problems.flag_bad_mass(CMSSMNoFV_info::Cha, eigenvalue_error > precision *
      Abs(MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
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
      problems.flag_running_tachyon(CMSSMNoFV_info::VWm);
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



std::ostream& operator<<(std::ostream& ostr, const CMSSMNoFV_mass_eigenstates_decoupling_scheme& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
