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

// File generated at Wed 16 Oct 2019 19:02:45

#include "NUHMSSMNoFVHimalaya_two_scale_high_scale_constraint.hpp"
#include "NUHMSSMNoFVHimalaya_two_scale_model.hpp"
#include "NUHMSSMNoFVHimalaya_info.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#include "version.hpp"
#endif

#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME NUHMSSMNoFVHimalaya<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto TYu = MODELPARAMETER(TYu); \
      const auto TYd = MODELPARAMETER(TYd); \
      const auto TYe = MODELPARAMETER(TYe); \
      const auto Mu = MODELPARAMETER(Mu); \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto vu = MODELPARAMETER(vu); \
      const auto vd = MODELPARAMETER(vd); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
      const auto MGlu = MODELPARAMETER(MGlu); \
      const auto MAh = MODELPARAMETER(MAh); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(Mu); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(vd); \
      pars.vu = Re(vu); \
      pars.mq2 = Re(UpMatrixL); \
      pars.md2 = Re(DownMatrixR); \
      pars.mu2 = Re(UpMatrixR); \
      pars.ml2 = Re(ElectronMatrixL); \
      pars.me2 = Re(ElectronMatrixR); \
      pars.Au = Re(TYu); \
      pars.Ad = Re(TYd); \
      pars.Ae = Re(TYe); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = 0; \
      pars.M2 = 0; \
      pars.MG = MGlu; \
      pars.MA = MAh; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = 1; \
      const double lambda_3L_uncertainty = 0; \
       \
                                                                        \
      double lambda_3L = 0.;                                            \
                                                                        \
      try {                                                             \
         const bool verbose = false;                                    \
         himalaya::HierarchyCalculator hc(pars, verbose);               \
                                                                        \
         const auto ho = hc.calculateDMh3L(false);                      \
                                                                        \
         lambda_3L =                                                    \
            lambda_3L_eft * (                                           \
               ho.getDLambda(3)                                         \
               + msbar_scheme*ho.getDLambdaDRbarPrimeToMSbarShift(3)    \
               + lambda_3L_uncertainty*ho.getDLambdaUncertainty(3)      \
            );                                                          \
                                                                        \
         VERBOSE_MSG("Himalaya top (hierarchy, Dlambda_3L) = ("         \
                     << ho.getSuitableHierarchy() << ", "               \
                     << lambda_3L <<")");                               \
      } catch (const std::exception& e) {                               \
         model->get_problems().flag_bad_mass(NUHMSSMNoFVHimalaya_info::hh); \
         WARNING(e.what());                                             \
         VERBOSE_MSG(pars);                                             \
      }                                                                 \
                                                                        \
      return lambda_3L;                                                 \
   }()
#else
#define FSHimalayaMh23L [] () {                                         \
      throw HimalayaError("The 3-loop corrections to lambda "           \
                          "require Himalaya 2.0.0 (or higher), but "    \
                          "FlexibleSUSY has not been configured with "  \
                          "this Himalaya version!");                    \
      return 0.;                                                        \
   }()
#endif

NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::NUHMSSMNoFVHimalaya_high_scale_constraint(
   NUHMSSMNoFVHimalaya<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   update_scale();

   const auto AeIN = INPUTPARAMETER(AeIN);
   const auto AmuonIN = INPUTPARAMETER(AmuonIN);
   const auto AtauIN = INPUTPARAMETER(AtauIN);
   const auto AdIN = INPUTPARAMETER(AdIN);
   const auto AsIN = INPUTPARAMETER(AsIN);
   const auto AbIN = INPUTPARAMETER(AbIN);
   const auto AuIN = INPUTPARAMETER(AuIN);
   const auto AcIN = INPUTPARAMETER(AcIN);
   const auto AtIN = INPUTPARAMETER(AtIN);
   const auto MuIN = INPUTPARAMETER(MuIN);
   const auto mA2IN = INPUTPARAMETER(mA2IN);
   const auto mq11IN = INPUTPARAMETER(mq11IN);
   const auto mq22IN = INPUTPARAMETER(mq22IN);
   const auto mq33IN = INPUTPARAMETER(mq33IN);
   const auto ml11IN = INPUTPARAMETER(ml11IN);
   const auto ml22IN = INPUTPARAMETER(ml22IN);
   const auto ml33IN = INPUTPARAMETER(ml33IN);
   const auto md11IN = INPUTPARAMETER(md11IN);
   const auto md22IN = INPUTPARAMETER(md22IN);
   const auto md33IN = INPUTPARAMETER(md33IN);
   const auto mu11IN = INPUTPARAMETER(mu11IN);
   const auto mu22IN = INPUTPARAMETER(mu22IN);
   const auto mu33IN = INPUTPARAMETER(mu33IN);
   const auto me11IN = INPUTPARAMETER(me11IN);
   const auto me22IN = INPUTPARAMETER(me22IN);
   const auto me33IN = INPUTPARAMETER(me33IN);
   const auto M1 = INPUTPARAMETER(M1);
   const auto M2 = INPUTPARAMETER(M2);
   const auto M3 = INPUTPARAMETER(M3);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_TYe(0,0,Re(AeIN*Ye(0,0)));
   MODEL->set_TYe(1,1,Re(AmuonIN*Ye(1,1)));
   MODEL->set_TYe(2,2,Re(AtauIN*Ye(2,2)));
   MODEL->set_TYd(0,0,Re(AdIN*Yd(0,0)));
   MODEL->set_TYd(1,1,Re(AsIN*Yd(1,1)));
   MODEL->set_TYd(2,2,Re(AbIN*Yd(2,2)));
   MODEL->set_TYu(0,0,Re(AuIN*Yu(0,0)));
   MODEL->set_TYu(1,1,Re(AcIN*Yu(1,1)));
   MODEL->set_TYu(2,2,Re(AtIN*Yu(2,2)));
   MODEL->set_Mu(Re(MuIN));
   MODEL->set_BMu(Re(mA2IN/(vd/vu + vu/vd)));
   MODEL->set_mq2(0,0,Re(Sqr(mq11IN)));
   MODEL->set_mq2(1,1,Re(Sqr(mq22IN)));
   MODEL->set_mq2(2,2,Re(Sqr(mq33IN)));
   MODEL->set_ml2(0,0,Re(Sqr(ml11IN)));
   MODEL->set_ml2(1,1,Re(Sqr(ml22IN)));
   MODEL->set_ml2(2,2,Re(Sqr(ml33IN)));
   MODEL->set_md2(0,0,Re(Sqr(md11IN)));
   MODEL->set_md2(1,1,Re(Sqr(md22IN)));
   MODEL->set_md2(2,2,Re(Sqr(md33IN)));
   MODEL->set_mu2(0,0,Re(Sqr(mu11IN)));
   MODEL->set_mu2(1,1,Re(Sqr(mu22IN)));
   MODEL->set_mu2(2,2,Re(Sqr(mu33IN)));
   MODEL->set_me2(0,0,Re(Sqr(me11IN)));
   MODEL->set_me2(1,1,Re(Sqr(me22IN)));
   MODEL->set_me2(2,2,Re(Sqr(me33IN)));
   MODEL->set_MassB(Re(M1));
   MODEL->set_MassWB(Re(M2));
   MODEL->set_MassG(Re(M3));


   check_non_perturbative();
}

bool NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yu = MODELPARAMETER(Yu);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::g3);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Ye2_2);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(NUHMSSMNoFVHimalaya_info::Yu2_2);
   }


   return problem;
}

double NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const NUHMSSMNoFVHimalaya_input_parameters& NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

NUHMSSMNoFVHimalaya<Two_scale>* NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<NUHMSSMNoFVHimalaya<Two_scale>*>(model_);
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto MuIN = INPUTPARAMETER(MuIN);
   const auto AtIN = INPUTPARAMETER(AtIN);
   const auto mq33IN = INPUTPARAMETER(mq33IN);
   const auto mu33IN = INPUTPARAMETER(mu33IN);

   initial_scale_guess = IF(Qin != 0, Qin, Sqrt(Sqrt((51200*AtIN*MuIN*TanBeta -
      25600*Sqr(MuIN) + (-25600*Sqr(AtIN) + (25600 + Sqr(mq33IN))*(25600 + Sqr(
      mu33IN)))*Sqr(TanBeta))/Sqr(TanBeta))));

   scale = initial_scale_guess;
}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto mq33IN = INPUTPARAMETER(mq33IN);
   const auto mu33IN = INPUTPARAMETER(mu33IN);
   const auto MuIN = INPUTPARAMETER(MuIN);
   const auto AtIN = INPUTPARAMETER(AtIN);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);

   scale = IF(Qin != 0, Qin, 0.7071067811865476*Sqrt(Sqrt((2*Sqr(mq33IN)*Sqr(
      TanBeta)*(2*Sqr(mu33IN) + Sqr(vu)*Sqr(Yu(2,2))) + Sqr(vu)*Sqr(Yu(2,2))*(4*
      AtIN*MuIN*TanBeta - 2*Sqr(MuIN) + Sqr(TanBeta)*(-2*Sqr(AtIN) + 2*Sqr(mu33IN)
      + Sqr(vu)*Sqr(Yu(2,2)))))/Sqr(TanBeta))));


}

void NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("NUHMSSMNoFVHimalaya_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
