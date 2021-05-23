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


#include "MRSSM2_slha_io.hpp"
#include "MRSSM2_input_parameters.hpp"
#include "MRSSM2_mass_eigenstates.hpp"
#include "MRSSM2_model_slha.hpp"
#include "MRSSM2_observables.hpp"
#include "MRSSM2_physical.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "observable_problems.hpp"
#include "observable_problems_format_slha.hpp"
#include "numerics2.hpp"
#include "spectrum_generator_problems.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"
#include "config.h"
#include "spectrum_generator_settings.hpp"

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <string>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) input.p
#define EXTRAPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

MRSSM2_slha_io::MRSSM2_slha_io()
   : slha_io()
   , print_imaginary_parts_of_majorana_mixings(false)
{
}

void MRSSM2_slha_io::clear()
{
   slha_io.clear();
}

void MRSSM2_slha_io::set_print_imaginary_parts_of_majorana_mixings(bool flag)
{
   print_imaginary_parts_of_majorana_mixings = flag;
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
void MRSSM2_slha_io::fill(MRSSM2_slha& model) const
{
   fill(static_cast<MRSSM2_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM2_slha_io::set_extpar(const MRSSM2_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(0, input.Ms, "Ms");
   slha_io.set_block(extpar);

}

/**
 * Stores the IMMINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM2_slha_io::set_imminpar(const MRSSM2_input_parameters& input)
{

}

/**
 * Stores the IMEXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM2_slha_io::set_imextpar(const MRSSM2_input_parameters& input)
{

}

/**
 * Stores the MODSEL input parameters in the SLHA object.
 *
 * @param modsel struct of MODSEL parameters
 */
void MRSSM2_slha_io::set_modsel(const SLHA_io::Modsel& modsel)
{
   slha_io.set_modsel(modsel);
}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM2_slha_io::set_minpar(const MRSSM2_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   slha_io.set_block(minpar);

}

/**
 * Stores all input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM2_slha_io::set_input(const MRSSM2_input_parameters& input)
{
   set_minpar(input);
   set_extpar(input);
   set_imminpar(input);
   set_imextpar(input);

   {
      std::ostringstream block;
      block << "Block HMIXIN" << '\n'
            << FORMAT_ELEMENT(101, (INPUTPARAMETER(BMuInput)), "BMuInput")
            << FORMAT_ELEMENT(301, (INPUTPARAMETER(LamSDInput)), "LamSDInput")
            << FORMAT_ELEMENT(302, (INPUTPARAMETER(LamSUInput)), "LamSUInput")
            << FORMAT_ELEMENT(303, (INPUTPARAMETER(LamTDInput)), "LamTDInput")
            << FORMAT_ELEMENT(304, (INPUTPARAMETER(LamTUInput)), "LamTUInput")
            << FORMAT_ELEMENT(201, (INPUTPARAMETER(MuDInput)), "MuDInput")
            << FORMAT_ELEMENT(202, (INPUTPARAMETER(MuUInput)), "MuUInput")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("MSD2IN", INPUTPARAMETER(md2Input), "md2Input");
   {
      std::ostringstream block;
      block << "Block MSOFTIN" << '\n'
            << FORMAT_ELEMENT(300, (INPUTPARAMETER(MDBSInput)), "MDBSInput")
            << FORMAT_ELEMENT(302, (INPUTPARAMETER(MDGocInput)), "MDGocInput")
            << FORMAT_ELEMENT(301, (INPUTPARAMETER(MDWBTInput)), "MDWBTInput")
            << FORMAT_ELEMENT(111, (INPUTPARAMETER(moc2Input)), "moc2Input")
            << FORMAT_ELEMENT(50, (INPUTPARAMETER(mRd2Input)), "mRd2Input")
            << FORMAT_ELEMENT(51, (INPUTPARAMETER(mRu2Input)), "mRu2Input")
            << FORMAT_ELEMENT(110, (INPUTPARAMETER(mT2Input)), "mT2Input")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("MSE2IN", INPUTPARAMETER(me2Input), "me2Input");
   slha_io.set_block("MSL2IN", INPUTPARAMETER(ml2Input), "ml2Input");
   slha_io.set_block("MSQ2IN", INPUTPARAMETER(mq2Input), "mq2Input");
   {
      std::ostringstream block;
      block << "Block NMSSMRUNIN" << '\n'
            << FORMAT_ELEMENT(10, (INPUTPARAMETER(mS2Input)), "mS2Input")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("MSU2IN", INPUTPARAMETER(mu2Input), "mu2Input");

}

/**
 * Stores the additional physical input (FlexibleSUSYInput block) in
 * the SLHA object.
 *
 * @param input class of input
 */
void MRSSM2_slha_io::set_physical_input(const Physical_input& input)
{
   slha_io.set_physical_input(input);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void MRSSM2_slha_io::set_settings(const Spectrum_generator_settings& settings)
{
   slha_io.set_settings(settings);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void MRSSM2_slha_io::set_FlexibleDecay_settings(const FlexibleDecay_settings& settings)
{
   slha_io.set_FlexibleDecay_settings(settings);
}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void MRSSM2_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void MRSSM2_slha_io::set_spinfo(const Spectrum_generator_problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void MRSSM2_slha_io::set_spinfo(const Problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the given problems and warnings in the SPINFO block in the
 * SLHA object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void MRSSM2_slha_io::set_spinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      spinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      spinfo << FORMAT_SPINFO(4, s);

   spinfo << FORMAT_SPINFO(5, MRSSM2_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void MRSSM2_slha_io::set_mass(const MRSSM2_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(404, LOCALPHYSICAL(MSRdp), "SRdp")
      << FORMAT_MASS(403, LOCALPHYSICAL(MSRum), "SRum")
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(3000022, LOCALPHYSICAL(MsigmaO), "sigmaO")
      << FORMAT_MASS(3000021, LOCALPHYSICAL(MphiO), "phiO")
      << FORMAT_MASS(2000024, LOCALPHYSICAL(MCha2(0)), "Cha2(1)")
      << FORMAT_MASS(2000037, LOCALPHYSICAL(MCha2(1)), "Cha2(2)")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha1(0)), "Cha1(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha1(1)), "Cha1(2)")
      << FORMAT_MASS(401, LOCALPHYSICAL(MRh(0)), "Rh(1)")
      << FORMAT_MASS(402, LOCALPHYSICAL(MRh(1)), "Rh(2)")
      << FORMAT_MASS(1000012, LOCALPHYSICAL(MSv(0)), "Sv(1)")
      << FORMAT_MASS(1000014, LOCALPHYSICAL(MSv(1)), "Sv(2)")
      << FORMAT_MASS(1000016, LOCALPHYSICAL(MSv(2)), "Sv(3)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHpm(1)), "Hpm(2)")
      << FORMAT_MASS(47, LOCALPHYSICAL(MHpm(2)), "Hpm(3)")
      << FORMAT_MASS(57, LOCALPHYSICAL(MHpm(3)), "Hpm(4)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(45, LOCALPHYSICAL(Mhh(2)), "hh(3)")
      << FORMAT_MASS(55, LOCALPHYSICAL(Mhh(3)), "hh(4)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(1)), "Ah(2)")
      << FORMAT_MASS(46, LOCALPHYSICAL(MAh(2)), "Ah(3)")
      << FORMAT_MASS(56, LOCALPHYSICAL(MAh(3)), "Ah(4)")
      << FORMAT_MASS(1000001, LOCALPHYSICAL(MSd(0)), "Sd(1)")
      << FORMAT_MASS(1000003, LOCALPHYSICAL(MSd(1)), "Sd(2)")
      << FORMAT_MASS(1000005, LOCALPHYSICAL(MSd(2)), "Sd(3)")
      << FORMAT_MASS(2000001, LOCALPHYSICAL(MSd(3)), "Sd(4)")
      << FORMAT_MASS(2000003, LOCALPHYSICAL(MSd(4)), "Sd(5)")
      << FORMAT_MASS(2000005, LOCALPHYSICAL(MSd(5)), "Sd(6)")
      << FORMAT_MASS(1000011, LOCALPHYSICAL(MSe(0)), "Se(1)")
      << FORMAT_MASS(1000013, LOCALPHYSICAL(MSe(1)), "Se(2)")
      << FORMAT_MASS(1000015, LOCALPHYSICAL(MSe(2)), "Se(3)")
      << FORMAT_MASS(2000011, LOCALPHYSICAL(MSe(3)), "Se(4)")
      << FORMAT_MASS(2000013, LOCALPHYSICAL(MSe(4)), "Se(5)")
      << FORMAT_MASS(2000015, LOCALPHYSICAL(MSe(5)), "Se(6)")
      << FORMAT_MASS(1000002, LOCALPHYSICAL(MSu(0)), "Su(1)")
      << FORMAT_MASS(1000004, LOCALPHYSICAL(MSu(1)), "Su(2)")
      << FORMAT_MASS(1000006, LOCALPHYSICAL(MSu(2)), "Su(3)")
      << FORMAT_MASS(2000002, LOCALPHYSICAL(MSu(3)), "Su(4)")
      << FORMAT_MASS(2000004, LOCALPHYSICAL(MSu(4)), "Su(5)")
      << FORMAT_MASS(2000006, LOCALPHYSICAL(MSu(5)), "Su(6)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void MRSSM2_slha_io::set_mixing_matrices(const MRSSM2_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("U1MIX", LOCALPHYSICAL(UM1), "UM1");
   slha_io.set_block("U2MIX", LOCALPHYSICAL(UM2), "UM2");
   slha_io.set_block("V1MIX", LOCALPHYSICAL(UP1), "UP1");
   slha_io.set_block("V2MIX", LOCALPHYSICAL(UP2), "UP2");
   slha_io.set_block("PSEUDOSCALARMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", LOCALPHYSICAL(ZD), "ZD");
   slha_io.set_block("SELMIX", LOCALPHYSICAL(ZE), "ZE");
   slha_io.set_block("SCALARMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("RHMIX", LOCALPHYSICAL(ZHR), "ZHR");
   slha_io.set_block("N1MIX", LOCALPHYSICAL(ZN1), "ZN1");
   slha_io.set_block("N2MIX", LOCALPHYSICAL(ZN2), "ZN2");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");
   slha_io.set_block("USQMIX", LOCALPHYSICAL(ZU), "ZU");
   slha_io.set_block("SNUMIX", LOCALPHYSICAL(ZV), "ZV");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UELMIX", LOCALPHYSICAL(ZEL), "ZEL");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(ZER), "ZER");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(ZDL), "ZDL");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(ZDR), "ZDR");
      slha_io.set_block("UULMIX", LOCALPHYSICAL(ZUL), "ZUL");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(ZUR), "ZUR");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
   }

}

void MRSSM2_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void MRSSM2_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
void MRSSM2_slha_io::set_model_parameters(const MRSSM2_slha& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
            << FORMAT_ELEMENT(310, (MODELPARAMETER(vT)), "vT")
            << FORMAT_ELEMENT(201, (MODELPARAMETER(MuD)), "MuD")
            << FORMAT_ELEMENT(202, (MODELPARAMETER(MuU)), "MuU")
            << FORMAT_ELEMENT(203, (MODELPARAMETER(BMuD)), "BMuD")
            << FORMAT_ELEMENT(204, (MODELPARAMETER(BMuU)), "BMuU")
            << FORMAT_ELEMENT(301, (MODELPARAMETER(LamSD)), "LamSD")
            << FORMAT_ELEMENT(302, (MODELPARAMETER(LamSU)), "LamSU")
            << FORMAT_ELEMENT(303, (MODELPARAMETER(LamTD)), "LamTD")
            << FORMAT_ELEMENT(304, (MODELPARAMETER(LamTU)), "LamTU")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2_slha), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2_slha), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2_slha), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2_slha), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2_slha), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(110, (MODELPARAMETER(mT2)), "mT2")
            << FORMAT_ELEMENT(111, (MODELPARAMETER(moc2)), "moc2")
            << FORMAT_ELEMENT(300, (MODELPARAMETER(MDBS)), "MDBS")
            << FORMAT_ELEMENT(301, (MODELPARAMETER(MDWBT)), "MDWBT")
            << FORMAT_ELEMENT(302, (MODELPARAMETER(MDGoc)), "MDGoc")
            << FORMAT_ELEMENT(50, (MODELPARAMETER(mRd2)), "mRd2")
            << FORMAT_ELEMENT(51, (MODELPARAMETER(mRu2)), "mRu2")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block NMSSMRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(10, (MODELPARAMETER(mS2)), "mS2")
            << FORMAT_ELEMENT(5, (MODELPARAMETER(vS)), "vS")
      ;
      slha_io.set_block(block);
   }


}

void MRSSM2_slha_io::set_model_parameters(const standard_model::Standard_model& model)
{
   {
      std::ostringstream block;
      block << "Block SMGAUGE Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_g1() * standard_model_info::normalization_g1), "gY")
            << FORMAT_ELEMENT(2, (model.get_g2()), "g2")
            << FORMAT_ELEMENT(3, (model.get_g3()), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("SMYu", ToMatrix(model.get_Yu()), "Yu", model.get_scale());
   slha_io.set_block("SMYd", ToMatrix(model.get_Yd()), "Yd", model.get_scale());
   slha_io.set_block("SMYe", ToMatrix(model.get_Ye()), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SMSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_mu2()), "mu2")
            << FORMAT_ELEMENT(2, (model.get_Lambdax()), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SMHMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (model.get_v()), "v")
      ;
      slha_io.set_block(block);
   }
}

void MRSSM2_slha_io::set_mass(const standard_model::Standard_model_physical& physical)
{
   std::ostringstream mass;

   mass << "Block SMMASS\n"
      << FORMAT_MASS(24, physical.MVWp, "VWp")
      << FORMAT_MASS(21, physical.MVG, "VG")
      << FORMAT_MASS(12, physical.MFv(0), "Fv(1)")
      << FORMAT_MASS(14, physical.MFv(1), "Fv(2)")
      << FORMAT_MASS(16, physical.MFv(2), "Fv(3)")
      << FORMAT_MASS(25, physical.Mhh, "hh")
      << FORMAT_MASS(1, physical.MFd(0), "Fd(1)")
      << FORMAT_MASS(3, physical.MFd(1), "Fd(2)")
      << FORMAT_MASS(5, physical.MFd(2), "Fd(3)")
      << FORMAT_MASS(2, physical.MFu(0), "Fu(1)")
      << FORMAT_MASS(4, physical.MFu(1), "Fu(2)")
      << FORMAT_MASS(6, physical.MFu(2), "Fu(3)")
      << FORMAT_MASS(11, physical.MFe(0), "Fe(1)")
      << FORMAT_MASS(13, physical.MFe(1), "Fe(2)")
      << FORMAT_MASS(15, physical.MFe(2), "Fe(3)")
      << FORMAT_MASS(22, physical.MVP, "VP")
      << FORMAT_MASS(23, physical.MVZ, "VZ")
      ;

   slha_io.set_block(mass);
}

void MRSSM2_slha_io::set_mixing_matrices(const standard_model::Standard_model_physical& physical)
{
   slha_io.set_block("SMUULMIX", physical.Vu, "Vu");
   slha_io.set_block("SMUDLMIX", physical.Vd, "Vd");
   slha_io.set_block("SMUURMIX", physical.Uu, "Uu");
   slha_io.set_block("SMUDRMIX", physical.Ud, "Ud");
   slha_io.set_block("SMUELMIX", physical.Ve, "Ve");
   slha_io.set_block("SMUERMIX", physical.Ue, "Ue");
}

void MRSSM2_slha_io::set_spectrum(const standard_model::Standard_model& model)
{
   const auto& physical = model.get_physical();

   set_model_parameters(model);
   set_mass(physical);
   set_mixing_matrices(physical);
}

/**
 * Stores the decays calculation information in the DCINFO block
 * in the SLHA object.
 *
 * @param problems struct with decays calculation problems
 */
void MRSSM2_slha_io::set_dcinfo(const FlexibleDecay_problems& problems)
{
   set_dcinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores given problems and warnings in the DCINFO block in the SLHA
 * object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void MRSSM2_slha_io::set_dcinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream dcinfo;
   dcinfo << "Block DCINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      dcinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      dcinfo << FORMAT_SPINFO(4, s);

   dcinfo << FORMAT_SPINFO(5, MRSSM2_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(dcinfo);
}

/**
 * Sort decays of every particle according to their width
 *
 */
std::vector<Decay> sort_decays_list(const Decays_list& decays_list) {
   std::vector<Decay> decays_list_as_vector;
   for (const auto& el : decays_list) {
      decays_list_as_vector.push_back(el.second);
   }
   std::sort(
      decays_list_as_vector.begin(),
      decays_list_as_vector.end(),
      [](const auto& d1, const auto& d2) {
         return d1.get_width() > d2.get_width();
      }
   );
   return decays_list_as_vector;
}

/**
 * Stores the branching ratios for a given particle in the SLHA
 * object.
 *
 * @param decays struct containing individual particle decays
 */
void MRSSM2_slha_io::set_decay_block(const Decays_list& decays_list, FlexibleDecay_settings const& flexibledecay_settings)
{
   const auto pdg = decays_list.get_particle_id();
   const auto width = decays_list.get_total_width();
   const std::string name = MRSSM2_info::get_particle_name_from_pdg(pdg);
   if (std::isnan(width) || std::isinf(width)) {
      throw std::runtime_error("Total width of " + name + " is " + std::to_string(width));
   }

   std::ostringstream decay;

   decay << "DECAY "
         << FORMAT_TOTAL_WIDTH(pdg, width, name + " decays");

   if (!is_zero(width, 1e-100)) {
      constexpr double NEGATIVE_BR_TOLERANCE = 1e-11;
      const double MIN_BR_TO_PRINT = flexibledecay_settings.get(FlexibleDecay_settings::min_br_to_print);
      std::vector<Decay> sorted_decays_list = sort_decays_list(decays_list);
      for (const auto& channel : sorted_decays_list) {
         auto const partial_width = channel.get_width();
         auto branching_ratio = partial_width / width;
         if (partial_width < 0 && !is_zero(branching_ratio, NEGATIVE_BR_TOLERANCE)) {
            std::stringstream ss;
            ss << std::scientific << partial_width;
            throw std::runtime_error("Error in " + channel.get_proc_string() + ": partial width is negative (" + ss.str() + " GeV).");
         }
         else if (partial_width < 0 && is_zero(branching_ratio, NEGATIVE_BR_TOLERANCE)) {
            branching_ratio = 0;
         }
         if (branching_ratio < MIN_BR_TO_PRINT) continue;
         const auto final_state = channel.get_final_state_particle_ids();
         std::string comment = "BR(" + name + " ->";
         for (auto id : final_state) {
            comment += " " + MRSSM2_info::get_particle_name_from_pdg(id);
         }
         comment += ")";

         decay << format_decay(branching_ratio, final_state, comment);
      }
   }

   slha_io.set_block(decay);
}

/**
 * Stores the particle decay branching ratios in the SLHA object.
 *
 * @param decays struct containing decays data
 */
void MRSSM2_slha_io::set_decays(const MRSSM2_decay_table& decay_table, FlexibleDecay_settings const& flexibledecay_settings)
{
   for (const auto& particle : decay_table) {
      set_decay_block(particle, flexibledecay_settings);
   }
}
void MRSSM2_slha_io::fill_decays_data(const MRSSM2_decays& decays, FlexibleDecay_settings const& flexibledecay_settings)
{
   const auto& decays_problems = decays.get_problems();
   const bool decays_error = decays_problems.have_problem();

   set_dcinfo(decays_problems);

   if (!decays_error) {
      set_decays(decays.get_decay_table(), flexibledecay_settings);
   }
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
void MRSSM2_slha_io::set_spectrum(const MRSSM2_slha& model)
{
   const MRSSM2_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
void MRSSM2_slha_io::set_extra(
   const MRSSM2_slha& model,
   const MRSSM2_scales& scales,
   const MRSSM2_observables& observables,
   const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings)
{
   const MRSSM2_physical physical(model.get_physical_slha());

   if (observables.problems.have_problem()) {
      std::ostringstream block;
      block << "Block OBSINFO\n";
      slha_format_problems_and_warnings(observables.problems,
                                        std::ostream_iterator<std::string>(block));
      slha_io.set_block(block);
   }

   if (spectrum_generator_settings.get(Spectrum_generator_settings::calculate_observables)) {
      {
         std::ostringstream block;
         block << "Block FlexibleSUSYLowEnergy Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
               << FORMAT_ELEMENT(21, (OBSERVABLES.a_muon), "Delta(g-2)_muon/2 FlexibleSUSY")
               << FORMAT_ELEMENT(23, (OBSERVABLES.edm_Fe_0), "electric dipole moment of Fe(0) [1/GeV]")
               << FORMAT_ELEMENT(24, (OBSERVABLES.edm_Fe_1), "electric dipole moment of Fe(1) [1/GeV]")
               << FORMAT_ELEMENT(25, (OBSERVABLES.edm_Fe_2), "electric dipole moment of Fe(2) [1/GeV]")
               << FORMAT_ELEMENT(26, (OBSERVABLES.Fe_to_Fe_VP), "BR(Fe1 -> Fe0 VP)")
         ;
         slha_io.set_block(block);
      }
   }

}

/**
 * Write SLHA object to given output.  If output == "-", then the SLHA
 * object is written to std::cout.  Otherwise, output is interpreted
 * as a file name
 *
 * @param output "-" for cout, or file name
 */
void MRSSM2_slha_io::write_to(const std::string& output) const
{
   if (output == "-")
      write_to_stream(std::cout);
   else
      write_to_file(output);
}

void MRSSM2_slha_io::write_to_file(const std::string& file_name) const
{
   slha_io.write_to_file(file_name);
}

void MRSSM2_slha_io::write_to_stream() const
{
   write_to_stream(std::cout);
}

void MRSSM2_slha_io::write_to_stream(std::ostream& ostr) const
{
   slha_io.write_to_stream(ostr);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double MRSSM2_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void MRSSM2_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
}

/**
 * Read SLHA object from source
 *
 * calls SLHA_io::read_from_source()
 *
 * @param source source name
 */
void MRSSM2_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void MRSSM2_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR,
 * EXTPAR and IMEXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void MRSSM2_slha_io::fill(MRSSM2_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor = [&input] (int key, double value) {
      return fill_minpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor extpar_processor = [&input] (int key, double value) {
      return fill_extpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imminpar_processor = [&input] (int key, double value) {
      return fill_imminpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imextpar_processor = [&input] (int key, double value) {
      return fill_imextpar_tuple(input, key, value);
   };

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);
   slha_io.read_block("IMMINPAR", imminpar_processor);
   slha_io.read_block("IMEXTPAR", imextpar_processor);

   input.BMuInput = slha_io.read_entry("HMIXIN", 101);
   input.LamSDInput = slha_io.read_entry("HMIXIN", 301);
   input.LamSUInput = slha_io.read_entry("HMIXIN", 302);
   input.LamTDInput = slha_io.read_entry("HMIXIN", 303);
   input.LamTUInput = slha_io.read_entry("HMIXIN", 304);
   slha_io.read_block("MSD2IN", input.md2Input);
   input.MDBSInput = slha_io.read_entry("MSOFTIN", 300);
   input.MDGocInput = slha_io.read_entry("MSOFTIN", 302);
   input.MDWBTInput = slha_io.read_entry("MSOFTIN", 301);
   slha_io.read_block("MSE2IN", input.me2Input);
   slha_io.read_block("MSL2IN", input.ml2Input);
   input.moc2Input = slha_io.read_entry("MSOFTIN", 111);
   slha_io.read_block("MSQ2IN", input.mq2Input);
   input.mRd2Input = slha_io.read_entry("MSOFTIN", 50);
   input.mRu2Input = slha_io.read_entry("MSOFTIN", 51);
   input.mS2Input = slha_io.read_entry("NMSSMRUNIN", 10);
   input.mT2Input = slha_io.read_entry("MSOFTIN", 110);
   slha_io.read_block("MSU2IN", input.mu2Input);
   input.MuDInput = slha_io.read_entry("HMIXIN", 201);
   input.MuUInput = slha_io.read_entry("HMIXIN", 202);

}

/**
 * Reads DR-bar parameters from a SLHA output file.
 *
 * @param model model class to be filled
 */
void MRSSM2_slha_io::fill_drbar_parameters(MRSSM2_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      Eigen::Matrix<double,3,3> Yu;
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      Eigen::Matrix<double,3,3> Yd;
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      Eigen::Matrix<double,3,3> Ye;
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_BMu(slha_io.read_entry("HMIX", 101));
   {
      Eigen::Matrix<double,3,3> mq2;
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      Eigen::Matrix<double,3,3> me2;
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      Eigen::Matrix<double,3,3> ml2;
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      Eigen::Matrix<double,3,3> mu2;
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      Eigen::Matrix<double,3,3> md2;
      slha_io.read_block("MSD2", md2);
      model.set_md2(md2);
   }
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));
   model.set_mS2(slha_io.read_entry("NMSSMRUN", 10));
   model.set_mT2(slha_io.read_entry("MSOFT", 110));
   model.set_moc2(slha_io.read_entry("MSOFT", 111));
   model.set_vS(slha_io.read_entry("NMSSMRUN", 5));
   model.set_vT(slha_io.read_entry("HMIX", 310));
   model.set_MDBS(slha_io.read_entry("MSOFT", 300));
   model.set_MDWBT(slha_io.read_entry("MSOFT", 301));
   model.set_MDGoc(slha_io.read_entry("MSOFT", 302));
   model.set_MuD(slha_io.read_entry("HMIX", 201));
   model.set_MuU(slha_io.read_entry("HMIX", 202));
   model.set_BMuD(slha_io.read_entry("HMIX", 203));
   model.set_BMuU(slha_io.read_entry("HMIX", 204));
   model.set_LamSD(slha_io.read_entry("HMIX", 301));
   model.set_LamSU(slha_io.read_entry("HMIX", 302));
   model.set_LamTD(slha_io.read_entry("HMIX", 303));
   model.set_LamTU(slha_io.read_entry("HMIX", 304));
   model.set_mRd2(slha_io.read_entry("MSOFT", 50));
   model.set_mRu2(slha_io.read_entry("MSOFT", 51));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 *
 * @param model model class to be filled
 */
void MRSSM2_slha_io::fill(MRSSM2_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   MRSSM2_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of extra physical input parameters from SLHA object
 * (FlexibleSUSYInput block)
 *
 * @param input struct of physical non-SLHA input parameters
 */
void MRSSM2_slha_io::fill(Physical_input& input) const
{
   slha_io.fill(input);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void MRSSM2_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void MRSSM2_slha_io::fill(FlexibleDecay_settings& settings) const
{
   slha_io.fill(settings);
}

void MRSSM2_slha_io::fill_minpar_tuple(MRSSM2_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 3: input.TanBeta = value; break;
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void MRSSM2_slha_io::fill_extpar_tuple(MRSSM2_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 0: input.Ms = value; break;
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

void MRSSM2_slha_io::fill_imminpar_tuple(MRSSM2_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMMINPAR: " << key); break;
   }

}

void MRSSM2_slha_io::fill_imextpar_tuple(MRSSM2_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMEXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file to be filled.
 */
void MRSSM2_slha_io::fill_physical(MRSSM2_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(ZD);
      slha_io.read_block("DSQMIX", ZD);
      LOCALPHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZU);
      slha_io.read_block("USQMIX", ZU);
      LOCALPHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZE);
      slha_io.read_block("SELMIX", ZE);
      LOCALPHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZV);
      slha_io.read_block("SNUMIX", ZV);
      LOCALPHYSICAL(ZV) = ZV;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZEL);
      slha_io.read_block("UELMIX", ZEL);
      LOCALPHYSICAL(ZEL) = ZEL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZER);
      slha_io.read_block("UERMIX", ZER);
      LOCALPHYSICAL(ZER) = ZER;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDL);
      slha_io.read_block("UDLMIX", ZDL);
      LOCALPHYSICAL(ZDL) = ZDL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZDR);
      slha_io.read_block("UDRMIX", ZDR);
      LOCALPHYSICAL(ZDR) = ZDR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZUL);
      slha_io.read_block("UULMIX", ZUL);
      LOCALPHYSICAL(ZUL) = ZUL;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZUR);
      slha_io.read_block("UURMIX", ZUR);
      LOCALPHYSICAL(ZUR) = ZUR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZHR);
      slha_io.read_block("RHMIX", ZHR);
      LOCALPHYSICAL(ZHR) = ZHR;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN1);
      slha_io.read_block("N1MIX", ZN1);
      LOCALPHYSICAL(ZN1) = ZN1;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN2);
      slha_io.read_block("N2MIX", ZN2);
      LOCALPHYSICAL(ZN2) = ZN2;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP1);
      slha_io.read_block("V1MIX", UP1);
      LOCALPHYSICAL(UP1) = UP1;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM1);
      slha_io.read_block("U1MIX", UM1);
      LOCALPHYSICAL(UM1) = UM1;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP2);
      slha_io.read_block("V2MIX", UP2);
      LOCALPHYSICAL(UP2) = UP2;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM2);
      slha_io.read_block("U2MIX", UM2);
      LOCALPHYSICAL(UM2) = UM2;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MSRdp) = slha_io.read_entry("MASS", 404);
   LOCALPHYSICAL(MSRum) = slha_io.read_entry("MASS", 403);
   LOCALPHYSICAL(MsigmaO) = slha_io.read_entry("MASS", 3000022);
   LOCALPHYSICAL(MphiO) = slha_io.read_entry("MASS", 3000021);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   LOCALPHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 1000003);
   LOCALPHYSICAL(MSd)(2) = slha_io.read_entry("MASS", 1000005);
   LOCALPHYSICAL(MSd)(3) = slha_io.read_entry("MASS", 2000001);
   LOCALPHYSICAL(MSd)(4) = slha_io.read_entry("MASS", 2000003);
   LOCALPHYSICAL(MSd)(5) = slha_io.read_entry("MASS", 2000005);
   LOCALPHYSICAL(MSv)(0) = slha_io.read_entry("MASS", 1000012);
   LOCALPHYSICAL(MSv)(1) = slha_io.read_entry("MASS", 1000014);
   LOCALPHYSICAL(MSv)(2) = slha_io.read_entry("MASS", 1000016);
   LOCALPHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   LOCALPHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 1000004);
   LOCALPHYSICAL(MSu)(2) = slha_io.read_entry("MASS", 1000006);
   LOCALPHYSICAL(MSu)(3) = slha_io.read_entry("MASS", 2000002);
   LOCALPHYSICAL(MSu)(4) = slha_io.read_entry("MASS", 2000004);
   LOCALPHYSICAL(MSu)(5) = slha_io.read_entry("MASS", 2000006);
   LOCALPHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   LOCALPHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 1000013);
   LOCALPHYSICAL(MSe)(2) = slha_io.read_entry("MASS", 1000015);
   LOCALPHYSICAL(MSe)(3) = slha_io.read_entry("MASS", 2000011);
   LOCALPHYSICAL(MSe)(4) = slha_io.read_entry("MASS", 2000013);
   LOCALPHYSICAL(MSe)(5) = slha_io.read_entry("MASS", 2000015);
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   LOCALPHYSICAL(Mhh)(3) = slha_io.read_entry("MASS", 55);
   LOCALPHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 46);
   LOCALPHYSICAL(MAh)(3) = slha_io.read_entry("MASS", 56);
   LOCALPHYSICAL(MRh)(0) = slha_io.read_entry("MASS", 401);
   LOCALPHYSICAL(MRh)(1) = slha_io.read_entry("MASS", 402);
   LOCALPHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MHpm)(2) = slha_io.read_entry("MASS", 47);
   LOCALPHYSICAL(MHpm)(3) = slha_io.read_entry("MASS", 57);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   LOCALPHYSICAL(MCha1)(0) = slha_io.read_entry("MASS", 1000024);
   LOCALPHYSICAL(MCha1)(1) = slha_io.read_entry("MASS", 1000037);
   LOCALPHYSICAL(MCha2)(0) = slha_io.read_entry("MASS", 2000024);
   LOCALPHYSICAL(MCha2)(1) = slha_io.read_entry("MASS", 2000037);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double MRSSM2_slha_io::read_scale() const
{
   static const std::array<std::string, 12> drbar_blocks =
      { "gauge", "Yu", "Yd", "Ye", "HMIX", "MSQ2", "MSE2", "MSL2", "MSU2", "MSD2",
   "MSOFT", "NMSSMRUN" };

   double scale = 0.;

   for (const auto& block: drbar_blocks) {
      const double block_scale = slha_io.read_scale(block);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
