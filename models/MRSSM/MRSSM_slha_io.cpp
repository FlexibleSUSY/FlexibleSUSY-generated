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

// File generated at Tue 5 Sep 2017 10:45:16

#include "MRSSM_slha_io.hpp"
#include "MRSSM_input_parameters.hpp"
#include "MRSSM_info.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p

using namespace softsusy;

namespace flexiblesusy {

char const * const MRSSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "HMIX", "MSQ2", "MSE2", "MSL2", "MSU2", "MSD2",
   "MSOFT", "NMSSMRUN" }
;

MRSSM_slha_io::MRSSM_slha_io()
   : slha_io()
   , print_imaginary_parts_of_majorana_mixings(false)
{
}

void MRSSM_slha_io::clear()
{
   slha_io.clear();
}

void MRSSM_slha_io::set_print_imaginary_parts_of_majorana_mixings(bool flag)
{
   print_imaginary_parts_of_majorana_mixings = flag;
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM_slha_io::set_extpar(const MRSSM_input_parameters& input)
{

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void MRSSM_slha_io::set_minpar(const MRSSM_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   slha_io.set_block(minpar);

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void MRSSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void MRSSM_slha_io::set_spinfo(const Problems<MRSSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::vector<std::string> warnings_vec, problems_vec;

   if (problems.have_warning()) {
      std::ostringstream ss;
      problems.print_warnings(ss);
      warnings_vec = { ss.str() };
   }

   if (problems.have_problem()) {
      std::ostringstream ss;
      problems.print_problems(ss);
      problems_vec = { ss.str() };
   }

   set_spinfo(problems_vec, warnings_vec);
}

/**
 * Stores the given problems and warnings in the SPINFO block in the
 * SLHA object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void MRSSM_slha_io::set_spinfo(
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

   spinfo << FORMAT_SPINFO(5, MRSSM_info::model_name)
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
void MRSSM_slha_io::set_mass(const MRSSM_physical& physical,
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
void MRSSM_slha_io::set_mixing_matrices(const MRSSM_physical& physical,
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

void MRSSM_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void MRSSM_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Write SLHA object to given output.  If output == "-", then the SLHA
 * object is written to std::cout.  Otherwise, output is interpreted
 * as a file name
 *
 * @param output "-" for cout, or file name
 */
void MRSSM_slha_io::write_to(const std::string& output)
{
   if (output == "-")
      write_to_stream(std::cout);
   else
      write_to_file(output);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double MRSSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void MRSSM_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
}

/**
 * Read SLHA object from source
 *
 * calls SLHA_io::read_from_source()
 *
 * @param source source name
 */
void MRSSM_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
   slha_io.read_modsel();
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void MRSSM_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void MRSSM_slha_io::fill(MRSSM_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(&MRSSM_slha_io::fill_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(&MRSSM_slha_io::fill_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   input.BMuDInput = slha_io.read_entry("HMIXIN", 203);
   input.BMuInput = slha_io.read_entry("HMIXIN", 101);
   input.BMuUInput = slha_io.read_entry("HMIXIN", 204);
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
   slha_io.read_block("MSU2IN", input.mu2Input);
   input.MuDInput = slha_io.read_entry("HMIXIN", 201);
   input.MuInput = slha_io.read_entry("HMIXIN", 1);
   input.MuUInput = slha_io.read_entry("HMIXIN", 202);
   input.vSInput = slha_io.read_entry("NMSSMRUNIN", 5);
   input.vTInput = slha_io.read_entry("HMIXIN", 310);

}

/**
 * Reads DR-bar parameters from a SLHA output file.
 *
 * @param model model class to be filled
 */
void MRSSM_slha_io::fill_drbar_parameters(MRSSM_mass_eigenstates& model) const
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
void MRSSM_slha_io::fill(MRSSM_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   MRSSM_physical physical_hk;
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
void MRSSM_slha_io::fill(Physical_input& input) const
{
   slha_io.fill(input);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void MRSSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

void MRSSM_slha_io::fill_minpar_tuple(MRSSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 3: input.TanBeta = value; break;
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void MRSSM_slha_io::fill_extpar_tuple(MRSSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file to be filled.
 */
void MRSSM_slha_io::fill_physical(MRSSM_physical& physical) const
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
double MRSSM_slha_io::read_scale() const
{
   double scale = 0.;

   for (unsigned i = 0; i < NUMBER_OF_DRBAR_BLOCKS; i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
