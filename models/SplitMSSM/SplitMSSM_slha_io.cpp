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

// File generated at Fri 8 Jan 2016 11:56:24

#include "SplitMSSM_slha_io.hpp"
#include "SplitMSSM_input_parameters.hpp"
#include "SplitMSSM_info.hpp"
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
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p

using namespace softsusy;

namespace flexiblesusy {

char const * const SplitMSSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "SM", "HMIX", "MSOFT", "SplitMSSM" }
;

SplitMSSM_slha_io::SplitMSSM_slha_io()
   : slha_io()
{
}

void SplitMSSM_slha_io::clear()
{
   slha_io.clear();
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void SplitMSSM_slha_io::set_extpar(const SplitMSSM_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(0, input.MSUSY, "MSUSY");
   extpar << FORMAT_ELEMENT(1, input.M1Input, "M1Input");
   extpar << FORMAT_ELEMENT(2, input.M2Input, "M2Input");
   extpar << FORMAT_ELEMENT(3, input.M3Input, "M3Input");
   extpar << FORMAT_ELEMENT(4, input.MuInput, "MuInput");
   extpar << FORMAT_ELEMENT(5, input.mAInput, "mAInput");
   extpar << FORMAT_ELEMENT(6, input.MEWSB, "MEWSB");
   extpar << FORMAT_ELEMENT(7, input.AtInput, "AtInput");
   extpar << FORMAT_ELEMENT(25, input.TanBeta, "TanBeta");
   slha_io.set_block(extpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void SplitMSSM_slha_io::set_minpar(const SplitMSSM_input_parameters& input)
{

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void SplitMSSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void SplitMSSM_slha_io::set_spinfo(const Problems<SplitMSSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_warning()) {
      std::ostringstream warnings;
      problems.print_warnings(warnings);
      spinfo << FORMAT_SPINFO(3, warnings.str());
   }

   if (problems.have_problem()) {
      std::ostringstream problems_str;
      problems.print_problems(problems_str);
      spinfo << FORMAT_SPINFO(4, problems_str.str());
   }

   spinfo << FORMAT_SPINFO(5, SplitMSSM_info::model_name)
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
void SplitMSSM_slha_io::set_mass(const SplitMSSM_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWp), "VWp")
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha(2)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(25, LOCALPHYSICAL(Mhh), "hh")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
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
void SplitMSSM_slha_io::set_mixing_matrices(const SplitMSSM_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("NMIX", LOCALPHYSICAL(ZN), "ZN");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(Ue), "Ue");
   }

}

void SplitMSSM_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void SplitMSSM_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Write SLHA object to file.
 *
 * @param file_name file name
 */
void SplitMSSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double SplitMSSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void SplitMSSM_slha_io::read_from_file(const std::string& file_name)
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
void SplitMSSM_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
   slha_io.read_modsel();
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void SplitMSSM_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void SplitMSSM_slha_io::fill(SplitMSSM_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(&SplitMSSM_slha_io::fill_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(&SplitMSSM_slha_io::fill_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   slha_io.read_block("MSQ2IN", input.msq2);
   slha_io.read_block("MSU2IN", input.msu2);
   slha_io.read_block("MSD2IN", input.msd2);
   slha_io.read_block("MSL2IN", input.msl2);
   slha_io.read_block("MSE2IN", input.mse2);

}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
void SplitMSSM_slha_io::fill_drbar_parameters(SplitMSSM_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_mu2(slha_io.read_entry("SM", 1));
   model.set_Lambdax(slha_io.read_entry("SM", 2));
   model.set_v(slha_io.read_entry("HMIX", 3));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_gYu(slha_io.read_entry("SplitMSSM", 1));
   model.set_g2u(slha_io.read_entry("SplitMSSM", 2));
   model.set_gYd(slha_io.read_entry("SplitMSSM", 3));
   model.set_g2d(slha_io.read_entry("SplitMSSM", 4));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 */
void SplitMSSM_slha_io::fill(SplitMSSM_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   SplitMSSM_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void SplitMSSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

void SplitMSSM_slha_io::fill_minpar_tuple(SplitMSSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void SplitMSSM_slha_io::fill_extpar_tuple(SplitMSSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 0: input.MSUSY = value; break;
   case 1: input.M1Input = value; break;
   case 2: input.M2Input = value; break;
   case 3: input.M3Input = value; break;
   case 4: input.MuInput = value; break;
   case 5: input.mAInput = value; break;
   case 6: input.MEWSB = value; break;
   case 7: input.AtInput = value; break;
   case 25: input.TanBeta = value; break;
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
void SplitMSSM_slha_io::fill_physical(SplitMSSM_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(Vu);
      slha_io.read_block("UULMIX", Vu);
      LOCALPHYSICAL(Vu) = Vu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vd);
      slha_io.read_block("UDLMIX", Vd);
      LOCALPHYSICAL(Vd) = Vd;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Uu);
      slha_io.read_block("UURMIX", Uu);
      LOCALPHYSICAL(Uu) = Uu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ud);
      slha_io.read_block("UDRMIX", Ud);
      LOCALPHYSICAL(Ud) = Ud;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ve);
      slha_io.read_block("UELMIX", Ve);
      LOCALPHYSICAL(Ve) = Ve;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ue);
      slha_io.read_block("UERMIX", Ue);
      LOCALPHYSICAL(Ue) = Ue;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN);
      slha_io.read_block("NMIX", ZN);
      LOCALPHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      LOCALPHYSICAL(UP) = UP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      LOCALPHYSICAL(UM) = UM;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   LOCALPHYSICAL(Mhh) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   LOCALPHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   LOCALPHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   LOCALPHYSICAL(MVWp) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double SplitMSSM_slha_io::read_scale() const
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
