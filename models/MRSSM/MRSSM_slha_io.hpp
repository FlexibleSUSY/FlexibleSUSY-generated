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

// File generated at Tue 24 Feb 2015 17:31:34

#ifndef MRSSM_SLHA_IO_H
#define MRSSM_SLHA_IO_H

#include "MRSSM_two_scale_model_slha.hpp"
#include "MRSSM_info.hpp"
#include "MRSSM_physical.hpp"
#include "slha_io.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_POLE_MASS(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(PHYSICAL(p))>::type>::type p;
#define SM(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct MRSSM_input_parameters;
class Spectrum_generator_settings;

struct MRSSM_scales {
   MRSSM_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class MRSSM_slha_io {
public:
   MRSSM_slha_io();
   ~MRSSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(MRSSM_input_parameters&) const;
   template <class T> void fill(MRSSM_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const MRSSM_input_parameters&);
   template <class T> void set_extra(const MRSSM_slha<T>&, const MRSSM_scales&);
   void set_minpar(const MRSSM_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const MRSSM_slha<T>&);
   template <class T> void set_spectrum(const MRSSM<T>&);
   void set_spinfo(const Problems<MRSSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(MRSSM_input_parameters&, int, double);
   static void fill_extpar_tuple(MRSSM_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   static void convert_to_slha_convention(MRSSM_physical&);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const MRSSM_slha<T>&, const QedQcd&, const MRSSM_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const MRSSM_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const MRSSM_slha<T>&, const QedQcd&, const MRSSM_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 12;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const MRSSM_physical&, bool);
   void set_mixing_matrices(const MRSSM_physical&, bool);
   template <class T> void set_model_parameters(const MRSSM_slha<T>&);
   double read_scale() const;
   template <class T> void fill_drbar_parameters(MRSSM_slha<T>&) const;
   template <class T> void fill_physical(MRSSM_slha<T>&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void MRSSM_slha_io::fill(MRSSM_slha<T>& model) const
{
   fill_drbar_parameters(model);
   fill_physical(model);
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
template <class T>
void MRSSM_slha_io::fill_drbar_parameters(MRSSM_slha<T>& model) const
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
   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_BMu(slha_io.read_entry("HMIX", 101));
   {
      DEFINE_PARAMETER(mq2);
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      DEFINE_PARAMETER(me2);
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      DEFINE_PARAMETER(ml2);
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      DEFINE_PARAMETER(mu2);
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      DEFINE_PARAMETER(md2);
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
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
template <class T>
void MRSSM_slha_io::fill_physical(MRSSM_slha<T>& model) const
{
   {
      DEFINE_PARAMETER(ZD);
      slha_io.read_block("DSQMIX", ZD);
      PHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PARAMETER(ZU);
      slha_io.read_block("USQMIX", ZU);
      PHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PARAMETER(ZE);
      slha_io.read_block("SELMIX", ZE);
      PHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PARAMETER(ZV);
      slha_io.read_block("SNUMIX", ZV);
      PHYSICAL(ZV) = ZV;
   }
   {
      DEFINE_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      PHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      PHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      PHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PARAMETER(ZEL);
      slha_io.read_block("UELMIX", ZEL);
      PHYSICAL(ZEL) = ZEL;
   }
   {
      DEFINE_PARAMETER(ZER);
      slha_io.read_block("UERMIX", ZER);
      PHYSICAL(ZER) = ZER;
   }
   {
      DEFINE_PARAMETER(ZDL);
      slha_io.read_block("UDLMIX", ZDL);
      PHYSICAL(ZDL) = ZDL;
   }
   {
      DEFINE_PARAMETER(ZDR);
      slha_io.read_block("UDRMIX", ZDR);
      PHYSICAL(ZDR) = ZDR;
   }
   {
      DEFINE_PARAMETER(ZUL);
      slha_io.read_block("UULMIX", ZUL);
      PHYSICAL(ZUL) = ZUL;
   }
   {
      DEFINE_PARAMETER(ZUR);
      slha_io.read_block("UURMIX", ZUR);
      PHYSICAL(ZUR) = ZUR;
   }
   {
      DEFINE_PARAMETER(ZHR);
      slha_io.read_block("RHMIX", ZHR);
      PHYSICAL(ZHR) = ZHR;
   }
   {
      DEFINE_PARAMETER(ZRP);
      slha_io.read_block("RPMIX", ZRP);
      PHYSICAL(ZRP) = ZRP;
   }
   {
      DEFINE_PARAMETER(ZN1);
      slha_io.read_block("N1MIX", ZN1);
      PHYSICAL(ZN1) = ZN1;
   }
   {
      DEFINE_PARAMETER(ZN2);
      slha_io.read_block("N2MIX", ZN2);
      PHYSICAL(ZN2) = ZN2;
   }
   {
      DEFINE_PARAMETER(UP1);
      slha_io.read_block("V1MIX", UP1);
      PHYSICAL(UP1) = UP1;
   }
   {
      DEFINE_PARAMETER(UM1);
      slha_io.read_block("U1MIX", UM1);
      PHYSICAL(UM1) = UM1;
   }
   {
      DEFINE_PARAMETER(UP2);
      slha_io.read_block("V2MIX", UP2);
      PHYSICAL(UP2) = UP2;
   }
   {
      DEFINE_PARAMETER(UM2);
      slha_io.read_block("U2MIX", UM2);
      PHYSICAL(UM2) = UM2;
   }

   PHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   PHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   PHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   PHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   PHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   PHYSICAL(MSOc) = slha_io.read_entry("MASS", 3000021);
   PHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   PHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   PHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   PHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 1000003);
   PHYSICAL(MSd)(2) = slha_io.read_entry("MASS", 1000005);
   PHYSICAL(MSd)(3) = slha_io.read_entry("MASS", 2000001);
   PHYSICAL(MSd)(4) = slha_io.read_entry("MASS", 2000003);
   PHYSICAL(MSd)(5) = slha_io.read_entry("MASS", 2000005);
   PHYSICAL(MSv)(0) = slha_io.read_entry("MASS", 1000012);
   PHYSICAL(MSv)(1) = slha_io.read_entry("MASS", 1000014);
   PHYSICAL(MSv)(2) = slha_io.read_entry("MASS", 1000016);
   PHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   PHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 1000004);
   PHYSICAL(MSu)(2) = slha_io.read_entry("MASS", 1000006);
   PHYSICAL(MSu)(3) = slha_io.read_entry("MASS", 2000002);
   PHYSICAL(MSu)(4) = slha_io.read_entry("MASS", 2000004);
   PHYSICAL(MSu)(5) = slha_io.read_entry("MASS", 2000006);
   PHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   PHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 1000013);
   PHYSICAL(MSe)(2) = slha_io.read_entry("MASS", 1000015);
   PHYSICAL(MSe)(3) = slha_io.read_entry("MASS", 2000011);
   PHYSICAL(MSe)(4) = slha_io.read_entry("MASS", 2000013);
   PHYSICAL(MSe)(5) = slha_io.read_entry("MASS", 2000015);
   PHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   PHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   PHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   PHYSICAL(Mhh)(3) = slha_io.read_entry("MASS", 55);
   PHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   PHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 46);
   PHYSICAL(MAh)(3) = slha_io.read_entry("MASS", 56);
   PHYSICAL(MRh)(0) = slha_io.read_entry("MASS", 401);
   PHYSICAL(MRh)(1) = slha_io.read_entry("MASS", 402);
   PHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   PHYSICAL(MHpm)(2) = slha_io.read_entry("MASS", 47);
   PHYSICAL(MHpm)(3) = slha_io.read_entry("MASS", 57);
   PHYSICAL(MRpm)(0) = slha_io.read_entry("MASS", 410);
   PHYSICAL(MRpm)(1) = slha_io.read_entry("MASS", 411);
   PHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   PHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   PHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   PHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   PHYSICAL(MCha1)(0) = slha_io.read_entry("MASS", 1000024);
   PHYSICAL(MCha1)(1) = slha_io.read_entry("MASS", 1000037);
   PHYSICAL(MCha2)(0) = slha_io.read_entry("MASS", 2000024);
   PHYSICAL(MCha2)(1) = slha_io.read_entry("MASS", 2000037);
   PHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   PHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   PHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   PHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   PHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   PHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   PHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   PHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   PHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   PHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

template <class T>
void MRSSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const MRSSM_slha<T>& model,
   const QedQcd& qedqcd, const MRSSM_scales& scales)
{
   MRSSM_slha_io slha_io;
   const MRSSM_input_parameters& input = model.get_input();
   const Problems<MRSSM_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
SLHAea::Coll MRSSM_slha_io::fill_slhaea(
   const MRSSM_slha<T>& model, const QedQcd& qedqcd)
{
   MRSSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll MRSSM_slha_io::fill_slhaea(
   const MRSSM_slha<T>& model, const QedQcd& qedqcd,
   const MRSSM_scales& scales)
{
   SLHAea::Coll slhaea;
   MRSSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void MRSSM_slha_io::set_model_parameters(const MRSSM_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", MODELPARAMETER(Yu), "Yu", model.get_scale());
   slha_io.set_block("Yd", MODELPARAMETER(Yd), "Yd", model.get_scale());
   slha_io.set_block("Ye", MODELPARAMETER(Ye), "Ye", model.get_scale());
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
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2), "md2", model.get_scale());
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

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void MRSSM_slha_io::set_extra(
   const MRSSM_slha<T>& model, const MRSSM_scales& scales)
{
   const MRSSM_physical physical(model.get_physical_slha());


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void MRSSM_slha_io::set_spectrum(const MRSSM<T>& model)
{
   const MRSSM_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void MRSSM_slha_io::set_spectrum(const MRSSM_slha<T>& model)
{
   const MRSSM_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);
}

} // namespace flexiblesusy

#endif
