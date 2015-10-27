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

// File generated at Tue 27 Oct 2015 15:30:07

#ifndef MSSMNoFVatMGUT_SLHA_IO_H
#define MSSMNoFVatMGUT_SLHA_IO_H

#include "MSSMNoFVatMGUT_two_scale_model_slha.hpp"
#include "MSSMNoFVatMGUT_info.hpp"
#include "MSSMNoFVatMGUT_physical.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct MSSMNoFVatMGUT_input_parameters;
class Spectrum_generator_settings;

struct MSSMNoFVatMGUT_scales {
   MSSMNoFVatMGUT_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class MSSMNoFVatMGUT_slha_io {
public:
   MSSMNoFVatMGUT_slha_io();
   ~MSSMNoFVatMGUT_slha_io() {}

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(MSSMNoFVatMGUT_input_parameters&) const;
   void fill(MSSMNoFVatMGUT_mass_eigenstates&) const;
   template <class T> void fill(MSSMNoFVatMGUT_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_extpar(const MSSMNoFVatMGUT_input_parameters&);
   template <class T> void set_extra(const MSSMNoFVatMGUT_slha<T>&, const MSSMNoFVatMGUT_scales&);
   void set_minpar(const MSSMNoFVatMGUT_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const MSSMNoFVatMGUT_slha<T>&);
   template <class T> void set_spectrum(const MSSMNoFVatMGUT<T>&);
   void set_spinfo(const Problems<MSSMNoFVatMGUT_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(MSSMNoFVatMGUT_input_parameters&, int, double);
   static void fill_extpar_tuple(MSSMNoFVatMGUT_input_parameters&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const MSSMNoFVatMGUT_slha<T>&, const softsusy::QedQcd&, const MSSMNoFVatMGUT_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const MSSMNoFVatMGUT_slha<T>&, const softsusy::QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const MSSMNoFVatMGUT_slha<T>&, const softsusy::QedQcd&, const MSSMNoFVatMGUT_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 14;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const MSSMNoFVatMGUT_physical&, bool);
   void set_mixing_matrices(const MSSMNoFVatMGUT_physical&, bool);
   template <class T> void set_model_parameters(const MSSMNoFVatMGUT_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(MSSMNoFVatMGUT_mass_eigenstates&) const;
   void fill_physical(MSSMNoFVatMGUT_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void MSSMNoFVatMGUT_slha_io::fill(MSSMNoFVatMGUT_slha<T>& model) const
{
   fill(static_cast<MSSMNoFVatMGUT_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class T>
void MSSMNoFVatMGUT_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const MSSMNoFVatMGUT_slha<T>& model,
   const softsusy::QedQcd& qedqcd, const MSSMNoFVatMGUT_scales& scales)
{
   MSSMNoFVatMGUT_slha_io slha_io;
   const MSSMNoFVatMGUT_input_parameters& input = model.get_input();
   const Problems<MSSMNoFVatMGUT_info::NUMBER_OF_PARTICLES>& problems
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
SLHAea::Coll MSSMNoFVatMGUT_slha_io::fill_slhaea(
   const MSSMNoFVatMGUT_slha<T>& model, const softsusy::QedQcd& qedqcd)
{
   MSSMNoFVatMGUT_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll MSSMNoFVatMGUT_slha_io::fill_slhaea(
   const MSSMNoFVatMGUT_slha<T>& model, const softsusy::QedQcd& qedqcd,
   const MSSMNoFVatMGUT_scales& scales)
{
   SLHAea::Coll slhaea;
   MSSMNoFVatMGUT_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void MSSMNoFVatMGUT_slha_io::set_model_parameters(const MSSMNoFVatMGUT_slha<T>& model)
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
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block("Te", MODELPARAMETER(TYe_slha), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd_slha), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu_slha), "TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
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
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
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
void MSSMNoFVatMGUT_slha_io::set_extra(
   const MSSMNoFVatMGUT_slha<T>& model, const MSSMNoFVatMGUT_scales& scales)
{
   const MSSMNoFVatMGUT_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(1, (SCALES(SUSYScale)), "SUSYScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block ALPHA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_NUMBER((ArcSin(Pole(ZH(1,1)))), "ArcSin(Pole(ZH(2,2)))")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(vu)/MODELPARAMETER(vd)), "vu/vd")
            << FORMAT_ELEMENT(3, (Sqrt(Sqr(MODELPARAMETER(vd)) + Sqr(MODELPARAMETER(vu)))), "Sqrt(Sqr(vd) + Sqr(vu))")
            << FORMAT_ELEMENT(4, (Sqr(MODELPARAMETER(MAh)(1))), "Sqr(MAh(2))")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMu)), "BMu")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Au Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYu)(0,0)/MODELPARAMETER(Yu)(0,0)), "TYu(1,1)/Yu(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYu)(1,1)/MODELPARAMETER(Yu)(1,1)), "TYu(2,2)/Yu(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYu)(2,2)/MODELPARAMETER(Yu)(2,2)), "TYu(3,3)/Yu(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ad Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYd)(0,0)/MODELPARAMETER(Yd)(0,0)), "TYd(1,1)/Yd(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYd)(1,1)/MODELPARAMETER(Yd)(1,1)), "TYd(2,2)/Yd(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYd)(2,2)/MODELPARAMETER(Yd)(2,2)), "TYd(3,3)/Yd(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block Ae Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_MIXING_MATRIX(1, 1, (MODELPARAMETER(TYe)(0,0)/MODELPARAMETER(Ye)(0,0)), "TYe(1,1)/Ye(1,1)")
            << FORMAT_MIXING_MATRIX(2, 2, (MODELPARAMETER(TYe)(1,1)/MODELPARAMETER(Ye)(1,1)), "TYe(2,2)/Ye(2,2)")
            << FORMAT_MIXING_MATRIX(3, 3, (MODELPARAMETER(TYe)(2,2)/MODELPARAMETER(Ye)(2,2)), "TYe(3,3)/Ye(3,3)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(31, (SignedAbsSqrt(MODELPARAMETER(ml2)(0,0))), "SignedAbsSqrt(ml2(1,1))")
            << FORMAT_ELEMENT(32, (SignedAbsSqrt(MODELPARAMETER(ml2)(1,1))), "SignedAbsSqrt(ml2(2,2))")
            << FORMAT_ELEMENT(33, (SignedAbsSqrt(MODELPARAMETER(ml2)(2,2))), "SignedAbsSqrt(ml2(3,3))")
            << FORMAT_ELEMENT(34, (SignedAbsSqrt(MODELPARAMETER(me2)(0,0))), "SignedAbsSqrt(me2(1,1))")
            << FORMAT_ELEMENT(35, (SignedAbsSqrt(MODELPARAMETER(me2)(1,1))), "SignedAbsSqrt(me2(2,2))")
            << FORMAT_ELEMENT(36, (SignedAbsSqrt(MODELPARAMETER(me2)(2,2))), "SignedAbsSqrt(me2(3,3))")
            << FORMAT_ELEMENT(41, (SignedAbsSqrt(MODELPARAMETER(mq2)(0,0))), "SignedAbsSqrt(mq2(1,1))")
            << FORMAT_ELEMENT(42, (SignedAbsSqrt(MODELPARAMETER(mq2)(1,1))), "SignedAbsSqrt(mq2(2,2))")
            << FORMAT_ELEMENT(43, (SignedAbsSqrt(MODELPARAMETER(mq2)(2,2))), "SignedAbsSqrt(mq2(3,3))")
            << FORMAT_ELEMENT(44, (SignedAbsSqrt(MODELPARAMETER(mu2)(0,0))), "SignedAbsSqrt(mu2(1,1))")
            << FORMAT_ELEMENT(45, (SignedAbsSqrt(MODELPARAMETER(mu2)(1,1))), "SignedAbsSqrt(mu2(2,2))")
            << FORMAT_ELEMENT(46, (SignedAbsSqrt(MODELPARAMETER(mu2)(2,2))), "SignedAbsSqrt(mu2(3,3))")
            << FORMAT_ELEMENT(47, (SignedAbsSqrt(MODELPARAMETER(md2)(0,0))), "SignedAbsSqrt(md2(1,1))")
            << FORMAT_ELEMENT(48, (SignedAbsSqrt(MODELPARAMETER(md2)(1,1))), "SignedAbsSqrt(md2(2,2))")
            << FORMAT_ELEMENT(49, (SignedAbsSqrt(MODELPARAMETER(md2)(2,2))), "SignedAbsSqrt(md2(3,3))")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block MASS Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1000021, (Pole(MGlu)), "Pole(MGlu)")
            << FORMAT_ELEMENT(1000012, (Pole(MSveL)), "Pole(MSveL)")
            << FORMAT_ELEMENT(1000014, (Pole(MSvmL)), "Pole(MSvmL)")
            << FORMAT_ELEMENT(1000016, (Pole(MSvtL)), "Pole(MSvtL)")
            << FORMAT_ELEMENT(1000024, (Pole(MCha(0))), "Pole(MCha(1))")
            << FORMAT_ELEMENT(1000037, (Pole(MCha(1))), "Pole(MCha(2))")
            << FORMAT_ELEMENT(25, (Pole(Mhh(0))), "Pole(Mhh(1))")
            << FORMAT_ELEMENT(35, (Pole(Mhh(1))), "Pole(Mhh(2))")
            << FORMAT_ELEMENT(37, (Pole(MHpm(1))), "Pole(MHpm(2))")
            << FORMAT_ELEMENT(36, (Pole(MAh(1))), "Pole(MAh(2))")
            << FORMAT_ELEMENT(1000001, (AbsSqr(Pole(ZD(0,0)))*Pole(MSd(0)) + AbsSqr(Pole(ZD(1,0)))*Pole(MSd(1))), "AbsSqr(Pole(ZD(1,1)))*Pole(MSd(1)) + AbsSqr(Pole(ZD(2,1)))*Pole(MSd(2))")
            << FORMAT_ELEMENT(2000001, (AbsSqr(Pole(ZD(0,1)))*Pole(MSd(0)) + AbsSqr(Pole(ZD(1,1)))*Pole(MSd(1))), "AbsSqr(Pole(ZD(1,2)))*Pole(MSd(1)) + AbsSqr(Pole(ZD(2,2)))*Pole(MSd(2))")
            << FORMAT_ELEMENT(1000003, (AbsSqr(Pole(ZS(0,0)))*Pole(MSs(0)) + AbsSqr(Pole(ZS(1,0)))*Pole(MSs(1))), "AbsSqr(Pole(ZS(1,1)))*Pole(MSs(1)) + AbsSqr(Pole(ZS(2,1)))*Pole(MSs(2))")
            << FORMAT_ELEMENT(2000003, (AbsSqr(Pole(ZS(0,1)))*Pole(MSs(0)) + AbsSqr(Pole(ZS(1,1)))*Pole(MSs(1))), "AbsSqr(Pole(ZS(1,2)))*Pole(MSs(1)) + AbsSqr(Pole(ZS(2,2)))*Pole(MSs(2))")
            << FORMAT_ELEMENT(1000005, (Pole(MSb(0))), "Pole(MSb(1))")
            << FORMAT_ELEMENT(2000005, (Pole(MSb(1))), "Pole(MSb(2))")
            << FORMAT_ELEMENT(1000011, (AbsSqr(Pole(ZE(0,0)))*Pole(MSe(0)) + AbsSqr(Pole(ZE(1,0)))*Pole(MSe(1))), "AbsSqr(Pole(ZE(1,1)))*Pole(MSe(1)) + AbsSqr(Pole(ZE(2,1)))*Pole(MSe(2))")
            << FORMAT_ELEMENT(2000011, (AbsSqr(Pole(ZE(0,1)))*Pole(MSe(0)) + AbsSqr(Pole(ZE(1,1)))*Pole(MSe(1))), "AbsSqr(Pole(ZE(1,2)))*Pole(MSe(1)) + AbsSqr(Pole(ZE(2,2)))*Pole(MSe(2))")
            << FORMAT_ELEMENT(1000013, (AbsSqr(Pole(ZM(0,0)))*Pole(MSm(0)) + AbsSqr(Pole(ZM(1,0)))*Pole(MSm(1))), "AbsSqr(Pole(ZM(1,1)))*Pole(MSm(1)) + AbsSqr(Pole(ZM(2,1)))*Pole(MSm(2))")
            << FORMAT_ELEMENT(2000013, (AbsSqr(Pole(ZM(0,1)))*Pole(MSm(0)) + AbsSqr(Pole(ZM(1,1)))*Pole(MSm(1))), "AbsSqr(Pole(ZM(1,2)))*Pole(MSm(1)) + AbsSqr(Pole(ZM(2,2)))*Pole(MSm(2))")
            << FORMAT_ELEMENT(1000015, (Pole(MStau(0))), "Pole(MStau(1))")
            << FORMAT_ELEMENT(2000015, (Pole(MStau(1))), "Pole(MStau(2))")
            << FORMAT_ELEMENT(1000002, (AbsSqr(Pole(ZU(0,0)))*Pole(MSu(0)) + AbsSqr(Pole(ZU(1,0)))*Pole(MSu(1))), "AbsSqr(Pole(ZU(1,1)))*Pole(MSu(1)) + AbsSqr(Pole(ZU(2,1)))*Pole(MSu(2))")
            << FORMAT_ELEMENT(2000002, (AbsSqr(Pole(ZU(0,1)))*Pole(MSu(0)) + AbsSqr(Pole(ZU(1,1)))*Pole(MSu(1))), "AbsSqr(Pole(ZU(1,2)))*Pole(MSu(1)) + AbsSqr(Pole(ZU(2,2)))*Pole(MSu(2))")
            << FORMAT_ELEMENT(1000004, (AbsSqr(Pole(ZC(0,0)))*Pole(MSc(0)) + AbsSqr(Pole(ZC(1,0)))*Pole(MSc(1))), "AbsSqr(Pole(ZC(1,1)))*Pole(MSc(1)) + AbsSqr(Pole(ZC(2,1)))*Pole(MSc(2))")
            << FORMAT_ELEMENT(2000004, (AbsSqr(Pole(ZC(0,1)))*Pole(MSc(0)) + AbsSqr(Pole(ZC(1,1)))*Pole(MSc(1))), "AbsSqr(Pole(ZC(1,2)))*Pole(MSc(1)) + AbsSqr(Pole(ZC(2,2)))*Pole(MSc(2))")
            << FORMAT_ELEMENT(1000006, (Pole(MSt(0))), "Pole(MSt(1))")
            << FORMAT_ELEMENT(2000006, (Pole(MSt(1))), "Pole(MSt(2))")
            << FORMAT_ELEMENT(1000022, (Pole(MChi(0))), "Pole(MChi(1))")
            << FORMAT_ELEMENT(1000023, (Pole(MChi(1))), "Pole(MChi(2))")
            << FORMAT_ELEMENT(1000025, (Pole(MChi(2))), "Pole(MChi(3))")
            << FORMAT_ELEMENT(1000035, (Pole(MChi(3))), "Pole(MChi(4))")
            << FORMAT_ELEMENT(24, (Pole(MVWm)), "Pole(MVWm)")
            << FORMAT_ELEMENT(23, (Pole(MVZ)), "Pole(MVZ)")
            << FORMAT_ELEMENT(1, (Pole(MFd)), "Pole(MFd)")
            << FORMAT_ELEMENT(2, (Pole(MFu)), "Pole(MFu)")
            << FORMAT_ELEMENT(3, (Pole(MFs)), "Pole(MFs)")
            << FORMAT_ELEMENT(4, (Pole(MFc)), "Pole(MFc)")
            << FORMAT_ELEMENT(5, (Pole(MFb)), "Pole(MFb)")
            << FORMAT_ELEMENT(6, (Pole(MFt)), "Pole(MFt)")
            << FORMAT_ELEMENT(11, (Pole(MFe)), "Pole(MFe)")
            << FORMAT_ELEMENT(13, (Pole(MFm)), "Pole(MFm)")
            << FORMAT_ELEMENT(15, (Pole(MFtau)), "Pole(MFtau)")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void MSSMNoFVatMGUT_slha_io::set_spectrum(const MSSMNoFVatMGUT<T>& model)
{
   const MSSMNoFVatMGUT_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void MSSMNoFVatMGUT_slha_io::set_spectrum(const MSSMNoFVatMGUT_slha<T>& model)
{
   const MSSMNoFVatMGUT_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#undef Pole
#undef PHYSICAL
#undef PHYSICAL_SLHA
#undef LOCALPHYSICAL
#undef MODELPARAMETER
#undef LowEnergyConstant
#undef SCALES

#endif
