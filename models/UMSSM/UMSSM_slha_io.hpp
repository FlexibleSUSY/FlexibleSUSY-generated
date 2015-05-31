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

// File generated at Sun 31 May 2015 12:31:25

#ifndef UMSSM_SLHA_IO_H
#define UMSSM_SLHA_IO_H

#include "UMSSM_two_scale_model_slha.hpp"
#include "UMSSM_info.hpp"
#include "UMSSM_physical.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"

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

struct UMSSM_input_parameters;
class Spectrum_generator_settings;

struct UMSSM_scales {
   UMSSM_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class UMSSM_slha_io {
public:
   UMSSM_slha_io();
   ~UMSSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(UMSSM_input_parameters&) const;
   void fill(UMSSM_mass_eigenstates&) const;
   template <class T> void fill(UMSSM_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const UMSSM_input_parameters&);
   template <class T> void set_extra(const UMSSM_slha<T>&, const UMSSM_scales&);
   void set_minpar(const UMSSM_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const UMSSM_slha<T>&);
   template <class T> void set_spectrum(const UMSSM<T>&);
   void set_spinfo(const Problems<UMSSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(UMSSM_input_parameters&, int, double);
   static void fill_extpar_tuple(UMSSM_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const UMSSM_slha<T>&, const QedQcd&, const UMSSM_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const UMSSM_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const UMSSM_slha<T>&, const QedQcd&, const UMSSM_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 16;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const UMSSM_physical&, bool);
   void set_mixing_matrices(const UMSSM_physical&, bool);
   template <class T> void set_model_parameters(const UMSSM_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(UMSSM_mass_eigenstates&) const;
   void fill_physical(UMSSM_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void UMSSM_slha_io::fill(UMSSM_slha<T>& model) const
{
   fill(static_cast<UMSSM_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class T>
void UMSSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const UMSSM_slha<T>& model,
   const QedQcd& qedqcd, const UMSSM_scales& scales)
{
   UMSSM_slha_io slha_io;
   const UMSSM_input_parameters& input = model.get_input();
   const Problems<UMSSM_info::NUMBER_OF_PARTICLES>& problems
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
SLHAea::Coll UMSSM_slha_io::fill_slhaea(
   const UMSSM_slha<T>& model, const QedQcd& qedqcd)
{
   UMSSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll UMSSM_slha_io::fill_slhaea(
   const UMSSM_slha<T>& model, const QedQcd& qedqcd,
   const UMSSM_scales& scales)
{
   SLHAea::Coll slhaea;
   UMSSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void UMSSM_slha_io::set_model_parameters(const UMSSM_slha<T>& model)
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
            << FORMAT_ELEMENT(4, (MODELPARAMETER(MassU)), "MassU")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block NMSSMRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(10, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(5, (MODELPARAMETER(vS)), "vS")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block GAUGE Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(5, (MODELPARAMETER(gp)), "gp")
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
void UMSSM_slha_io::set_extra(
   const UMSSM_slha<T>& model, const UMSSM_scales& scales)
{
   const UMSSM_physical physical(model.get_physical_slha());


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void UMSSM_slha_io::set_spectrum(const UMSSM<T>& model)
{
   const UMSSM_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void UMSSM_slha_io::set_spectrum(const UMSSM_slha<T>& model)
{
   const UMSSM_physical physical(model.get_physical_slha());
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
