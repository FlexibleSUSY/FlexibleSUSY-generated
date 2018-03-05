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

// File generated at Mon 5 Mar 2018 18:37:29

#ifndef NUTSMSSM_SLHA_IO_H
#define NUTSMSSM_SLHA_IO_H

#include "NUTSMSSM_mass_eigenstates.hpp"
#include "NUTSMSSM_model_slha.hpp"
#include "NUTSMSSM_info.hpp"
#include "NUTSMSSM_observables.hpp"
#include "NUTSMSSM_physical.hpp"
#include "problems.hpp"
#include "spectrum_generator_problems.hpp"
#include "standard_model_two_scale_model.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <tuple>
#include <utility>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define EXTRAPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct NUTSMSSM_input_parameters;
class Spectrum_generator_settings;

template <class T>
class NUTSMSSM;

struct NUTSMSSM_scales {
   double HighScale{0.}, SUSYScale{0.}, LowScale{0.};
   double pole_mass_scale{0.};
};

class NUTSMSSM_slha_io {
public:
   NUTSMSSM_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(NUTSMSSM_input_parameters&) const;
   void fill(NUTSMSSM_mass_eigenstates&) const;
   template <class Model> void fill(NUTSMSSM_slha<Model>&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   template <class Model> void set_extra(const NUTSMSSM_slha<Model>&, const NUTSMSSM_scales&, const NUTSMSSM_observables&);
   void set_input(const NUTSMSSM_input_parameters&);
   void set_modsel(const SLHA_io::Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class... Ts> void set_spectrum(const std::tuple<Ts...>&);
   template <class Model> void set_spectrum(const NUTSMSSM_slha<Model>&);
   template <class T> void set_spectrum(const NUTSMSSM<T>&);
   void set_spectrum(const standard_model::Standard_model&);
   void set_spinfo(const Spectrum_generator_problems&);
   void set_spinfo(const Problems&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&) const;
   void write_to_file(const std::string& file_name) const { slha_io.write_to_file(file_name); }
   void write_to_stream(std::ostream& ostr = std::cout) const { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(NUTSMSSM_input_parameters&, int, double);
   static void fill_extpar_tuple(NUTSMSSM_input_parameters&, int, double);
   static void fill_imminpar_tuple(NUTSMSSM_input_parameters&, int, double);
   static void fill_imextpar_tuple(NUTSMSSM_input_parameters&, int, double);

   template <class Model>
   static void fill_slhaea(SLHAea::Coll&, const NUTSMSSM_slha<Model>&, const softsusy::QedQcd&, const NUTSMSSM_scales&, const NUTSMSSM_observables&);

   template <class Model>
   static SLHAea::Coll fill_slhaea(const NUTSMSSM_slha<Model>&, const softsusy::QedQcd&, const NUTSMSSM_scales&, const NUTSMSSM_observables&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;

   void set_extpar(const NUTSMSSM_input_parameters&);
   void set_imminpar(const NUTSMSSM_input_parameters&);
   void set_imextpar(const NUTSMSSM_input_parameters&);
   void set_minpar(const NUTSMSSM_input_parameters&);
   void set_mass(const NUTSMSSM_physical&, bool);
   void set_mass(const standard_model::Standard_model_physical&);
   void set_mixing_matrices(const NUTSMSSM_physical&, bool);
   void set_mixing_matrices(const standard_model::Standard_model_physical&);
   template <class Model> void set_model_parameters(const NUTSMSSM_slha<Model>&);
   void set_model_parameters(const standard_model::Standard_model&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(NUTSMSSM_mass_eigenstates&) const;
   void fill_physical(NUTSMSSM_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class Model>
void NUTSMSSM_slha_io::fill(NUTSMSSM_slha<Model>& model) const
{
   fill(static_cast<NUTSMSSM_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class Model>
void NUTSMSSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const NUTSMSSM_slha<Model>& model,
   const softsusy::QedQcd& qedqcd, const NUTSMSSM_scales& scales,
   const NUTSMSSM_observables& observables)
{
   NUTSMSSM_slha_io slha_io;
   const NUTSMSSM_input_parameters& input = model.get_input();
   const auto& problems = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_input(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales, observables);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class Model>
SLHAea::Coll NUTSMSSM_slha_io::fill_slhaea(
   const NUTSMSSM_slha<Model>& model, const softsusy::QedQcd& qedqcd,
   const NUTSMSSM_scales& scales, const NUTSMSSM_observables& observables)
{
   SLHAea::Coll slhaea;
   NUTSMSSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales, observables);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class Model>
void NUTSMSSM_slha_io::set_model_parameters(const NUTSMSSM_slha<Model>& model)
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
            << FORMAT_ELEMENT(30, (MODELPARAMETER(L1)), "L1")
            << FORMAT_ELEMENT(31, (MODELPARAMETER(MS)), "MS")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(LL1)), "LL1")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(BMS)), "BMS")
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
   {
      std::ostringstream block;
      block << "Block NMSSMRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Kappa)), "Kappa")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(TKappa)), "TKappa")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(10, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(5, (MODELPARAMETER(vS)), "vS")
      ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block Phases Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (Re(MODELPARAMETER(PhaseGlu))), "Re(PhaseGlu)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block IMPhases Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (Im(MODELPARAMETER(PhaseGlu))), "Im(PhaseGlu)")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
template <class Model>
void NUTSMSSM_slha_io::set_extra(
   const NUTSMSSM_slha<Model>& model, const NUTSMSSM_scales& scales,
   const NUTSMSSM_observables& observables)
{
   const NUTSMSSM_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput" << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(1, (SCALES(SUSYScale)), "SUSYScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block NMSSMRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Kappa)), "Kappa")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(TLambdax)/MODELPARAMETER(Lambdax)), "TLambdax/Lambdax")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(TKappa)/MODELPARAMETER(Kappa)), "TKappa/Kappa")
            << FORMAT_ELEMENT(5, (0.7071067811865475*MODELPARAMETER(vS)*MODELPARAMETER(Lambdax)), "0.7071067811865475*vS*Lambdax")
            << FORMAT_ELEMENT(6, (MODELPARAMETER(L1)), "L1")
            << FORMAT_ELEMENT(7, (MODELPARAMETER(LL1)), "LL1")
            << FORMAT_ELEMENT(8, (MODELPARAMETER(MS)), "MS")
            << FORMAT_ELEMENT(9, (MODELPARAMETER(BMS)), "BMS")
            << FORMAT_ELEMENT(10, (MODELPARAMETER(ms2)), "ms2")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices of
 * all given models in the SLHA object.
 *
 * @todo Use generic lambda instead of Set_spectrum in C++14
 *
 * @param models model classes
 */
template <class... Ts>
void NUTSMSSM_slha_io::set_spectrum(const std::tuple<Ts...>& models)
{
   Set_spectrum<NUTSMSSM_slha_io> ss(this);
   boost::fusion::for_each(models, ss);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void NUTSMSSM_slha_io::set_spectrum(const NUTSMSSM<T>& model)
{
   set_spectrum(NUTSMSSM_slha<NUTSMSSM<T> >(model));
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class Model>
void NUTSMSSM_slha_io::set_spectrum(const NUTSMSSM_slha<Model>& model)
{
   const NUTSMSSM_physical physical(model.get_physical_slha());
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
#undef MODEL
#undef MODELPARAMETER
#undef EXTRAPARAMETER
#undef OBSERVABLES
#undef LowEnergyConstant
#undef SCALES

#endif
