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

#include "MSSMNoFVatMGUT_model_slha.hpp"
#include "ckm.hpp"
#include "linalg2.hpp"
#include "mixings.hpp"
#include "pmns.hpp"
#include "wrappers.hpp"

#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) this->p
#define PHYSICAL_SLHA(p) physical_slha.p
#define PHYSICAL_SLHA_REAL(p) Re(physical_slha.p)

namespace flexiblesusy {

MSSMNoFVatMGUT_slha::MSSMNoFVatMGUT_slha(const MSSMNoFVatMGUT_input_parameters& input_,
                                   bool do_convert_masses_to_slha)
   : MSSMNoFVatMGUT_mass_eigenstates(input_)
   , convert_masses_to_slha(do_convert_masses_to_slha)
{
}

/**
 * Copy constructor.  Copies from base class (model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 * @param do_convert_masses_to_slha whether to convert majorana
 *    fermion masses to SLHA convention (allow them to be negative)
 */
MSSMNoFVatMGUT_slha::MSSMNoFVatMGUT_slha(const MSSMNoFVatMGUT_mass_eigenstates& model_, bool do_convert_masses_to_slha)
   : MSSMNoFVatMGUT_mass_eigenstates(model_)
   , convert_masses_to_slha(do_convert_masses_to_slha)
{
   convert_to_slha();
}

void MSSMNoFVatMGUT_slha::clear()
{
   MSSMNoFVatMGUT_mass_eigenstates::clear();
   physical_slha.clear();
}

void MSSMNoFVatMGUT_slha::calculate_spectrum()
{
   MSSMNoFVatMGUT_mass_eigenstates::calculate_spectrum();
   convert_to_slha();
}

void MSSMNoFVatMGUT_slha::convert_to_slha()
{
   physical_slha = this->get_physical();

   if (convert_masses_to_slha)
      physical_slha.convert_to_slha();

   convert_yukawa_couplings_to_slha();
   calculate_ckm_matrix();
   calculate_pmns_matrix();
   convert_trilinear_couplings_to_slha();
   convert_soft_squared_masses_to_slha();
}

void MSSMNoFVatMGUT_slha::calculate_ckm_matrix()
{

}

void MSSMNoFVatMGUT_slha::calculate_pmns_matrix()
{
   pmns << 1, 0, 0, 0, 1, 0, 0, 0, 1;

}

/**
 * Convert Yukawa couplings to SLHA convention
 */
void MSSMNoFVatMGUT_slha::convert_yukawa_couplings_to_slha()
{
   Yu_slha = MODELPARAMETER(Yu).diagonal().real();
   Yd_slha = MODELPARAMETER(Yd).diagonal().real();
   Ye_slha = MODELPARAMETER(Ye).diagonal().real();

}

/**
 * Convert trilinear couplings to SLHA convention
 */
void MSSMNoFVatMGUT_slha::convert_trilinear_couplings_to_slha()
{
   TYu_slha = (MODELPARAMETER(TYu)).real();
   TYd_slha = (MODELPARAMETER(TYd)).real();
   TYe_slha = (MODELPARAMETER(TYe)).real();

}

/**
 * Convert soft-breaking squared mass parameters to SLHA convention
 */
void MSSMNoFVatMGUT_slha::convert_soft_squared_masses_to_slha()
{
   mq2_slha = (MODELPARAMETER(mq2)).real();
   mu2_slha = (MODELPARAMETER(mu2)).real();
   md2_slha = (MODELPARAMETER(md2)).real();
   ml2_slha = (MODELPARAMETER(ml2)).real();
   me2_slha = (MODELPARAMETER(me2)).real();

}

const MSSMNoFVatMGUT_physical& MSSMNoFVatMGUT_slha::get_physical_slha() const
{
   return physical_slha;
}

MSSMNoFVatMGUT_physical& MSSMNoFVatMGUT_slha::get_physical_slha()
{
   return physical_slha;
}

void MSSMNoFVatMGUT_slha::print(std::ostream& ostr) const
{
   MSSMNoFVatMGUT_mass_eigenstates::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
}

void MSSMNoFVatMGUT_slha::set_convert_masses_to_slha(bool flag)
{
   convert_masses_to_slha = flag;
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFVatMGUT_slha& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
