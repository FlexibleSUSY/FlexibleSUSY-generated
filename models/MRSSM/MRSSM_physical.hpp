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

// File generated at Tue 7 Jul 2015 12:34:45

#ifndef MRSSM_PHYSICAL_H
#define MRSSM_PHYSICAL_H

#include "linalg2.hpp"
#include <Eigen/Core>

#include <iosfwd>
#include <string>

namespace flexiblesusy {

struct MRSSM_physical {
   MRSSM_physical();
   void clear();
   void convert_to_hk();   ///< converts pole masses to HK convention
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   void print(std::ostream&) const;

   double MGlu;
   Eigen::Array<double,3,1> MFv;
   double MSRdp;
   double MSRum;
   double MsigmaO;
   double MphiO;
   double MVZ;
   Eigen::Array<double,6,1> MSd;
   Eigen::Array<double,3,1> MSv;
   Eigen::Array<double,6,1> MSu;
   Eigen::Array<double,6,1> MSe;
   Eigen::Array<double,4,1> Mhh;
   Eigen::Array<double,4,1> MAh;
   Eigen::Array<double,2,1> MRh;
   Eigen::Array<double,4,1> MHpm;
   Eigen::Array<double,4,1> MChi;
   Eigen::Array<double,2,1> MCha1;
   Eigen::Array<double,2,1> MCha2;
   Eigen::Array<double,3,1> MFe;
   Eigen::Array<double,3,1> MFd;
   Eigen::Array<double,3,1> MFu;
   double MVG;
   double MVP;
   double MVWm;

   Eigen::Matrix<double,6,6> ZD;
   Eigen::Matrix<double,3,3> ZV;
   Eigen::Matrix<double,6,6> ZU;
   Eigen::Matrix<double,6,6> ZE;
   Eigen::Matrix<double,4,4> ZH;
   Eigen::Matrix<double,4,4> ZA;
   Eigen::Matrix<double,2,2> ZHR;
   Eigen::Matrix<double,4,4> ZP;
   Eigen::Matrix<std::complex<double>,4,4> ZN1;
   Eigen::Matrix<std::complex<double>,4,4> ZN2;
   Eigen::Matrix<std::complex<double>,2,2> UM1;
   Eigen::Matrix<std::complex<double>,2,2> UP1;
   Eigen::Matrix<std::complex<double>,2,2> UM2;
   Eigen::Matrix<std::complex<double>,2,2> UP2;
   Eigen::Matrix<std::complex<double>,3,3> ZEL;
   Eigen::Matrix<std::complex<double>,3,3> ZER;
   Eigen::Matrix<std::complex<double>,3,3> ZDL;
   Eigen::Matrix<std::complex<double>,3,3> ZDR;
   Eigen::Matrix<std::complex<double>,3,3> ZUL;
   Eigen::Matrix<std::complex<double>,3,3> ZUR;

};

std::ostream& operator<<(std::ostream&, const MRSSM_physical&);

} // namespace flexiblesusy

#endif
