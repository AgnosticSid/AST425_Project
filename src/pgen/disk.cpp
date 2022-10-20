//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk.cpp
//! \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//! spherical polar coordinates.  Initial conditions are in vertical hydrostatic eqm.

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <valarray>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
Real Ell_r(const Real rad, const Real phi, const Real z, int i);
Real Vel_cir(const Real rad, const Real phi, const Real z);
Real V_r(const Real rad, const Real phi, const Real z);
Real V_theta(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
Real dfloor;
Real Omega0;
Real a_0, a_f, e_0, d_e, n_fil, width;
} // namespace

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// void Gravity(MeshBlock *pmb, const Real time, const Real dt,
//              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//              AthenaArray<Real> &cons_scalar);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters for the shape of the ellipse
  a_0 = pin->GetOrAddReal("problem","a_0",0.0);
  a_f = pin->GetOrAddReal("problem","a_f",0.0);
  e_0 = pin->GetOrAddReal("problem","e_0",0.0);
  d_e = pin->GetOrAddReal("problem","d_e",0.0);
  n_fil = pin->GetOrAddReal("problem","n_fil",0.0);
  width = pin->GetOrAddReal("problem","width",0.0);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  }
  // EnrollUserExplicitSourceFunction(Gravity);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel, v_r, v_theta;
  Real x1, x2, x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        den = DenProfileCyl(rad,phi,z);
        vel = VelProfileCyl(rad,phi,z);
        v_r = V_r(rad,phi,z);
        v_theta = V_theta(rad,phi,z);
        if (porb->orbital_advection_defined)
          vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM1,k,j,i) = den*v_r;
          phydro->u(IM2,k,j,i) = den*v_theta;
          //std::cout << den << " , " << v_theta << " , " << vel << std::endl;
          // if ((rad >= (ell_r-0.1)) && (rad <= (ell_r+0.1))) {
          //   std::cout << phi << " , " << v_theta*rad << std::endl;
          // }
          // if ((rad >= Ell_r(rad,phi,z)-0.1) && (rad <= Ell_r(rad,phi,z)+0.1)) {
          //   std::cout << rad << " , " << 1/((std::pow(1/0.6,3/2)*std::pow(rad*215.032,-5/2))/107) << std::endl;
          // }
          // std::cout << rad << " , " << phi << " , " << z << " , " << den*rad*v_theta << std::endl;
          phydro->u(IM3,k,j,i) = 0.0;

        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = den*vel;
        }

        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates

Real Ell_r(const Real rad, const Real phi, const Real z, int i) {
  Real ell_r;
  Real a_ini = a_0;
  Real a_fin = a_f;
  Real nfil = n_fil;
  Real da = (a_f-a_0)/nfil;
  Real a_i = a_0 + i*da;
  Real e_ini = e_0;
  Real de = d_e;
  Real e_i = e_0; //+ i*de;
  ell_r = (a_i*(1-std::pow(e_i,2)))/(1+e_i*std::cos(phi));
  return ell_r;
}

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real dentem;
  Real den_i;
  Real p_over_r = p0_over_r0;
  // if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);

  // std::vector< double > arr;
  // for (int i=0; i<=10; i++) {
  //   arr.push_back(Ell_r(rad,phi,z,i));
  // }
  // for (double i: arr)
  //   std::cout << i << std::endl;

  dentem = 1.0e-1;
  // Real nfil = n_fil;
  // int i = 0;
  // while (i <= nfil) {
  //   if ((rad >= Ell_r(rad,phi,z,i)-width) && (rad <= Ell_r(rad,phi,z,i)+width))
  //     den_i = 6.4 - ((6.4-0.3)/nfil)*i;
  //     dentem = den_i*std::exp(-std::pow((rad-Ell_r(rad,phi,z,i)),2)); // /(std::sqrt(2)*0.04)
  //   i++;
  // }

  if ((rad >= Ell_r(rad,phi,z,0)-width) && (rad <= Ell_r(rad,phi,z,0)+width)) {
    dentem = 6.4*std::exp(-std::pow((rad-Ell_r(rad,phi,z,0)),2));
  }

  // else if ((rad >= Ell_r(rad,phi,z,0)+width)) {
  //   dentem = 0.1;
  // }
  // else {
  //   dentem = 3.1;
  // }
  
  // else {
  //   den = 1.0;
  // }

  // int i = 0;
  // do {
  //   if ((rad >= Ell_r(rad,phi,z,i)-0.1) && (rad <= Ell_r(rad,phi,z,i)+0.1)) {
  //     den = 6.4;
  //   }
  //   else {
  //     den = 1.0;
  //   }
  //   i++;
  // } while (i <= 10);

  // while (i <= 10) {
  // // for (int i=0; i<=10; i++) {
  //   if ((rad >= Ell_r(rad,phi,z,i)-0.01) && (rad <= Ell_r(rad,phi,z,i)+0.01)) {
  //     den = 6.4;
  //   }
  //   else {
  //     den = 1.0;
  //   }
  //   i++;
  // // } 
  // }
  // if ((rad >= Ell_r(rad,phi,z)-0.1) && (rad <= Ell_r(rad,phi,z)+0.1)) {
  //   den = 6.4;
  // }
  // if ((phi >= (0)-0.3) && (phi <= (0)+0.3)) {
  //   if ((rad >= Ell_r(rad,phi,z)-0.1) && (rad <= Ell_r(rad,phi,z)+0.1)) {
  //     den = 6.4;
  //   }
  //   else {
  //     den = 1.0;
  //   }
  // }
    // std::cout << Ell_r(rad,phi,z) << std::endl;
  // else if ((rad >= Ell_r_1(rad,phi,z)-0.1) && (rad <= Ell_r_1(rad,phi,z)+0.1)) {
  //   den = 6.4;
  //   // std::cout << Ell_r(rad,phi,z) << std::endl;
  // }
  // else {
  //   den = 1.0;
  // }
  // }
  // else {
  //   den = 1.0;
  // }
  // else {
  //   den = dfloor;
  // }
  //}
  // if ((phi >= -0.3) && (phi <= +0.3)) {
  //   if ((rad >= Ell_r(rad,phi,z)-0.1) && (rad <= Ell_r(rad,phi,z)+0.1)) {
  //     den = 6.4;
  //   }
  //   else {
  //     den = dfloor;
  // }
  // }
  // Real denmid = rho0*std::pow(rad/r0,dslope);
  // den = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));

  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! computes rotational velocity in cylindrical coordinates

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}

// Real Vel_cir(const Real rad, const Real phi, const Real z) {
//   Real vel_cir;
//   Real a = a_0;
//   Real P = std::sqrt(std::pow(a,3)/gm0);
//   vel_cir = a/P;
//   return vel_cir;
// }

Real V_r(const Real rad, const Real phi, const Real z) {
  Real v_r;
  v_r = 0.0;
  Real nfil = n_fil;
  int i = 0;
  while (i <= nfil) {
    if ((rad >= Ell_r(rad,phi,z,i)-width) && (rad <= Ell_r(rad,phi,z,i)+width)) {
      Real a_ini = a_0;
      Real a_fin = a_f;
      Real nfil = n_fil;
      Real da = (a_f-a_0)/nfil;
      Real a_i = a_0 + i*da;
      Real e_ini = e_0;
      Real de = d_e;
      Real e_i = e_0; //+ i*de; 
      v_r = VelProfileCyl(a_i,phi,z)*(e_i*std::sin(phi))/(std::sqrt((1-std::pow(e_i,2))));
    }
    i++;
  }

  // if ((rad >= Ell_r(rad,phi,z,0)-width) && (rad <= Ell_r(rad,phi,z,0)+width)) {
  //   Real a = a_0;
  //   Real P = std::sqrt(std::pow(a,3)/gm0);
  //   Real e = e_0;
  //   v_r = VelProfileCyl(a,phi,z)*(e*std::sin(phi))/(std::sqrt((1-std::pow(e,2))));
  // }

  // else if ((rad >= Ell_r_1(rad,phi,z)-0.1) && (rad <= Ell_r_1(rad,phi,z)+0.1)) {
  //   Real a = a_0_1;
  //   Real P = std::sqrt(std::pow(a,3)/gm0);
  //   Real e = e_0;
  //   v_r = VelProfileCyl(a,phi,z)*(e*std::sin(phi))/(std::sqrt((1-std::pow(e,2))));
  // } 

  // else {
  //   v_r = 0.0;
  // }
  return v_r;
}

Real V_theta(const Real rad, const Real phi, const Real z) {
  Real v_theta;
  v_theta = VelProfileCyl(rad,phi,z);
  Real nfil = n_fil;
  int i = 0;
  while (i <= nfil) {
    if ((rad >= Ell_r(rad,phi,z,i)-width) && (rad <= Ell_r(rad,phi,z,i)+width)) {
      Real a_ini = a_0;
      Real a_fin = a_f;
      Real nfil = n_fil;
      Real da = (a_f-a_0)/nfil;
      Real a_i = a_0 + i*da;
      Real e_ini = e_0;
      Real de = d_e;
      Real e_i = e_0; //+ i*de; 
      v_theta = VelProfileCyl(a_i,phi,z)*(1+e_i*std::cos(phi))/(std::sqrt((1-std::pow(e_i,2))));
    }
    i++;
  }

  // if ((rad >= Ell_r(rad,phi,z,0)-width) && (rad <= Ell_r(rad,phi,z,0)+width)) {
  //   Real a = a_0;
  //   Real P = std::sqrt(std::pow(a,3)/gm0);
  //   Real e = e_0;
  //   v_theta = VelProfileCyl(a,phi,z)*(1+e*std::cos(phi))/(std::sqrt((1-std::pow(e,2))));
  // }

  // else if ((rad >= Ell_r_1(rad,phi,z)-0.1) && (rad <= Ell_r_1(rad,phi,z)+0.1)) {
  //   Real a = a_0_1;
  //   Real P = std::sqrt(std::pow(a,3)/gm0);
  //   Real e = e_0;
  //   v_theta = VelProfileCyl(a,phi,z)*(1+e*std::cos(phi))/(std::sqrt((1-std::pow(e,2))));
  // }

  // else {
  //   v_theta = VelProfileCyl(rad,phi,z);
  // } //if ((phi >= M_PI-0.3) && (phi <= M_PI+0.3)){
    //  std::cout << "v_theta is equal to " << v_theta << "\n";
  //}
  return v_theta;
}

} // namespace

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = v_theta;
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = v_theta;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = v_theta;
          prim(IM3,k,jl-j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = 0.0;
          prim(IM3,k,jl-j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = v_theta;
          prim(IM3,k,ju+j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = 0.0;
          prim(IM3,k,ju+j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = v_theta;
          prim(IM3,kl-k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = 0.0;
          prim(IM3,kl-k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  Real v_theta;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          v_theta = V_theta(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = v_theta;
          prim(IM3,ku+k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = 0.0;
          prim(IM3,ku+k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  }
}

// Real AngularMom_tot(MeshBlock *pmb, int iout)
// {
//   Real angmom_tot = 0;
//   int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
//   AthenaArray<Real> angmomtot;
//   for (int k=ks; k<=ke; ++k) {
//     x3 = pcoord->x3v(k);
//     for (int j=js; j<=je; ++j) {
//       x2 = pcoord->x2v(j);
//       for (int i=is; i<=ie; ++i) {
//         x1 = pcoord->x1v(i);
//         GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
//         // compute initial conditions in cylindrical coordinates
//         den = DenProfileCyl(rad,phi,z);
//         vel = VelProfileCyl(rad,phi,z);
//         v_r = V_r(rad,phi,z);
//         v_theta = V_theta(rad,phi,z);
//         if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
//           phydro->u(IM1,k,j,i) = den*v_theta*x1f(k,j,i);
//         }  
//       }  
//     }
//   }
// }