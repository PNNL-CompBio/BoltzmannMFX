//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include "bmx_pc.H"
#include <bmx.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fface (const amrex::Real& overlap, const amrex::Real& rad)
{
  amrex::Real del = overlap / rad;
  return 0.25*( 3.0 - del)*del*del;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fedge (const amrex::Real& ahat,
                   const amrex::Real& bhat)
{

  amrex::Real xhat_sq = 1.0 - ahat*ahat - bhat*bhat;
  if ( xhat_sq <= 0.0 ) {
    return 0.0;
  } else {

    amrex::Real xhat = sqrt(xhat_sq);

    return (1./(4.*M_PI)) * ( 2.0*ahat*bhat*xhat -
                            ( 3.*ahat - ahat*ahat*ahat) * std::atan(xhat/bhat) -
                            ( 3.*bhat - bhat*bhat*bhat) * std::atan(xhat/ahat) +
                            ( 2. * std::atan(xhat*ahat/bhat)) +
                            ( 2. * std::atan(xhat*bhat/ahat)));
  }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fcorner (const amrex::Real& ahat,
                     const amrex::Real& bhat,
                     const amrex::Real& chat,
                     const amrex::Real& fe12)
{

  if ( ahat*ahat + bhat*bhat + chat*chat >= 1.0 ) {
    return 0.0;

  } else {

    amrex::Real Ahat = sqrt(1.0 - ahat*ahat - chat*chat);
    amrex::Real Bhat = sqrt(1.0 - bhat*bhat - chat*chat);

    return 0.5*fe12 - (1./(8.*M_PI)) *
      (6.0*ahat*bhat*chat - 2.*ahat*Ahat*chat - 2.*bhat*Bhat*chat -
       (3.*bhat - bhat*bhat*bhat)*std::atan(chat/Bhat) -
       (3.*ahat - ahat*ahat*ahat)*std::atan(chat/Ahat) +
       (3.*chat - chat*chat*chat)*(std::atan(Ahat/ahat) - std::atan(bhat/Bhat)) +
       2.*(std::atan(chat*(ahat/Ahat)) + std::atan(chat*(bhat/Bhat))));

  }

}

struct OneToOneDeposition
{
  // These three variables are set when the OneToOneDeposition struct is passed
  // to the ParticleToMesh routine.
  int start_part_comp;
  int start_mesh_comp;
  int   num_comp;

  AMREX_GPU_DEVICE
  void operator()
                  (const BMXParticleContainer::ParticleType& p,
                   amrex::Array4<amrex::Real> const& rho,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) const noexcept
  {
      amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0];
      amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1];
      amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2];

      int i = static_cast<int>(amrex::Math::floor(lx));
      int j = static_cast<int>(amrex::Math::floor(ly));
      int k = static_cast<int>(amrex::Math::floor(lz));

      for (int icomp = 0; icomp < num_comp; icomp++)
          rho(i,j,k,start_mesh_comp+icomp) += p.rdata(start_part_comp+icomp);
  }
};

struct TrilinearDeposition
{
  // These three variables are set when the TrilinearDeposition struct is passed
  // to the ParticleToMesh routine.
  int start_part_comp;
  int start_mesh_comp;
  int   num_comp;

  AMREX_GPU_DEVICE
  void operator()
                  (const BMXParticleContainer::ParticleType& p,
                   amrex::Array4<amrex::Real> const& rho,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) const noexcept
  {
      amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
      amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
      amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

      int i = static_cast<int>(amrex::Math::floor(lx));
      int j = static_cast<int>(amrex::Math::floor(ly));
      int k = static_cast<int>(amrex::Math::floor(lz));

      amrex::Real xint = lx - i;
      amrex::Real yint = ly - j;
      amrex::Real zint = lz - k;

      amrex::Real wx[] = {1.-xint, xint};
      amrex::Real wy[] = {1.-yint, yint};
      amrex::Real wz[] = {1.-zint, zint};

      GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;
      amrex::Real total_weight = 0.0;

      for (int ii = 0; ii <= 1; ++ii)
        for (int jj = 0; jj <= 1; ++jj)
          for (int kk = 0; kk <= 1; ++kk){
            {
              weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
              total_weight += weights[ii][jj][kk];
            }
          }

      for (int ii = 0; ii <= 1; ++ii)
        for (int jj = 0; jj <= 1; ++jj)
          for (int kk = 0; kk <= 1; ++kk)
            weights[ii][jj][kk] /= total_weight;

      for (int kk = 0; kk <= 1; ++kk) {
          for (int jj = 0; jj <= 1; ++jj) {
              for (int ii = 0; ii <= 1; ++ii) {
                  for (int icomp = 0; icomp < num_comp; icomp++)
                  {
                      amrex::Gpu::Atomic::AddNoRet(&rho(i+ii-1, j+jj-1, k+kk-1, start_mesh_comp+icomp),
			                            weights[ii][jj][kk]*p.rdata(start_part_comp+icomp));
                  }
              }
          }
      }
  }
};

struct TrilinearDPVMSquareDeposition
{
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  void operator() (amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                   const amrex::RealVect& pos,
                   const amrex::Real radius,
                   int& i,
                   int& j,
                   int& k,
                   amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,2>,2>,2>& weights,
                   const amrex::Real deposition_scale_factor) const
  {
    amrex::Real x = (pos[0] - plo[0]) * dxi[0] + 0.5;
    amrex::Real y = (pos[1] - plo[1]) * dxi[1] + 0.5;
    amrex::Real z = (pos[2] - plo[2]) * dxi[2] + 0.5;

    i = static_cast<int>(amrex::Math::floor(x));
    j = static_cast<int>(amrex::Math::floor(y));
    k = static_cast<int>(amrex::Math::floor(z));

    amrex::GpuArray<amrex::Real,2> wx;
    amrex::GpuArray<amrex::Real,2> wy;
    amrex::GpuArray<amrex::Real,2> wz;

    amrex::Real rp  = amrex::max(dx[0], radius * deposition_scale_factor);

    amrex::Real i4rp3 = 1./(4.*rp*rp*rp);

    amrex::Real hh;

    if ( i == static_cast<int>(amrex::Math::floor(x-0.5)) ) {
      hh = (i*dx[0]) - pos[0] + rp;

      wx[0] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wx[1] = 1.0 - wx[0];

    } else {
      hh = pos[0] + rp - (i*dx[0]);

      wx[1] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wx[0] = 1.0 - wx[1];
    }

    if ( j == static_cast<int>(amrex::Math::floor(y-0.5)) ) {
      hh = (j*dx[1]) - pos[1] + rp;

      wy[0] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wy[1] = 1.0 - wy[0];

    } else {
      hh = pos[1] + rp - (j*dx[1]);

      wy[1] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wy[0] = 1.0 - wy[1];
    }


    if ( k == static_cast<int>(amrex::Math::floor(z-0.5)) ) {
      hh = (k*dx[2]) - pos[2] + rp;

      wz[0] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wz[1] = 1.0 - wz[0];

    } else {
      hh = pos[2] + rp - (k*dx[2]);

      wz[1] = ( hh < 0.0 ) ? 0.0 :
        (hh*hh*(3.*rp - hh)) * i4rp3;

      wz[0] = 1.0 - wz[1];
    }


    amrex::Real total_weight = 0.0;

    for (int ii = 0; ii <= 1; ++ii)
      for (int jj = 0; jj <= 1; ++jj)
        for (int kk = 0; kk <= 1; ++kk){
          {
            weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
            total_weight += weights[ii][jj][kk];
          }
        }

    for (int ii = 0; ii <= 1; ++ii)
      for (int jj = 0; jj <= 1; ++jj)
        for (int kk = 0; kk <= 1; ++kk)
          weights[ii][jj][kk] /= total_weight;
  }
};


struct TrueDPVMDeposition
{
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  void operator() (amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                   const amrex::RealVect& pos,
                   const amrex::Real radius,
                   int& i,
                   int& j,
                   int& k,
                   amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,2>,2>,2>& weights,
                   const amrex::Real deposition_scale_factor) const
  {
    amrex::Real x = (pos[0] - plo[0]) * dxi[0] + 0.5;
    amrex::Real y = (pos[1] - plo[1]) * dxi[1] + 0.5;
    amrex::Real z = (pos[2] - plo[2]) * dxi[2] + 0.5;

    i = static_cast<int>(amrex::Math::floor(x));
    j = static_cast<int>(amrex::Math::floor(y));
    k = static_cast<int>(amrex::Math::floor(z));

    amrex::GpuArray<amrex::Real,2> wx;
    amrex::GpuArray<amrex::Real,2> wy;
    amrex::GpuArray<amrex::Real,2> wz;

    amrex::Real rp  = amrex::max(dx[0], radius * deposition_scale_factor);

    // This gives us the i,j,k of the cell containing the centriod
    // of the particle.
    const int ip = static_cast<int>(amrex::Math::floor(x-0.5));
    const int jp = static_cast<int>(amrex::Math::floor(y-0.5));
    const int kp = static_cast<int>(amrex::Math::floor(z-0.5));

    amrex::RealVect overlap(0.);

    if ( i == ip ) {

      // particle sits on low side of i index
      overlap[0] = amrex::max(0.0, (i*dx[0]) - pos[0] + rp);

      wx[0] = (overlap[0] == 0.0) ? 0.0 : fface(overlap[0], rp);
      wx[1] = 1.0 - wx[0];

    } else {

      // particle sits on high side of i-1 index
      overlap[0] = amrex::max(0.0, pos[0] + rp - (i*dx[0]));

      wx[1] = (overlap[0] == 0.0) ? 0.0 : fface(overlap[0], rp);
      wx[0] = 1.0 - wx[1];

    }

    if ( j == jp ) {

      // particle sits on low side of j index
      overlap[1] = amrex::max(0.0, (j*dx[1]) - pos[1] + rp);

      wy[0] = (overlap[1] == 0.0) ? 0.0 : fface(overlap[1], rp);
      wy[1] = 1.0 - wy[0];

    } else {

      // particle sits on high side of j-1 index
      overlap[1] = amrex::max(0.0, pos[1] + rp - (j*dx[1]));

      wy[1] = (overlap[1] == 0.0) ? 0.0 : fface(overlap[1], rp);
      wy[0] = 1.0 - wy[1];

    }


    if ( k == kp ) {

      // particle sits on low side of k index
      overlap[2] = amrex::max(0.0, (k*dx[2]) - pos[2] + rp);

      wz[0] = (overlap[2] == 0.0) ? 0.0 : fface(overlap[2], rp);
      wz[1] = 1.0 - wz[0];

    } else {

      // particle sits on high side of k-1 index
      overlap[2] = amrex::max(0.0, pos[2] + rp - (k*dx[2]));

      wz[1] = (overlap[2] == 0.0) ? 0.0 : fface(overlap[2], rp);
      wz[0] = 1.0 - wz[1];

    }

    int cuts = 0;
    for ( int ii=0; ii<=2; ++ii) {
      if( overlap[ii] != 0. ) cuts += 1;
    }

    const int i0 = i-ip;
    const int i1 = amrex::Math::abs(1-i0);
    const int j0 = j-jp;
    const int j1 = amrex::Math::abs(1-j0);
    const int k0 = k-kp;
    const int k1 = amrex::Math::abs(1-k0);

    // if ( cuts == 0 || cuts == 1 ) we don't need to do anything as the
    // calculation of the initial 'weights' makes that all work out.

    if ( cuts == 2 ) {

      // This has two sub-cases that have to be addressed. The first is when
      // we only intersect the two faces and need to deposit into three cells;
      // the second is when we intersect the corner and need to deposit
      // into four cells.

      if ( overlap[0] == 0.0 ) {

        weights[i0][0][0] = 0.0;
        weights[i0][0][1] = 0.0;
        weights[i0][1][0] = 0.0;
        weights[i0][1][1] = 0.0;

        amrex::Real ahat = 1.0 - overlap[1] / rp;
        amrex::Real bhat = 1.0 - overlap[2] / rp;

        weights[i1][j0][k0] =  fedge(ahat, bhat);
        weights[i1][j1][k0] =  wy[j0] - weights[i1][j0][k0];
        weights[i1][j0][k1] =  wz[k0] - weights[i1][j0][k0];
        weights[i1][j1][k1] =  1.0 - weights[i1][j0][k0] - weights[i1][j0][k1] - weights[i1][j1][k0];

      } else if ( overlap[1] == 0.0 ) {

        weights[0][j0][0] = 0.0;
        weights[0][j0][1] = 0.0;
        weights[1][j0][0] = 0.0;
        weights[1][j0][1] = 0.0;

        amrex::Real ahat = 1.0 - overlap[0] / rp;
        amrex::Real bhat = 1.0 - overlap[2] / rp;

        weights[i0][j1][k0] =  fedge(ahat, bhat);
        weights[i1][j1][k0] =  wx[i0] - weights[i0][j1][k0];
        weights[i0][j1][k1] =  wz[k0] - weights[i0][j1][k0];
        weights[i1][j1][k1] =  1.0   - weights[i0][j1][k0] - weights[i1][j1][k0] - weights[i0][j1][k1];

      } else {

        weights[0][0][k0] = 0.0;
        weights[0][1][k0] = 0.0;
        weights[1][0][k0] = 0.0;
        weights[1][1][k0] = 0.0;

        amrex::Real ahat = 1.0 - overlap[0] / rp;
        amrex::Real bhat = 1.0 - overlap[1] / rp;

        weights[i0][j0][k1] =  fedge(ahat, bhat);
        weights[i1][j0][k1] =  wx[i0] - weights[i0][j0][k1];
        weights[i0][j1][k1] =  wy[j0] - weights[i0][j0][k1];
        weights[i1][j1][k1] =  1.0   - weights[i0][j0][k1] - weights[i1][j0][k1] - weights[i0][j1][k1];

      }

    } else if ( cuts == 3 ) {

      amrex::Real ahat = 1.0 - overlap[0] / rp;
      amrex::Real bhat = 1.0 - overlap[1] / rp;
      amrex::Real chat = 1.0 - overlap[2] / rp;

      amrex::Real edge12 = fedge(ahat, bhat);
      amrex::Real edge13 = fedge(ahat, chat);
      amrex::Real edge23 = fedge(bhat, chat);

      amrex::Real corner = fcorner(ahat, bhat, chat, edge12);

      // corner
      weights[i0][j0][k0] =  corner;

      // face overlaps
      weights[i0][j1][k1] =  wx[i0] - (edge12 + edge13) + corner;
      weights[i1][j0][k1] =  wy[j0] - (edge12 + edge23) + corner;
      weights[i1][j1][k0] =  wz[k0] - (edge13 + edge23) + corner;

      // edge overlaps
      weights[i0][j0][k1] =  edge12 - corner;
      weights[i0][j1][k0] =  edge13 - corner;
      weights[i1][j0][k0] =  edge23 - corner;

      // centroid
      weights[i1][j1][k1] =  1.0 - (wx[i0] + wy[j0] + wz[k0] -
                                    edge12 - edge13 - edge23 + corner);

    } else {

      weights[0][0][0] = wx[0] * wy[0] * wz[0];
      weights[0][0][1] = wx[0] * wy[0] * wz[1];
      weights[0][1][0] = wx[0] * wy[1] * wz[0];
      weights[0][1][1] = wx[0] * wy[1] * wz[1];
      weights[1][0][0] = wx[1] * wy[0] * wz[0];
      weights[1][0][1] = wx[1] * wy[0] * wz[1];
      weights[1][1][0] = wx[1] * wy[1] * wz[0];
      weights[1][1][1] = wx[1] * wy[1] * wz[1];

    }
  }
};


struct CentroidDeposition
{
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  void operator() (amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& /*dx*/,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                   const amrex::RealVect& pos,
                   const amrex::Real /*radius*/,
                   int& i,
                   int& j,
                   int& k,
                   amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,2>,2>,2>& weights,
                   const amrex::Real /*deposition_scale_factor*/) const
  {
    amrex::Real x = (pos[0] - plo[0]) * dxi[0] + 0.5;
    amrex::Real y = (pos[1] - plo[1]) * dxi[1] + 0.5;
    amrex::Real z = (pos[2] - plo[2]) * dxi[2] + 0.5;

    i = static_cast<int>(amrex::Math::floor(x));
    j = static_cast<int>(amrex::Math::floor(y));
    k = static_cast<int>(amrex::Math::floor(z));

    amrex::GpuArray<amrex::Real,2> wx;
    amrex::GpuArray<amrex::Real,2> wy;
    amrex::GpuArray<amrex::Real,2> wz;

    // This gives us the i,j,k of the cell containing the centriod
    // of the particle.
    const int ip = static_cast<int>(amrex::Math::floor(x-0.5));
    const int jp = static_cast<int>(amrex::Math::floor(y-0.5));
    const int kp = static_cast<int>(amrex::Math::floor(z-0.5));

    // particle sits on low side of i index otherwise in i-1
    wx[0] = ( i == ip ) ? 0.0 : 1.0;
    wx[1] = 1.0 - wx[0];

    // particle sits on low side of j index otherwise in j-1
    wy[0] = ( j == jp ) ? 0.0 : 1.0;
    wy[1] = 1.0 - wy[0];

    // particle sits on low side of k index otherwise in k-1
    wz[0] = ( k == kp ) ? 0.0 : 1.0;
    wz[1] = 1.0 - wz[0];

    // Note that we are not scaling for covered cells. In theory, a particle
    // centroid should never end of in a covered cell.

    for (int ii = 0; ii <= 1; ++ii)
      for (int jj = 0; jj <= 1; ++jj)
        for (int kk = 0; kk <= 1; ++kk)
          weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
  }
};
