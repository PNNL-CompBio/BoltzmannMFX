#include <AMReX.H>
#include <AMReX_Particles.H>
#include <bmx_pc.H>

#include <bmx_deposition_K.H>
#include <bmx_dem_parms.H>
#include <bmx_chem_species_parms.H>
#include <bmx_fluid_parms.H>
#include <bmx_algorithm.H>

#ifdef NEW_CHEM
#include <bmx_chem.H>
#endif

using namespace amrex;

void BMXParticleContainer::
SolidsVolumeDeposition (int lev,
                  amrex::MultiFab & mf_to_be_filled)
{

  if (bmx::m_deposition_scheme == DepositionScheme::trilinear) {

    SolidsVolumeDeposition(TrilinearDeposition(),
                     lev, mf_to_be_filled);

  } else if (bmx::m_deposition_scheme == DepositionScheme::square_dpvm) {

    SolidsVolumeDeposition(TrilinearDPVMSquareDeposition(),
                     lev, mf_to_be_filled);

  } else if (bmx::m_deposition_scheme == DepositionScheme::true_dpvm) {

    SolidsVolumeDeposition(TrueDPVMDeposition(),
                     lev, mf_to_be_filled);

  } else if (bmx::m_deposition_scheme == DepositionScheme::centroid) {

    SolidsVolumeDeposition(CentroidDeposition(),
                     lev, mf_to_be_filled);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }


}

template <typename F>
void BMXParticleContainer::
SolidsVolumeDeposition (F WeightFunc, int lev,
                        amrex::MultiFab & mf_to_be_filled)
{
  BL_PROFILE("BMXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_to_be_filled;

    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = mf_to_be_filled[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      {
        auto volarr = fab_to_be_filled.array();

#ifdef _OPENMP
        const int ncomp = mf_to_be_filled.nComp();
        Box tile_box = bx;

        if(Gpu::notInLaunchRegion())
        {
          tile_box.grow(mf_to_be_filled.nGrow());
          local_fab_to_be_filled.resize(tile_box, ncomp);
          local_fab_to_be_filled.setVal<RunOn::Host>(0.0);
          volarr = local_fab_to_be_filled.array();
        }
#endif

        const amrex::Real deposition_scale_factor =
          bmx::m_deposition_scale_factor;

        amrex::ParallelFor(nrp,
          [pstruct,p_realarray,plo,dx,dxi,deposition_scale_factor,volarr,reg_cell_vol,WeightFunc]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            WeightFunc(plo, dx, dxi, p.pos(), p_realarray[realData::radius][ip], i, j, k, weights,
                deposition_scale_factor);

            amrex::Real pvol = p_realarray[realData::volume][ip] / reg_cell_vol;

            for (int kk = -1; kk <= 0; ++kk) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int ii = -1; ii <= 0; ++ii) {
                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];
                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                }
              }
            }
          });

#ifdef _OPENMP
        if(Gpu::notInLaunchRegion())
        {
          fab_to_be_filled.atomicAdd<RunOn::Host>(local_fab_to_be_filled,
              tile_box, tile_box, 0, 0, ncomp);
        }
#endif

      }
    }
  }
}

void BMXParticleContainer::
InterphaseTxfrDeposition (int lev,
                          amrex::MultiFab & txfr_mf, Real dt)
{
  if (bmx::m_deposition_scheme == DepositionScheme::trilinear) {

    InterphaseTxfrDeposition(TrilinearDeposition(), lev, txfr_mf, dt);

  } else if (bmx::m_deposition_scheme == DepositionScheme::square_dpvm) {

    InterphaseTxfrDeposition(TrilinearDPVMSquareDeposition(), lev, txfr_mf, dt);

  } else if (bmx::m_deposition_scheme == DepositionScheme::true_dpvm) {

    InterphaseTxfrDeposition(TrueDPVMDeposition(), lev, txfr_mf, dt);

  } else if (bmx::m_deposition_scheme == DepositionScheme::centroid) {

    InterphaseTxfrDeposition(CentroidDeposition(), lev, txfr_mf, dt);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }

}

template <typename F>
void BMXParticleContainer::
InterphaseTxfrDeposition (F WeightFunc, int lev,
                          amrex::MultiFab & txfr_mf, Real dt)
{
  BL_PROFILE("BMXParticleContainer::InterphaseTxfrDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      grid_vol = dx[0]*dx[1]*dx[2];

#ifdef NEW_CHEM
  BMXChemistry *bmxchem = BMXChemistry::instance();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_txfr;

    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      const long nrp = pti.numParticles();

      FArrayBox& txfr_fab = txfr_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      auto        txfr_arr = txfr_fab.array();

      const amrex::Real deposition_scale_factor = bmx::m_deposition_scale_factor;

#ifdef _OPENMP
        const int ncomp = txfr_mf.nComp();
        Box tile_box = box;

        if (Gpu::notInLaunchRegion())
        {
          tile_box.grow(txfr_mf.nGrow());
          local_txfr.resize(tile_box, ncomp);
          local_txfr.setVal<RunOn::Host>(0.0);
          txfr_arr = local_txfr.array();
        }
#endif

        amrex::Print() << "DEPOSITION OF " << nrp << " particles ... " << std::endl;

        amrex::ParallelFor(nrp,
#ifdef NEW_CHEM
          [pstruct,plo,dx,dxi,deposition_scale_factor,WeightFunc,txfr_arr,bmxchem,grid_vol,dt]
#else
          [pstruct,plo,dx,dxi,deposition_scale_factor,WeightFunc,txfr_arr]
#endif
           AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

#ifdef NEW_CHEM
            WeightFunc(plo, dx, dxi, p.pos(), p.rdata(realIdx::a_size), i, j, k, weights,
                       deposition_scale_factor);

            int nvals = p.idata(intIdx::num_reals);
            amrex::Real *p_vals = &p.rdata(realIdx::first_data);
            amrex::Real *chem_incr = p_vals + p.idata(intIdx::first_real_inc);
            amrex::Real cell_vol = p.rdata(realIdx::vol);
            amrex::Real cell_area = p.rdata(realIdx::area);
            bmxchem->xferParticleToMesh(grid_vol, grid_vol, cell_area, chem_incr, p_vals, dt);
            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  for (int nn = 0; nn < nvals; ++nn) {

                    amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];

                    amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,0), weight_vol*chem_incr[nn]);
                  }

                }
              }
            }
#else
            WeightFunc(plo, dx, dxi, p.pos(), p.rdata(realData::radius), i, j, k, weights,
                       deposition_scale_factor);

            amrex::Print() << "IJK FOR DEPOSITION " << IntVect(i,j,k) << std::endl;

            amrex::Real p_A = p.rdata(realData::consume_A);
            amrex::Real p_B = p.rdata(realData::consume_B);

            amrex::Print() << "P_A / P_B DEPOSITION " << p_A << " " << p_B << std::endl;

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {

                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];

                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,0), weight_vol*p_A);
                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,1), weight_vol*p_B);

                }
              }
            }
#endif
          });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion())
        {
          txfr_fab.atomicAdd<RunOn::Host>(local_txfr, tile_box, tile_box, 0, 0, ncomp);
        }
#endif

    }
  }
}
