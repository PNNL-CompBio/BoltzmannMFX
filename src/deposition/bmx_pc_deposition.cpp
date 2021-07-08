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
  const Geometry& gm  = Geom(lev);
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

      //const auto& particles = pti.GetArrayOfStructs();
      auto& particles = pti.GetArrayOfStructs();
      //const ParticleType* pstruct = particles().dataPtr();
      ParticleType* pstruct = particles().dataPtr();

      const long nrp = pti.numParticles();

      FArrayBox& txfr_fab = txfr_mf[pti];

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

//        amrex::Print() << "DEPOSITION OF " << nrp << " particles ... " << std::endl;

        amrex::ParallelFor(nrp,
#ifdef NEW_CHEM
          [pstruct,plo,dx,dxi,deposition_scale_factor,WeightFunc,txfr_arr,bmxchem,grid_vol,dt]
#else
          [pstruct,plo,dx,dxi,deposition_scale_factor,WeightFunc,txfr_arr]
#endif
           AMREX_GPU_DEVICE (int ip) noexcept
          {
            //const ParticleType& p = pstruct[ip];
            ParticleType& p = pstruct[ip];

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

            amrex::Real *cell_par = &p.rdata(0);
            bmxchem->xferParticleToMesh(grid_vol, cell_par, chem_incr, p_vals, dt);
            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  for (int nn = 0; nn < nvals; ++nn) {

                    amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];

                    amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,nn), weight_vol*chem_incr[nn]);
                  }

                }
              }
            }
            bmxchem->printCellConcentrations(p_vals, cell_par);
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
