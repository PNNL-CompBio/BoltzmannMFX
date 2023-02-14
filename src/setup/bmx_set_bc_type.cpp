//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>

#include <bmx_bc_parms.H>
#include <bmx_fluid_parms.H>


using namespace BC;

void
bmx::bmx_set_bc_type (int lev)
{
    const int und_  = bc_list.get_undefined();
    const int ig_   = bc_list.get_ig();

    const int l_chem_species = FLUID::nchem_species;

    // Set the defaults for BCRecs
    m_bcrec_chem_species.resize(l_chem_species);

    { // begin x direction

      const int dir = 0;

      Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
      Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();

      const int init_x = geom[lev].isPeriodic(0) ? und_ : ig_;

      Box domainx(geom[lev].Domain());
      domainx.grow(1,nghost);  // Add ghost cells to y
      domainx.grow(2,nghost);  // Add ghost cells to z

      { // x-lo side of the domain

        Box box_ilo = amrex::adjCellLo(domainx,0,1);

        IntVect ibx_lo(box_ilo.loVect());
        IntVect ibx_hi(box_ilo.hiVect());

        int xlo_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ilo, [bc_ilo_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ilo_type(i,j,k,0) = init_x;});

        set_bcrec_lo(lev, dir, xlo_type);

      } // end x-lo side of the domain


      { // x-hi side of the domain

        Box box_ihi = amrex::adjCellHi(domainx,0,1);

        IntVect ibx_lo(box_ihi.loVect());
        IntVect ibx_hi(box_ihi.hiVect());

        int xhi_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ihi, [bc_ihi_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ihi_type(i,j,k,0) = init_x;});

        set_bcrec_hi(lev, dir, xhi_type);

      } // end x-hi side of the domain
    } // end x-direction



    { // begin y-direction

      const int dir = 1;

      const int init_y = geom[lev].isPeriodic(1) ? und_ : ig_;

      Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
      Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();

      Box domainy(geom[lev].Domain());
      domainy.grow(0,nghost);  // Add ghost cells to x
      domainy.grow(2,nghost);  // Add ghost cells to z

      { // y-lo side of the domain

        Box box_jlo = amrex::adjCellLo(domainy,1,1);

        IntVect jbx_lo(box_jlo.loVect());
        IntVect jbx_hi(box_jlo.hiVect());

        int ylo_type = init_y;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_jlo, [bc_jlo_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jlo_type(i,j,k,0) = init_y;});

        set_bcrec_lo(lev, dir, ylo_type);

      }// end y-lo side of the domain


      { // y-hi side of the domain

        Box box_jhi = amrex::adjCellHi(domainy,1,1);

        IntVect jbx_lo(box_jhi.loVect());
        IntVect jbx_hi(box_jhi.hiVect());

        int yhi_type = init_y;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_jhi, [bc_jhi_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jhi_type(i,j,k,0) = init_y;});

        set_bcrec_hi(lev, dir, yhi_type);

      } // end y-hi side of the domain
    } // end y-direction




    { // begin z-direction

      const int dir = 2;

      const int init_z = geom[lev].isPeriodic(2) ? und_ : ig_;

      Array4<int> const& bc_klo_type = bc_klo[lev]->array();
      Array4<int> const& bc_khi_type = bc_khi[lev]->array();

      Box domainz(geom[lev].Domain());
      domainz.grow(0,nghost);  // Add ghost cells to x
      domainz.grow(1,nghost);  // Add ghost cells to y

      { // z-lo side of the domain

        Box box_klo = amrex::adjCellLo(domainz,2,1);

        IntVect kbx_lo(box_klo.loVect());
        IntVect kbx_hi(box_klo.hiVect());

        int zlo_type = init_z;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_klo, [bc_klo_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_klo_type(i,j,k,0) = init_z;});

        set_bcrec_lo(lev, dir, zlo_type);

      } // end z-lo side of the domain


      { // z-hi side of the domain

        Box box_khi = amrex::adjCellHi(domainz,2,1);

        IntVect kbx_lo(box_khi.loVect());
        IntVect kbx_hi(box_khi.hiVect());

        int zhi_type = init_z;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_khi, [bc_khi_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_khi_type(i,j,k,0) = init_z;});

        set_bcrec_hi(lev, dir, zhi_type);

      } // end z-hi side of the domain
    }// end z-direction


    if (l_chem_species > 0) {
      m_bcrec_chem_species_d.resize(l_chem_species);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_chem_species_d.data(), m_bcrec_chem_species.data(), sizeof(BCRec)*l_chem_species);
    }

    m_h_bc_X_k.resize(FLUID::nchem_species, Vector<Real>(bc.size()));

    if ( FLUID::solve and advect_fluid_chem_species) {
        for (int n = 0; n < FLUID::nchem_species; ++n) {
            Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_X_k[n].begin(), m_h_bc_X_k[n].end(),
                           m_bc_X_k[n].begin());
        }
    }
    Gpu::synchronize();
}


void bmx::set_bcrec_lo(const int lev, const int dir, const int /*l_type*/)
{
  // Scalar BC Recs
  if (geom[lev].isPeriodic(dir)) {

    for (auto& b : m_bcrec_chem_species) b.setLo(dir, BCType::int_dir);

  } else {

    for (auto& b : m_bcrec_chem_species) b.setLo(dir, BCType::foextrap);

  }
}

void bmx::set_bcrec_hi(const int lev, const int dir, const int /*l_type*/)
{
  // Scalar BC Recs
  if (geom[lev].isPeriodic(dir)) {

    for (auto& b : m_bcrec_chem_species) b.setHi(dir, BCType::int_dir);

  } else {

    for (auto& b : m_bcrec_chem_species) b.setHi(dir, BCType::foextrap);

  }
}
