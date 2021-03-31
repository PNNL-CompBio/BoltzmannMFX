#include <bmx.H>
#include <bmx_des_K.H>
#include <bmx_interp_K.H>
#include <bmx_filcc.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <bmx_mf_helpers.H>
#include <bmx_dem_parms.H>

#ifdef NEW_CHEM
#include <bmx_chem.H>
#endif

/**
 * @brief this function transfers data from particles to the continuum chemical
 * species fields defined on the AMR grid
 */
void
bmx::bmx_calc_txfr_fluid (Real time, Real dt)
{
  const Real strttime = ParallelDescriptor::second();

  for (int lev = 0; lev < nlev; lev++)
    m_leveldata[lev]->X_rhs->setVal(0);

  if (nlev > 2)
    amrex::Abort("For right now"
        " BMXParticleContainer::TrilinearDepositionFluidDragForce can only"
        " handle up to 2 levels");

  Vector< MultiFab* > txfr_ptr(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++) {

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 and OnSameGrids) {

      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      txfr_ptr[lev] = m_leveldata[lev]->X_rhs;

    } else if (lev == 0 and (not OnSameGrids)) {

      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                              pc->ParticleDistributionMap(lev),
                                              m_leveldata[lev]->X_rhs->nComp(),
                                              m_leveldata[lev]->X_rhs->nGrow());

    } else {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
      txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                   m_leveldata[lev]->X_rhs->nComp(), 1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (txfr_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    txfr_ptr[lev]->setVal(0.0, 0, m_leveldata[lev]->X_rhs->nComp(), txfr_ptr[lev]->nGrow());
  }

  const Geometry& gm = Geom(0);

  // Deposit the chem_species_rhs to the grid
  for (int lev = 0; lev < nlev; lev++) {
    pc->InterphaseTxfrDeposition(lev, *txfr_ptr[lev], dt); 
  }

  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    bmx_deposition_bcs(lev, *txfr_ptr[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    txfr_ptr[lev]->SumBoundary(gm.periodicity());
    txfr_ptr[lev]->setBndry(0.0);
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from txfr_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (txfr_ptr[0] != m_leveldata[0]->X_rhs) {
    m_leveldata[0]->X_rhs->copy(*txfr_ptr[0], 0, 0, m_leveldata[0]->X_rhs->nComp());
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (txfr_ptr[lev] != m_leveldata[lev]->X_rhs)
      delete txfr_ptr[lev];
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "BMXParticleContainer::TrilinearDepositionFluidDragForce time: " << stoptime << '\n';
  }
}

//
// Interpolate fluid chem_species onto particle locations
//
void
bmx::bmx_calc_txfr_particle (Real time, Real dt)
{
  using BMXParIter = BMXParticleContainer::BMXParIter;

#ifdef NEW_CHEM
  BMXChemistry *bmxchem = BMXChemistry::instance();
#endif
  //
  BL_PROFILE("bmx::bmx_calc_txfr_particle()");

  bmx_set_chem_species_bcs(time, get_X_k(), get_D_k());

  for (int lev = 0; lev < nlev; lev++)
  {
    Box domain(geom[lev].Domain());

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    // Pointer to Multifab for interpolation
    MultiFab* interp_ptr;

#ifdef NEW_CHEM
#if 0
    int nint, nreal, tot_int, tot_real;
    bmxchem->getVarArraySizes(&nint, &nreal, &tot_int, &tot_real);

    const int interp_ng    = 1;    // Only one layer needed for interpolation
    const int interp_ncomp = static_cast<const int>(nreal);

    if (m_leveldata[0]->X_k->nComp() != nreal)
      amrex::Abort("We are not interpolating the right number of components in calc_txfr_particle");
#endif
    const int interp_ng    = 1;    // Only one layer needed for interpolation
    const int interp_ncomp = 3;

    if (m_leveldata[0]->X_k->nComp() != 3)
      amrex::Abort("We are not interpolating the right number of components in calc_txfr_particle");
#else
    const int interp_ng    = 1;    // Only one layer needed for interpolation
    const int interp_ncomp = 2;

    if (m_leveldata[0]->X_k->nComp() != 2)
      amrex::Abort("We are not interpolating the right number of components in calc_txfr_particle");
#endif

    if (OnSameGrids)
    {
      // Store X_k for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_ncomp, interp_ng, MFInfo());

      // Copy 
      interp_ptr->copy(*m_leveldata[lev]->X_k, 0, 0,
                        m_leveldata[lev]->X_k->nComp(),
                        interp_ng, interp_ng);
      interp_ptr->FillBoundary(geom[lev].periodicity());

    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Store X_k for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_ncomp, interp_ng, MFInfo());

      // Copy 
      interp_ptr->copy(*m_leveldata[lev]->X_k, 0, 0,
                        m_leveldata[lev]->X_k->nComp(),
                        interp_ng, interp_ng);

      interp_ptr->FillBoundary(geom[lev].periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      Real grid_vol = dx[0]*dx[1]*dx[2];

      for (BMXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        const auto& interp_array = interp_ptr->array(pti);

        amrex::ParallelFor(np,
            [pstruct,interp_array,plo,dxi,bmxchem,grid_vol,dt]
            AMREX_GPU_DEVICE (int pid) noexcept
              {
              // Local array storing interpolated values
              GpuArray<Real, interp_ncomp> interp_loc;

              BMXParticleContainer::ParticleType& p = pstruct[pid];

              trilinear_interp(p.pos(), &interp_loc[0],
                               interp_array, plo, dxi, interp_ncomp);

#ifdef NEW_CHEM
              Real cell_vol = p.rdata(realIdx::vol);
              Real cell_area = p.rdata(realIdx::area);
              Real *p_vals = &p.rdata(realIdx::first_data);
              bmxchem->xferMeshToParticle(grid_vol, cell_vol, cell_area,
                        &interp_loc[0], p_vals, dt);
              bmxchem->updateChemistry(p_vals, dt);
#else
              // Interpolate values from mesh to particles
              p.rdata(realData::fluid_A) = interp_loc[0];
              p.rdata(realData::fluid_B) = interp_loc[1];

              amrex::Print() << "VALUE OF X ON PARTICLE " << 
                  p.rdata(realData::fluid_A) << " " << 
                  p.rdata(realData::fluid_B) << std::endl;;

              // The particle will consume (dt * 50%) of what the mesh value is
              p.rdata(realData::consume_A) = -0.5*interp_loc[0];
              p.rdata(realData::consume_B) = -0.5*interp_loc[1];
#endif
            });
      } // pti
    } // omp region

    delete interp_ptr;

  } // lev
}
