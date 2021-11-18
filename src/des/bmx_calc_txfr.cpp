#include <bmx.H>
#include <bmx_deposition_K.H>
#include <bmx_interp_K.H>
#include <bmx_fluid_parms.H>
#include <AMReX_AmrParticles.H>

#ifdef NEW_CHEM
#include <bmx_chem.H>
#include <bmx_chem_K.H>
#endif

/**
 * @brief this function transfers data from particles to the continuum chemical
 * species fields defined on the AMR grid
 */
void
bmx::bmx_calc_txfr_fluid (Real /*time*/, Real /*dt*/)
{
  BMXChemistry *bmxchem = BMXChemistry::instance();
  int inc_start = bmxchem->getIntData(intIdx::first_real_inc);
  int start_part_comp = realIdx::first_data + inc_start;
  int start_mesh_comp = 0;
  int        num_comp = FLUID::nchem_species;

  // Initialize to zero because the deposition routine will only change values
  // where there are particles (note this is the default)
  bool zero_out_input = true;

  // Here we don't divide the quantity on the mesh after deposition
  // (note the default is true)
  bool vol_weight     = false;

  if (bmx::m_cnc_deposition_scheme == DepositionScheme::one_to_one) {

     ParticleToMesh(*pc,get_X_rhs(),0,finest_level,
                    OneToOneDeposition{start_part_comp,start_mesh_comp,num_comp},
                    zero_out_input, vol_weight);

  } else if (bmx::m_cnc_deposition_scheme == DepositionScheme::trilinear) {

     ParticleToMesh(*pc,get_X_rhs(),0,finest_level,
                    TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp},
                    zero_out_input, vol_weight);

#if 0
  } else if (bmx::m_cnc_deposition_scheme == DepositionScheme::square_dpvm) {

    ParticleToMesh(*pc,get_X_rhsvf(),0,finest_level,DPVMSquareDeposition());

  } else if (bmx::m_cnc_deposition_scheme == DepositionScheme::true_dpvm) {

    ParticleToMesh(*pc,get_X_rhsvf(),0,finest_level,TrueDPVMDeposition());

  } else if (bmx::m_cnc_deposition_scheme == DepositionScheme::centroid) {

    ParticleToMesh(*pc,get_X_rhsvf(),0,finest_level,CentroidDeposition());
#endif

  } else {
    amrex::Abort("Don't know this deposition_scheme for concentrations!");
  }
}

/**
 * @brief this function interpolates fluid chem_species onto particle locations
 */
void
bmx::bmx_calc_txfr_particle (Real time, Real dt)
{
  using BMXParIter = BMXParticleContainer::BMXParIter;

#ifdef NEW_CHEM
  BMXChemistry *bmxchem = BMXChemistry::instance();
  std::vector<Real> chempar_vec;
  bmxchem->getChemParams(chempar_vec);
  Real *chempar = &chempar_vec[0];
#endif
  //
  BL_PROFILE("bmx::bmx_calc_txfr_particle()");

  bmx_set_chem_species_bcs(time, get_X_k(), get_D_k());

  long nparticles = 0;
  for (int lev = 0; lev <= finest_level; lev++)
  {
    Box domain(geom[lev].Domain());

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    nparticles += pc->NumberOfParticlesAtLevel(lev);

    // Pointer to Multifab for interpolation
    MultiFab* interp_ptr;

    // Pointer to Multifab for volume fraction
    MultiFab* interp_vptr;

    // Create Multifab with number of particles in grid cell
    MultiFab temp_npart(grids[lev], dmap[lev], 1, 0);
    temp_npart.setVal(0);
    pc->Increment(temp_npart, lev);
    MultiFab* interp_nptr;

#ifdef NEW_CHEM
    const int interp_ng    = 1;    // Only one layer needed for interpolation
    //const int interp_ncomp = bmxchem->getIntData(intIdx::num_reals);
    const int interp_ncomp = 3;

    if (m_leveldata[lev]->X_k->nComp() != interp_ncomp)
      amrex::Abort("We are not interpolating the right number of components in calc_txfr_particle");
#endif

    if (OnSameGrids)
    {
      // Store X_k for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_ncomp, interp_ng, MFInfo());

      // Copy 
      MultiFab::Copy(*interp_ptr,*m_leveldata[lev]->X_k, 0, 0,
                      m_leveldata[lev]->X_k->nComp(), interp_ng);
      interp_ptr->FillBoundary(geom[lev].periodicity());

      // Store vf_n for interpolation
      interp_vptr = new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo());

      // Copy 
      MultiFab::Copy(*interp_vptr,*m_leveldata[lev]->vf_n, 0, 0,
                      m_leveldata[lev]->vf_n->nComp(), 1);
      interp_vptr->FillBoundary(geom[lev].periodicity());

      // Store n_part for interpolation
      interp_nptr = new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo());

      // Copy 
      MultiFab::Copy(*interp_nptr,temp_npart, 0, 0, temp_npart.nComp(), 0);
      interp_nptr->FillBoundary(geom[lev].periodicity());

    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Store X_k for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_ncomp, interp_ng, MFInfo());

      // Copy 
      interp_ptr->ParallelCopy(*m_leveldata[lev]->X_k, 0, 0,
                                m_leveldata[lev]->X_k->nComp(),
                                interp_ng, interp_ng);

      interp_ptr->FillBoundary(geom[lev].periodicity());

      // Store vf_n for interpolation
      interp_vptr = new MultiFab(pba, pdm, 1, 1, MFInfo());

      // Copy 
      interp_vptr->ParallelCopy(*m_leveldata[lev]->vf_n, 0, 0,
                                m_leveldata[lev]->vf_n->nComp(),
                                1, 1);

      interp_vptr->FillBoundary(geom[lev].periodicity());

      // Store n_part for interpolation
      interp_nptr = new MultiFab(pba, pdm, 1, 1, MFInfo());

      // Copy 
      interp_nptr->ParallelCopy(temp_npart, 0, 0,
                                temp_npart.nComp(),
                                1, 1);

      interp_nptr->FillBoundary(geom[lev].periodicity());
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

      if (m_verbose != 0) {
        for (BMXParIter pti(*pc, lev); pti.isValid(); ++pti)
        {
          auto& particles = pti.GetArrayOfStructs();
          BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

          const int np = particles.size();

          for (int pid=0; pid<np; pid++) {
              // Local array storing interpolated values

              BMXParticleContainer::ParticleType& p = pstruct[pid];
              Real *cell_par = &p.rdata(0);
              Real *p_vals = &p.rdata(realIdx::first_data);
              bmxchem->printCellConcentrations((int)p.id(), p_vals, cell_par);
          }
        }
      }

      for (BMXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int np = particles.size();

        const auto& interp_array = interp_ptr->array(pti);

        const auto& interp_varray = interp_vptr->array(pti);

        const auto& interp_narray = interp_nptr->array(pti);

        int l_cnc_deposition_scheme;
        if (bmx::m_cnc_deposition_scheme == DepositionScheme::one_to_one) 
             l_cnc_deposition_scheme = 0;
        else if (bmx::m_cnc_deposition_scheme == DepositionScheme::trilinear) 
             l_cnc_deposition_scheme = 1;
        else 
           amrex::Abort("Dont know this depsoition scheme in calc_txfr_particle");

        int l_vf_deposition_scheme;
        if (bmx::m_vf_deposition_scheme == DepositionScheme::one_to_one) 
             l_vf_deposition_scheme = 0;
        else if (bmx::m_vf_deposition_scheme == DepositionScheme::trilinear) 
             l_vf_deposition_scheme = 1;
        else 
           amrex::Abort("Dont know this depsoition scheme in calc_txfr_particle");

        int nloop = m_nloop;
        amrex::ParallelFor(np,
            [pstruct,interp_array,interp_varray,interp_narray,plo,dxi,grid_vol,dt,
             nloop,chempar,l_cnc_deposition_scheme,l_vf_deposition_scheme]
            AMREX_GPU_DEVICE (int pid) noexcept
              {
              // Local array storing interpolated values
              GpuArray<Real, interp_ncomp> interp_loc;

              // Array storing volume fraction
              GpuArray<Real, 1> interp_vloc;

              // Array storing number of particles
              GpuArray<Real, 1> interp_nloc;

              BMXParticleContainer::ParticleType& p = pstruct[pid];

              if (l_cnc_deposition_scheme == 0)
              {
                  one_to_one_interp(p.pos(), &interp_loc[0],
                                    interp_array, plo, dxi, interp_ncomp);
              }
              else if (l_cnc_deposition_scheme == 1)
              {
                  trilinear_interp(p.pos(), &interp_loc[0],
                                   interp_array, plo, dxi, interp_ncomp);
              }
              if (l_vf_deposition_scheme == 0)
              {
                  one_to_one_interp(p.pos(), &interp_vloc[0],
                                    interp_varray, plo, dxi, 1);
              }
              else if (l_vf_deposition_scheme == 1)
              {
                  trilinear_interp(p.pos(), &interp_vloc[0],
                                   interp_varray, plo, dxi, 1);
              }
              one_to_one_interp(p.pos(), &interp_nloc[0],
                  interp_narray, plo, dxi, 1);

#ifndef AMREX_USE_GPU
              //std::cout<<"Number of particles in grid cell: "<<interp_nloc[0]<<std::endl;
#endif
#ifdef NEW_CHEM
              Real *cell_par = &p.rdata(0);
              Real *p_vals = &p.rdata(realIdx::first_data);
              if (interp_nloc[0] == 0.0) amrex::Abort("Number of particles is Zero!");
#if 0
              printf("   fluid volume fraction    : %16.8e\n",interp_vloc[0]);
              printf("   grid cell volume         : %16.8e\n",grid_vol);
              printf("   number of particles/cell : %16.8f\n",interp_nloc[0]);
              printf("   fluid volume per particle: %16.8e\n",grid_vol*interp_vloc[0]/interp_nloc[0]);
              printf("   time increment           : %16.8e\n",dt);
#endif
              xferMeshToParticleAndUpdateChem(grid_vol*interp_vloc[0], interp_nloc[0], cell_par,
                                              &interp_loc[0], p_vals, dt, nloop, chempar);
#endif
            });
      } // pti
    } // omp region

    delete interp_ptr;
    delete interp_vptr;
    delete interp_nptr;

  } // lev
  amrex::Print() << "TOTAL PARTICLES "<<nparticles<<std::endl;
}
