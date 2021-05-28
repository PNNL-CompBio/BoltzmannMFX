#include <bmx_init_fluid.H>

#include <bmx_calc_fluid_coeffs.H>
#include <bmx_calc_cell.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;

// Forward declarations
void set_ic_chem_species (const Box& sbx, const Box& domain,
                       const Real dx, const Real dy, const Real dz,
                       const GpuArray<Real, 3>& plo, FArrayBox& X_k_fab);

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid                                              !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

void init_fluid (const Box& sbx,
                 const Box& bx,
                 const Box& domain,
                 const MFIter& mfi,
                 LevelData& ld,
                 const Real dx,
                 const Real dy,
                 const Real dz,
                 const Real xlength,
                 const Real ylength,
                 const Real zlength,
                 const GpuArray<Real, 3>& plo,
                 const int advect_fluid_chem_species)
{
  // Fluid SPECIES Initialization
  if (advect_fluid_chem_species) {
    // Set the initial fluid chem_species mass fractions
    set_ic_chem_species(sbx, domain, dx, dy, dz, plo, (*ld.X_k)[mfi]);
  }

  // Initialize all the fluid and fluid chem_species parameters
  init_fluid_parameters(bx, domain, dx, dy, dz, mfi, ld, advect_fluid_chem_species);
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid_parameters                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void init_fluid_parameters (const Box& bx,
                            const Box& domain,
                            const Real dx,
                            const Real dy,
                            const Real dz,
                            const MFIter& mfi,
                            LevelData& ld,
                            const int advect_fluid_chem_species)
{
  // Initialize D_k
  if (advect_fluid_chem_species)  
      calc_D_k(bx, domain, dx, dy, dz, (*ld.D_k)[mfi]);

  Gpu::synchronize();
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid chem_species mass fractions initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_chem_species (const Box& sbx,
                       const Box& domain,
                       const Real dx,
                       const Real dy,
                       const Real dz,
                       const GpuArray<Real, 3>& plo,
                       FArrayBox& X_k_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& X_k = X_k_fab.array();

#ifdef NEW_CHEM
  const RealVect chem(&FLUID::init_conc[0]);
#endif

  const int nchem_species = X_k_fab.nComp();

  amrex::Print() << "SETTING INITIAL CONDITIONS FOR SPECIES " << std::endl;

  ParmParse pp("geometry");
  amrex::Vector<Real> bmin, bmax;
  pp.getarr("prob_lo",bmin);
  pp.getarr("prob_hi",bmax);

  // We set the coeffs at cell centers to D_k in the lower region and 0 above zhi
  Real zhi = FLUID::surface_location;

  ParallelFor(sbx, nchem_species, [=]
       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { 
          Real x   = (i+.5)*dx+bmin[0];
          Real y   = (j+.5)*dy+bmin[1];
          Real z   = (k+.5)*dz+bmin[2];

#if 0
          Real ra = 0.1;
          Real rb = 0.2;

          Real rsq = (x-.5)*(x-.5) +  (y-.5)*(y-.5) +  (z-.5)*(z-.5); 
          // non-zero only within fixed radius
          if (n == 0)
          {
              if (rsq < ra*ra) 
                  X_k(i,j,k,n) = 1.0;
              else
                  X_k(i,j,k,n) = 0.0;
          } else {
              if (rsq < rb*rb) 
                  X_k(i,j,k,n) = 2.0;
              else
                  X_k(i,j,k,n) = 0.0;
          }
#else
         // two layers
         if (z > zhi)
            X_k(i,j,k,n) = 0.0;
         else
            X_k(i,j,k,n) = chem[n];
#endif
      }); 
}
