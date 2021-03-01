#include <bmx_init_fluid.H>

#include <bmx_calc_fluid_coeffs.H>
#include <bmx_calc_cell.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>
// #include <bmx_ic_parms.H>

using namespace amrex;

// Forward declarations
void set_ic_chem_species_g (const Box& sbx, const Box& domain,
                       const Real dx, const Real dy, const Real dz,
                       const GpuArray<Real, 3>& plo, FArrayBox& X_gk_fab);

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
    set_ic_chem_species_g(sbx, domain, dx, dy, dz, plo, (*ld.X_gk)[mfi]);
  }

  // Initialize all the fluid and fluid chem_species parameters
  init_fluid_parameters(bx, mfi, ld, advect_fluid_chem_species);
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid_parameters                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void init_fluid_parameters (const Box& bx,
                            const MFIter& mfi,
                            LevelData& ld,
                            const int advect_fluid_chem_species)
{
  // Initialize D_gk
  if (advect_fluid_chem_species) {
    calc_D_gk(bx, (*ld.D_gk)[mfi]);
  }

  Gpu::synchronize();
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid chem_species mass fractions initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_chem_species_g (const Box& sbx,
                       const Box& domain,
                       const Real dx,
                       const Real dy,
                       const Real dz,
                       const GpuArray<Real, 3>& plo,
                       FArrayBox& X_gk_fab)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  Array4<Real> const& X_gk = X_gk_fab.array();

  const int nchem_species_g = X_gk_fab.nComp();

  amrex::Print() << "SETTING INITIAL CONDITIONS FOR SPECIES " << std::endl;

  ParallelFor(sbx, nchem_species_g, [=]
       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { 
          Real x   = (i+.5)*dx;
          Real y   = (j+.5)*dy;
          Real z   = (k+.5)*dz;

          Real ra = 0.1;
          Real rb = 0.2;

          Real rsq = (x-.5)*(x-.5) +  (y-.5)*(y-.5) +  (z-.5)*(z-.5); 
#if 0
          // non-zero only within fixed radius
          if (n == 0)
          {
              if (rsq < ra*ra) 
                  X_gk(i,j,k,n) = 1.0;
              else
                  X_gk(i,j,k,n) = 0.0;
          } else {
              if (rsq < rb*rb) 
                  X_gk(i,j,k,n) = 2.0;
              else
                  X_gk(i,j,k,n) = 0.0;
          }
#else
         // constant value everywhere
         X_gk(i,j,k,n) = 1.0;
#endif
      }); 
}
