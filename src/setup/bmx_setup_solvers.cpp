//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>

#include <bmx_diffusion_op.H>
#include <bmx_bc_parms.H>

void
bmx::bmx_init_solvers ()
{
    BL_PROFILE("bmx::bmx_init_solvers");

    diffusion_op.reset(new DiffusionOp(this, BC::diff_chem_species_lobc, BC::diff_chem_species_hibc, nghost));
}

void
bmx::bmx_setup_solvers ()
{
    BL_PROFILE("bmx::bmx_setup_solvers");

    diffusion_op->setup(this);
}
