#pragma once

#define M_PI       3.14159265358979323846
#define M_G        9.81
#define KELVIN_OFFSET 273.15
#ifndef ATMOSPHERIC_PRESSURE 
#define ATMOSPHERIC_PRESSURE 101325.0
#endif
#define TECHNICAL_ATMOSPHERE 98097.0



#include <pde_solvers/timeseries.h>

#include "core/ring_buffer.h"
#include "core/differential_equation.h"
#include "core/profile_structures.h"

#include "solvers/moc_solver.h"
#include "solvers/ode_solver.h"
#include "solvers/godunov_solver.h"
#include "solvers/quick_solver.h"


#include "pipe/oil.h"
#include "pipe/pipe_hydraulic_computations.h"
#include "pipe/pipe_hydraulic_struct.h"

#include "pipe/pipe_advection_pde.h"
#include "pipe/pipe_advection_solver.h"
#include "pipe/pipe_hydraulic_pde.h"
#include "pipe/pipe_condensate_pde.h"
#include "pipe/pipe_profile_utils.h"

#include "pipe/heat_transfer.h"
#include "pipe/pipe_heat_struct.h"
#include "pipe/pipe_dynamic_soil.h"
#include "pipe/pipe_dynamic_soil_multizone.h"
#include "pipe/pipe_heat_computations.h"
#include "pipe/pipe_heat_pde.h"
#include "pipe/pipe_heat_util.h"



#include "solvers/diffusion_solver.h" // нужно инклудить после объявления трубы и проч.
#include "tasks/pipe_heat_task_util.h"
#include "tasks/endogenous_values.h"
#include "tasks/pipe_solver_interface.h"

#include "tasks/isothermal_quasistatic_task.h"
#include "tasks/isothermal_quasistatic_ident.h"
#include "tasks/nonisothermal_quasistatic_task.h"
#include "tasks/nonisothermal_quasistatic_task_p.h"
#include "tasks/nonisothermal_quasistatic_ident.h"
#include "tasks/condensate_pipe_task.h"

namespace pde_solvers
{
;

}


