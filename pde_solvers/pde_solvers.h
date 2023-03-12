#pragma once

#define M_PI       3.14159265358979323846
#define M_G        9.81
#define KELVIN_OFFSET 273.15
#define ATMOSPHERIC_PRESSURE 101325.0
#define TECHNICAL_ATMOSPHERE 98097.0

#include "core/custom_buffer.h"
#include "core/differential_equation.h"
#include "core/profile_structures.h"

#include "solvers/moc_solver.h"
#include "solvers/ode_solver.h"
#include "solvers/godunov_solver.h"

#include "pipe/oil.h"
#include "pipe/pipe_hydraulic_computations.h"
#include "pipe/pipe_hydraulic_struct.h"
#include "pipe/pipe_hydraulic_pde.h"


