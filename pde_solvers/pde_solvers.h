#pragma once

#define M_PI       3.14159265358979323846
#define M_G        9.81
#define KELVIN_OFFSET 273.15
#ifndef ATMOSPHERIC_PRESSURE 
#define ATMOSPHERIC_PRESSURE 101325.0
#endif
#define TECHNICAL_ATMOSPHERE 98097.0

template <typename ContainerType, typename StreamType>
inline StreamType& send_stream(StreamType& o, const ContainerType& x)
{
    for (typename ContainerType::const_iterator i = x.begin(); i != x.end(); i++) {
        if (i != x.begin())
            o << "; ";
        o << std::fixed << (*i);
    }
    return o;
}

template <typename StreamType>
inline StreamType& operator<<(StreamType& o, const vector<double>& x) {
    return send_stream<vector<double>>(o, x);
}


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
#include "pipe/pipe_hydraulic_pde.h"
#include "pipe/pipe_advection_pde.h"

#include "solvers/diffusion_solver.h" // нужно инклудить после объявления трубы и проч.