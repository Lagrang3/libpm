#pragma once

#include "state.hpp"
#include <cassert>

/*
    advance of dt
    X: position
    V: velocity
*/
template <class T>
void advance(T& X, const T& V, const double dt)
{
    assert(X.size() == V.size());
    const size_t N = X.size();
    for (size_t i = 0; i < N; ++i)
        X[i] += dt * V[i];
}

/* Leapfrog kick-drift  */
template <class T>
void leapfrogKDK_evolve(T& particle_mesh, state& s, const double dt)
{
    /* dP/dt = -Gm/a \sum_i \vec x_i/x_i^3  */
    advance(s.P, s.F, s.Gm * dt * 0.5 / s.a);

    /* dX/dt = P/a^2 */
    advance(s.X, s.P, dt / s.a / s.a);

    // update force
    s.enforce_PBC();
    particle_mesh(s.F, s.X);

    /* dP/dt = -Gm/a \sum_i \vec x_i/x_i^3  */
    advance(s.P, s.F, s.Gm * dt * 0.5 / s.a);
}
