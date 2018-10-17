#pragma once

#include "planar3D.hh"

/// \brief Функция расчитывает расстояние от граничных элементов до фронта
double collectDistanceFromVelocity(
    const size_t i,
    const size_t j,
    double opening,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &velocities
);
