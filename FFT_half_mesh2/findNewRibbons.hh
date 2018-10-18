#pragma once

#include "planar3D.hh"

/// \brief Функция переопределяет ранее созданный список граничных элементов
void findNewRibbons(
    const size_t i,
    const size_t j,
    const double d,
    const double dMin,
    const double dCenter,
    std::vector<double> &Wt,
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<Ribbon> &ribbons,
    std::vector<Ribbon> &ribbonsOld,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    const double currentTime,
    const double opening
);
