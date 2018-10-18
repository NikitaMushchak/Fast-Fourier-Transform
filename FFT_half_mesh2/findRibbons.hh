#pragma once

#include "planar3D.hh"

/// \brief Функция определяет список граничных элементов
std::vector<Ribbon> findRibbons(
    std::vector< std::vector<Cell> > &mesh,
    const std::vector< std::vector<size_t> > activeElements,
    std::vector< std::vector<double> > &distances,
    double testDistance,
    const double dMin,
    const double dMax
);
