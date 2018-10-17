#pragma once

#include "planar3D.hh"

/// \brief Функция определяет тип производной в ячейках сетки
void setDerivatives(
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<size_t> > &activeElements
);
