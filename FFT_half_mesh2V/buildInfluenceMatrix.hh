#pragma once

#include "planar3D.hh"

/// \brief Функция заполняет матрицу коэффициентов влияния
void buildInfluenceMatrix(
    std::vector<
        std::vector<
            std::vector<
                std::vector<double>
            >
        >
    > & influenceMatrix,
	const size_t xSize,
	const size_t ySize
);
void buildInfluenceMatrixBig(
    std::vector<
        std::vector<
            std::vector<
                std::vector<double>
            >
        >
    > & influenceMatrix,
	const size_t xSize,
	const size_t ySize
);
