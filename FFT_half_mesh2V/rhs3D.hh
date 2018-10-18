#pragma once

#include "planar3D.hh"

/// \brief Функция умножает матрицу на вектор
inline void multiply(
    std::vector< std::vector<double> > &matrix,
    std::vector<double> &vector,
    std::vector<double> &result
);

/// \brief Функция расчёта давлений
void calculatePressure(
    std::vector<double> &pressure,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector<double> &Ws,
    std::vector<double> &sigma,
    const size_t size
);
void calculatePressureF(
    std::vector<double> &pressure,
     std::vector<double> &InfluenceVector1,
     std::vector<double> &InfluenceVector2,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<double> &opening,
    std::vector<double> &sigma,
    const size_t size
);
/// \brief Функция расчёта скорости фронта
void calculateVelocity(
    std::vector< std::vector<double> > &velocities,
    std::vector< std::vector<size_t> > &index,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &distances,
    std::vector<double> &Wt
);

/// \brief Функция расчёта производной раскрытия по времени
void calculateOpeningSpeed(
   std::vector<double> &Wt,
   std::vector<double> &dWdt,
   std::vector<double> &dWdts,
   std::vector<double> &pressure,
   std::vector< std::vector<Cell> > &mesh,
   std::vector< std::vector<size_t> > &index,
   std::vector< std::vector<size_t> > &activeElements,
double currentTime
);

/*/
* функция, рассчитывающая производную раскрытия по времени с учётом давления и
* раскрытия в соседних элементах, а также утечки жидкости в пласт, а также переноса проппанта
* TVD схема тоже явная
/*/

void calculateOpeningSpeedProp_TVD(
	std::vector<double> &Wt,
	std::vector<double> &dWdt,
	std::vector<double> &dWdts,
	std::vector<double> &pressure,
	std::vector< std::vector<Cell> > &mesh,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<size_t> > &activeElements,
	std::vector<double> &concentration,
	std::vector<double> &concentration_temp,
	double currentTime
);

/*/
* объявление функции, рассчитывающей производную раскрытия по времени с учётом
* давления и раскрытия в соседних элементах, утечки жидкости в пласт, а также переноса проппанта
/*/

void calculateOpeningSpeedProp(
	std::vector<double> &Wt,
	std::vector<double> &dWdt,
	std::vector<double> &dWdts,
	std::vector<double> &pressure,
	std::vector< std::vector<Cell> > &mesh,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<size_t> > &activeElements,
	std::vector<double> &concentration,
	std::vector<double> &concentration_temp,
	double currentTime
);
