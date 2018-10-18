#pragma once

#include "planar3D.hh"

/// \brief Функция строит профиль раскрытия по автомодельному решению
double getInitialOpening(
     const double x,
     const double y,
     const double radius,
     const std::vector<double> zP,
     const std::vector<double> W0
);

/// \brief Функция считывает начальные данные из файлов
bool setInitialData(
    const double bn,
    double &xStar0,
    std::vector< std::vector<double> > &A,
    const std::string pathToBarriersFile,
    std::vector< std::vector<double> > &barriers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
);

/// \brief Функция сохранят начальные параметры в файл формата JSON
void saveInitialDataJSON(
    const std::string filename,
    const double modelingTime,
    const std::size_t meshHeight,
    const std::size_t meshLength,
    const double cellHeight,
    const double cellLength,
    std::vector< std::vector<double> > &fluid,
    std::vector< std::vector<double> > &proppant,
    std::vector< std::vector<double> > &layers,
    const std::string tab = std::string("    ")
);
