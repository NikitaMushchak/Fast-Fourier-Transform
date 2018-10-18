#include "planar3D.hh"

/// \details Функция определяет начальные граничные элементы из списка 
/// активных элементов с учётом заданной сетки
std::vector<Ribbon> findRibbons(
    std::vector< std::vector<Cell> > &mesh,
    const std::vector< std::vector<size_t> > activeElements,
    std::vector< std::vector<double> > &distances,
    double testDistance,
    const double dMin,
    const double dMax
){
    std::vector<double> zeroVector(mesh[0].size(), 0.);
    
    distances.resize(mesh.size());
    std::fill(distances.begin(), distances.end(), zeroVector);
        
    testDistance += epsilon;
    
    const double testDistanceMin = std::pow(testDistance - dMax, 2);
    const double testDistanceMax = std::pow(testDistance - dMin, 2);
    
    std::vector<Ribbon> ribbons;
    
    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];
        
        const double d = std::pow(mesh[i][j].x, 2)
            + std::pow(mesh[i][j].y, 2);
        
        if(testDistanceMin < d && testDistanceMax > d){
            ribbons.push_back(Ribbon(i, j));
            distances[i][j] = testDistance - std::sqrt(d);
            mesh[i][j].type = RIBBON;
        }
    }
    
    return ribbons;
}
