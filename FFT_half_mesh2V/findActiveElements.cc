#include "planar3D.hh"

/// \details Функция определяет активные элементы в пределах указанного 
/// расстояния, опираясь на заданную сетку
std::vector< std::vector<size_t> > findActiveElements(
	std::vector< std::vector<Cell> > &mesh,
    double testDistance
){
    const size_t xSize = mesh.size();
    const size_t ySize = mesh[0].size();
    
    testDistance = std::pow(testDistance, 2);
    
    std::vector< std::vector<size_t> > activeElements;
    
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            if(std::pow(std::abs(mesh[i][j].x) + 0.25 * dx, 2)
                + std::pow(std::abs(mesh[i][j].y) + 0.25 * dy, 2)
                < testDistance
            ){
                activeElements.push_back(std::vector<size_t>{i, j});
            }
        }
    }
    
    return activeElements;
}
