#include <algorithm>

#include "planar3D.hh"

/// \details Функция определяет тип производной вдоль обоих орт (с учётом 
/// симметрии) для каждой ячейки сетки по списку активных элементов
void setDerivatives(
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<size_t> > &activeElements
){
    auto testSide = [activeElements](const size_t i, const size_t j) -> bool{
        const std::vector<size_t> element{i, j};
        
        return std::find(activeElements.begin(), activeElements.end(), element)
            != activeElements.end();
    };
    
    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];
        
        bool leftIsActive;
        if(i00 == i){
            leftIsActive = false;
        }else{
            leftIsActive = testSide(i - 1, j);
        }
        
        const bool rightIsActive = testSide(i + 1, j);
        const bool bottomIsActive = testSide(i, j - 1);
        const bool topIsActive = testSide(i, j + 1);
        
        if(leftIsActive && rightIsActive){
            mesh[i][j].xDerivative = REGULAR;
        }else{
            if(leftIsActive){
                mesh[i][j].xDerivative = RIGHT;
            }else{
                if(rightIsActive){
                    mesh[i][j].xDerivative = LEFT;
                }else{
                    mesh[i][j].xDerivative = ALONE;
                }
            }
        }
        
        if(bottomIsActive && topIsActive){
            mesh[i][j].yDerivative = REGULAR;
        }else{
            if(bottomIsActive){
                mesh[i][j].yDerivative = TOP;
            }else{
                if(topIsActive){
                    mesh[i][j].yDerivative = BOTTOM;
                }else{
                    mesh[i][j].yDerivative = ALONE;
                }
            }
        }
    }
}

