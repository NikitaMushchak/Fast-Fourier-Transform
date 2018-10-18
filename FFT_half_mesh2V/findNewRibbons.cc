#include <algorithm>

#include "planar3D.hh"
#include "collectDistanceFromVelocity.hh"

/// \details Функция определяет новые граничные элементы из списка активных 
/// элементов, рассчитывающей их скорость и время активации (для учёта утечек 
/// в пласт)
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
){
    const std::vector<size_t> element{i, j};
    
    const bool elementIsAMember = (
        std::find(activeElements.begin(), activeElements.end(), element) 
        != activeElements.end()
    );
        
    if(elementIsAMember && d > dMin && 1 > mesh[i][j].type){
        ribbons.push_back(Ribbon(i, j));
        
        mesh[i][j].type = RIBBON;
        
        if(0 == regime){
            distances[i][j] = 0.25 * sqrt(0.5 * M_PI) * E * opening / Kic;
            distances[i][j] *= distances[i][j];
        }else{
            distances[i][j] = collectDistanceFromVelocity(i, j, Wt[index[i][j]],
                ribbonsOld, velocities);
        }
        
        if(epsilon > distances[i][j]){
            distances[i][j] = epsilon;
        }
    }
    
    if(!elementIsAMember && d >= dCenter){
        activeElements.push_back(std::vector<size_t>{i, j});
        
        mesh[i][j].activationTime = currentTime;
    }
}

