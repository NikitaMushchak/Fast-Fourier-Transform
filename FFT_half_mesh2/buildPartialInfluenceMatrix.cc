#include "planar3D.hh"

/// \details Функция заполняет матрицу значениями из общей матрицы 
/// коэффициентов влияния, также заполняя значениями укороченные столбцы 
/// раскрытий и напряжений
void buildPartialInfluenceMatrix(
    std::vector< std::vector<double> > &newInfluenceMatrix,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<double> &Wk,
    std::vector<double> &sigmaCol,
    std::vector<double> &WkNew,
    std::vector< std::vector<double> > &partialInfluenceMatrix,
    std::vector<double> &dSigmaCol,
    std::vector< std::vector<size_t> > &index
){
    std::vector<double> zeroVector(activeElements.size(), 0);

    WkNew = zeroVector;
    dSigmaCol = zeroVector;

    partialInfluenceMatrix.resize(activeElements.size());
    std::fill(partialInfluenceMatrix.begin(), partialInfluenceMatrix.end(),
        zeroVector);

    for(size_t j = 0; j < activeElements.size(); ++j){
        const size_t ai = activeElements[j][0];
        const size_t aj = activeElements[j][1];
        const size_t pc1 = index[ai][aj];
        
        WkNew[j] = Wk[pc1];
        dSigmaCol[j] = sigmaCol[pc1];
        
        for(size_t i = 0; i < activeElements.size(); ++i){
            const size_t ai2 = activeElements[i][0];
            const size_t aj2 = activeElements[i][1];
            const size_t pc2 = index[ai2][aj2];
            
            partialInfluenceMatrix[i][j] = newInfluenceMatrix[pc2][pc1];
        }
    }
}
