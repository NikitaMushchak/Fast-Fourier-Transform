#include "ai.hh"
#include "planar3D.hh"

/// \details Функция заполняет матрицу коэффициентами влияния с учётом
/// симметрии вдоль одного из орт
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
){
    const double coeff = -1. / (8. * M_PI);

    std::vector<double> xIntegrationLimits;
    std::vector<double> yIntegrationLimits;
    std::vector<double> xCollocationPoints;
    std::vector<double> yCollocationPoints;

    std::vector< std::vector<double> > A(xSize);
    for(size_t i = 0; i < A.size(); ++i){
        A[i].resize(ySize);
    }

    for(size_t i = 0; i < ySize; ++i){
        yCollocationPoints.push_back((double) i + 1. - j00);
    }
    std::cout <<"yCollocation"<<std::endl;
    ai::printVector(yCollocationPoints);
    for(size_t i =0; i < yCollocationPoints.size() ;++i)
    {
         xCollocationPoints.push_back(yCollocationPoints[i]);
    }

    for(size_t i = 0; i <= ySize; ++i){
        yIntegrationLimits.push_back((double) i - j00 + 0.5);
    }
    std::cout<<"integr. limits"<<std::endl;
    ai::printVector(yIntegrationLimits);

    for(size_t i =0;i<yIntegrationLimits.size();++i)
    {
        xIntegrationLimits.push_back(yIntegrationLimits[i]);
    }
    influenceMatrix.resize(xSize);
    for(size_t i = 0; i < influenceMatrix.size(); ++i){
        influenceMatrix[i].resize(ySize);

        for(size_t j = 0; j < influenceMatrix[i].size(); ++j){
            influenceMatrix[i][j].resize(xSize);

            for(size_t k = 0; k < influenceMatrix[i][j].size(); ++k){
                influenceMatrix[i][j][k].resize(ySize);
            }
        }
    }

    for(size_t k = 0; k < xSize; ++k){ //
        for(size_t m = 0; m < ySize; ++m){ //
            const double x1 = xCollocationPoints[k];
            const double x2 = yCollocationPoints[m];

            for(size_t i = 0; i < xSize; ++i){//
                for(size_t j = 0; j < ySize; ++j){//
                    const double a = xIntegrationLimits[i];
                    const double b = xIntegrationLimits[i + 1];
                    const double c = yIntegrationLimits[j];
                    const double d = yIntegrationLimits[j + 1];

                    const double term1 = std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - c, 2)) / (x1 - b)
                        - std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - c, 2)) / (x1 - a);
                    const double term2 = std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - d, 2)) / (x1 - a)
                        - std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - d, 2)) / (x1 - b);
                    const double term3 = std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(c - x2, 2)) / (x1 + b)
                        - std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(c - x2, 2)) / (x1 + a);
                    const double term4 = std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(d - x2, 2)) / (x1 + a)
                        - std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(d - x2, 2)) / (x1 + b);

                    A[i][j] = coeff * (term1 / (x2 - c) + term2 / (x2 - d));
                        // + term3 / (c - x2) + term4 / (d - x2));
                }
            }

            influenceMatrix[k][m] = A;
        }
    }
}
//Unchecked!!!!!
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
){
    const double coeff = -1. / (8. * M_PI);

    std::vector<double> xIntegrationLimits;
    std::vector<double> yIntegrationLimits;
    std::vector<double> xCollocationPoints;
    std::vector<double> yCollocationPoints;

    std::vector< std::vector<double> > A(xSize);
    for(size_t i = 0; i < A.size(); ++i){
        A[i].resize(ySize);
    }
    //xCollocationPoints.push_back(i00 -0.5);
    for(size_t i = 0; i < xSize; ++i){
        xCollocationPoints.push_back((double) i - i00);
    }
    xCollocationPoints[0] = 0.25;

    for(size_t i = 0; i < ySize; ++i){
        yCollocationPoints.push_back((double) i + 1. - j00);
    }

    for(size_t i = 0; i <= xSize; ++i){
        xIntegrationLimits.push_back((double) i - i00 - 0.5);
    }
    xIntegrationLimits[0] = 0;

    for(size_t i = 0; i <= ySize; ++i){
        yIntegrationLimits.push_back((double) i - j00 + 0.5);
    }
    std::cout<<"y Intergat small"<<std::endl;
    ai::printVector(yIntegrationLimits);

    influenceMatrix.resize(xSize);
    for(size_t i = 0; i < influenceMatrix.size(); ++i){
        influenceMatrix[i].resize(ySize);

        for(size_t j = 0; j < influenceMatrix[i].size(); ++j){
            influenceMatrix[i][j].resize(xSize);

            for(size_t k = 0; k < influenceMatrix[i][j].size(); ++k){
                influenceMatrix[i][j][k].resize(ySize);
            }
        }
    }

    for(size_t k = 0; k < xSize; ++k){
        for(size_t m = 0; m < ySize; ++m){
            const double x1 = xCollocationPoints[k];
            const double x2 = yCollocationPoints[m];

            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    const double a = xIntegrationLimits[i];
                    const double b = xIntegrationLimits[i + 1];
                    const double c = yIntegrationLimits[j];
                    const double d = yIntegrationLimits[j + 1];

                    const double term1 = std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - c, 2)) / (x1 - b)
                        - std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - c, 2)) / (x1 - a);
                    const double term2 = std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - d, 2)) / (x1 - a)
                        - std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - d, 2)) / (x1 - b);
                    const double term3 = std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(c - x2, 2)) / (x1 + b)
                        - std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(c - x2, 2)) / (x1 + a);
                    const double term4 = std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(d - x2, 2)) / (x1 + a)
                        - std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(d - x2, 2)) / (x1 + b);

                    A[i][j] = coeff * (term1 / (x2 - c) + term2 / (x2 - d)
                        + term3 / (c - x2) + term4 / (d - x2));
                }
            }

            influenceMatrix[k][m] = A;
        }
    }
}
