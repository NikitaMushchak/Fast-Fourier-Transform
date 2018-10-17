#include "buildInfluenceMatrix.hh"
#include "ai.hh"
#include "rhs3D.hh"
#include "initialData.hh"
#include "findRibbons.hh"
#include "findNewRibbons.hh"
#include "setDerivatives.hh"
#include "findActiveElements.hh"
#include "buildInfluenceMatrix.hh"
#include "buildPartialInfluenceMatrix.hh"
#include "Circulant_tools2.hh"





void buildInfluenceVector(std::vector<std::vector<double>> &InfluenceVector,
                                            const double axMax,
                                          std::vector<double>& zP,
                                          std::vector<double> &openingAtTheStart,
                                        std::vector< std::vector<double> > &barriers)
{
std::vector<double> x;
std::vector<double> y;
// for(double i = 0; i <= axMax + epsilon; i += dx){
//     x.push_back(i);
// }
for(double i = -axMax - 0.5 * dx; i <= axMax + 0.5 * dx + epsilon; i += dx){
    x.push_back(i);
}

for(double i = -axMax - 0.5 * dy; i <= axMax + 0.5 * dy + epsilon; i += dy){
    y.push_back(i);
}

// i00 = 0;
i00 = floor(0.5 * x.size());
j00 = floor(0.5 * y.size());

const size_t xSize = x.size();
const size_t ySize = y.size();
std::cout<<"xSize= "<<xSize<<"ySize = "<<ySize<<std::endl;
std::vector< std::vector<Cell> > mesh(xSize);
for(size_t i = 0; i < xSize; ++i){
    mesh[i].resize(ySize);

    for(size_t j = 0; j < ySize; ++j){
        mesh[i][j].setCoordinates(x[i], y[j]);
    }
}

std::vector<double> zeroVectorX(xSize, 0.);
std::vector<double> zeroVectorY(ySize, 0.);
std::vector<size_t> zeroSizeTVectorY(ySize, 0);
std::vector<double> zeroVectorXY(xSize * ySize, 0.);
std::vector< std::vector<double> > zeroMatrixXY(xSize, zeroVectorY);

std::cout << std::endl;
std::cout << "Building big influence matrix... ";

std::vector<
    std::vector<
        std::vector<
            std::vector<double>
        >
    >
> influenceMatrix;
buildInfluenceMatrixBig(influenceMatrix, xSize, ySize);

std::cout<<"Inf matrix big size = "<<influenceMatrix.size()<<std::endl;
std::cout << "OK." << std::endl;

std::vector<double> Wk = zeroVectorXY;
std::vector<double> sigmaCol = zeroVectorXY;
std::vector< std::vector<size_t> > index(xSize, zeroSizeTVectorY);
std::vector< std::vector<double> > newInfluenceMatrix(xSize * ySize,
    zeroVectorXY);

size_t p = 0;
size_t m = 0;
for(size_t i = 0; i < xSize; ++i){
    for(size_t j = 0; j < ySize; ++j){
        p = 0;

        for(size_t k = 0; k < xSize; ++k){
            for(size_t l = 0; l < ySize; ++l){
                newInfluenceMatrix[m][p] = influenceMatrix[i][j][k][l];
                ++p;
            }
        }

        index[i][j] = m;
        Wk[m] = getInitialOpening(mesh[i][j].x, mesh[i][j].y, xStar0, zP,
            openingAtTheStart);
        ai::assignFromVectorByIntervalCondition(sigmaCol[m], y[j],
            barriers);

        ++m;
    }
}
// ai::saveMatrix("./infbig",newInfluenceMatrix);
std::vector<double> tempInfluenceVector;

//simC = circulant([newInfluenceMatrix(M:N,M); newInfluenceMatrix(1:M-1,M)],1);

//const size_t middle = (size_t) floor(0.5 * (newInfluenceMatrix.size() - 1));
//size_t addition = 0;
size_t middle = (size_t) floor(0.5 * newInfluenceMatrix.size())  -1- floor(0.5 * xSize);
//middle = 2047-32;
std::cout<<"midlle = "<< middle<<std::endl;
for(size_t i = middle; i < newInfluenceMatrix.size(); ++i){
   tempInfluenceVector.push_back(newInfluenceMatrix[i][middle]);
}
for(size_t i = 0; i < middle; ++i){
   tempInfluenceVector.push_back(newInfluenceMatrix[i][middle]);
}

//ai::saveVector("./inf_vector", tempInfluenceVector);

//std::vector<std::vector<double>> InfluenceVector;
InfluenceVector.resize(newInfluenceMatrix.size());
for(size_t i =0 ; i<newInfluenceMatrix.size();++i)
{
InfluenceVector[i].resize(2);
InfluenceVector[i][0]=tempInfluenceVector[i];
}
// ai::saveMatrix("./infvecCom",InfluenceVector);
fft2(InfluenceVector);
std::cout<<"Inf vector size = "<<InfluenceVector.size()<<std::endl;
newInfluenceMatrix.clear();
tempInfluenceVector.clear();


}
