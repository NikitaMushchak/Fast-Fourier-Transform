/*!
/// \file planar3D.cc
/// \author Egor Starobinskii, Sergey Khlopin
/// \date 17 Spt 2018


файл с функцией planar3D - главной функцией расчета и вызова функций различных режимов
Данный файл содержит в себе определения основных
классов, используемых в демонстрационной программе
*/

#include <vector>
#include <iostream>
#include <algorithm>

#include "ai.hh"
#include "rhs3D.hh"
#include "initialData.hh"
#include "findRibbons.hh"
#include "findNewRibbons.hh"
#include "setDerivatives.hh"
#include "findActiveElements.hh"
#include "buildInfluenceVector.hh"
#include "buildInfluenceMatrix.hh"
#include "buildPartialInfluenceMatrix.hh"

//#define OPENMP
/*
#ifdef OPENMP
    #include <omp.h>
    int getNumberOfThreads(){
        int counter = 0;
        #pragma omp parallel
        {
            counter = omp_get_num_threads();
        }
        return counter;
    }
#else*/
    int getNumberOfThreads(){
        return 1;
    }
    //#endif

///глобальные переменные, необходимые в разных фрагментах расчётной программы

int regime;

double xStar0;
/*!
double dx шаг вычислений по координате x
\brief шаг по координате х
*/
double dx;
double dy;
double bn;
double Q0;
double E;
double C;
double Kic;
double alpha;
double Amu;
double epsilon;
double dt_step;

size_t i00;
size_t j00;

size_t i01;
size_t j01;

std::string buildVersion;

/// \details Функция возвращает название текущего режима распространения
/// трещины
std::string regimeName(){
    switch(regime){
        case -1:
            return "leak-off dominated";
        case 0:
            return "toughness dominated";
        case 1:
            return "viscosity dominated";
        default:
            return "unknown";
    }
}

/*!
\brief __printLogo__ - функция отображения логотипа программы в консоли
\details Функция выводит текстовое лого программы
*/
void printLogo(){
    std::cout << " ____  _                       _____ ____     ____ _     ___ "
        << std::endl;
    std::cout << "|  _ \\| | __ _ _ __   __ _ _ _|___ /|  _ \\   / ___| |   |_ _|"
        << std::endl;
    std::cout << "| |_) | |/ _` | '_ \\ / _` | '__||_ \\| | | | | |   | |    | | "
        << std::endl;
    std::cout << "|  __/| | (_| | | | | (_| | |  ___) | |_| | | |___| |___ | | "
        << std::endl;
    std::cout << "|_|   |_|\\__,_|_| |_|\\__,_|_| |____/|____/   \\____|_____|___|"
        << std::endl;
    std::cout << std::endl << "Developed by REC \"Gazpromneft-Polytech\"."
        << std::endl << std::endl;
}

/*!
\brief __planar3D__ - основная расчётная функция

\detailed основная расчётная функция, считающая распространение трещины в режиме доминирующей вязкости по планарной модели и сохраняющая данные о результатах расчёта в виде текстовых файлов


\param[in] T1 Время окончания расчета
\param[in] mu Динамический коэффициент вязкости
\param[in] theta
\param[in] ts Массштаб времени (число реальных секунд в расчетном времени)
\param[in] cellSize	Число ячеек, которое приходится на радиус автомодельного решения
\param[in] meshSize	Число областей cellSize в расчетной области \f$meshSize=\frac{numCels}{2*\left( cellSize-1 \right)}\f$
\param[in] pathToBarriersFile Путь (относительный или абсолютный) к расположению файла барьеров
\param[in] pathToInjectionFile Путь (относительный или абсолютный) к расположению файла параметров закачки
\param[in] runningFromGUI
\param[in] saveSteps флаг сохранения значений промежуточных расчетов на каждом временном шаге
*/
int planar3D(
    const double T1,
    const double mu,
    const double theta,
    const double ts,
    const double cellSize,
    const double meshSize,
    const std::string pathToBarriersFile,
    const std::string pathToInjectionFile,
    const bool runningFromGUI,
    const bool saveSteps
){
    if(!runningFromGUI){
        printLogo();

        std::cout << "Build: "  << buildVersion<< "." << std::endl;
        std::cout << "AiLibrary version: " << ai::getVersion() << "."
            << std::endl;
        std::cout << std::endl;
        std::cout << "Incoming parameters:" << std::endl
            << "  Q = " << Q0 << ", mu = " << mu << ", "
            << "n = " << bn << ";" << std::endl
            << "  E\' = " << E << ", " << "C = " << C << ","
            << " Kic = " << Kic << ";" << std::endl
            << "  time = " << T1 << ", ts = " << ts << ";" << std::endl
            << "  mesh = " << meshSize << ", cell = " << cellSize << "."
            << std::endl;
        std::cout << "Initial regime: " << regimeName() << "." << std::endl;
        std::cout << std::endl;
    }

    if(!ai::folderExists("./Results")){
        std::cerr << "Cannot find ./Results, create directory and restart."
            << std::endl;

        return 21;
    }

    if(saveSteps){
        if(!ai::folderExists("./Results/Opening")){
            std::cerr << "Cannot find ./Results/Opening, create directory and "
                << "restart." << std::endl;

            return 21;
        }
    }

    std::vector< std::vector<double> > A;
    std::vector< std::vector<double> > barriers;
    std::vector< std::vector<double> > injection;
    if(
        !setInitialData(bn, xStar0, A, pathToBarriersFile, barriers,
            pathToInjectionFile, injection
        )
    ){
        return 22;
    }

    const double T0 = 1.;
    const double mud = theta * mu;
    const double wn = std::pow(1. / std::pow(ts, bn) * mud / E, 1. / (bn + 2.));
    const double gammaR = 1. / 3. * (1. + bn / (bn + 2.));
    const double dt = 0.0001 * std::pow(5. / cellSize, 2);
    dt_step = dt;                       //сохранили шаг по времени для использования в расчете проппанта (Света)

    C *= std::sqrt(ts) / wn;
    Q0 *= ts / wn;

    size_t injectionIndex = 0;

    double timeToChangeInjection = T1 + 10. * dt;

    if(injection.size() > injectionIndex){
        for(size_t i = 0; i < injection.size(); ++i){
            injection[i][1] *= ts / wn;
        }

        if(injection[injectionIndex][0] < dt){
            Q0 = injection[injectionIndex][1];
            ++injectionIndex;
        }

        if(injection.size() > injectionIndex){
            timeToChangeInjection = injection[injectionIndex][0];
        }
    }





    std::vector<double> zP;
    std::vector<double> openingAtTheStart;
    for(size_t i = 0; i < A.size(); ++i){
        A[i][1] *= std::pow(Q0, 1. / 3.);
        openingAtTheStart.push_back(
            A[i][1] * std::pow(T0, (1. - 2. * bn / (bn + 2.)) / 3.)
        );
        zP.push_back(A[i][0]);
    }

    xStar0 *= std::pow(Q0, 1. / 3.) * std::pow(T0, gammaR);

    double vStar = xStar0 * gammaR * std::pow(T0, gammaR - 1);

    const double axMax = meshSize * xStar0;

    dx = xStar0 / cellSize;
    dy = dx;

    for(size_t i = 0; i < barriers.size(); ++i){
        barriers[i][0] *= dx;
        barriers[i][0] += ai::sign(barriers[i][0]) * (xStar0 + epsilon);

        barriers[i][1] *= dx;
        barriers[i][1] += ai::sign(barriers[i][0]) * (xStar0 + epsilon);

        barriers[i][2] *= std::pow(10, 6) / (wn * E);
    }


    std::vector<double> InfluenceVector1;
    std::vector<double> InfluenceVector2;
        std::cout<<"Building Influence Vector";
        buildInfluenceVector(InfluenceVector1,
        InfluenceVector2,
          axMax,
           zP,
          openingAtTheStart,
          barriers);
        std::cout<<".... DONE"<<std::endl;




    std::vector<double> x;
    std::vector<double> y;
    for(double i = 0; i <= axMax + 0.5 * dx + epsilon; i += dx){
        x.push_back(i);
    }
    for(double i = -axMax -0.5 * dy; i <= axMax + 0.5 * dy+ epsilon; i += dy){
        y.push_back(i);
    }

    i00 = 0;
    j00 = floor(0.5 * y.size());

    i01 = i00;
    j01 = j00 - 1 ;

    const size_t xSize = x.size();
    const size_t ySize = y.size();

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
    std::cout << "Building influence matrix... ";

    std::vector<
        std::vector<
            std::vector<
                std::vector<double>
            >
        >
    > influenceMatrix;
    buildInfluenceMatrix(influenceMatrix, xSize, ySize);

    std::cout << "OK............" << std::endl;

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

    const double dMin1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    const double dCenter1 = dx;
    const double dMax1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(0.5 * dy, 2));
    const double dMin2 = std::sqrt(std::pow(1.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    const double dCenter2 = std::sqrt(2) * dx;

    std::cout << "Mesh: " << mesh.size() << "x" << mesh[0].size() << "."
        << std::endl;

    std::vector< std::vector<size_t> >
        activeElements = findActiveElements(mesh, xStar0);

    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;

    setDerivatives(mesh, activeElements);

    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];


        if(
			RIBBON != mesh[i][j].type
            && REGULAR == mesh[i][j].xDerivative
            && REGULAR == mesh[i][j].yDerivative
        ){
            mesh[i][j].type = CHANNEL;
        }

        if(
			RIBBON != mesh[i][j].type
            && LEFT == mesh[i][j].xDerivative
            && REGULAR == mesh[i][j].yDerivative
        ){
            mesh[i][j].type = CHANNEL;
        }

        if(OUTSIDE < mesh[i][j].type){
            mesh[i][j].activationTime = T0 - dt;
        }
    }

    std::vector< std::vector<double> > distances;
    std::vector<Ribbon> ribbons = findRibbons(mesh, activeElements,
        distances, xStar0, dMax1, 2. * dx);

    std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;

    std::vector<double> WkNew;
    std::vector<double> dSigmaCol;
    std::vector< std::vector<double> > partialInfluenceMatrix;

    std::cout << "Building partial influence matrix... ";

    buildPartialInfluenceMatrix(newInfluenceMatrix, activeElements, Wk,
        sigmaCol, WkNew, partialInfluenceMatrix, dSigmaCol, index);

    std::cout << "OK." << std::endl;

    double height = 0.;
    double length = 0.;
    std::vector< std::vector<double> > fracture;


    size_t stepToCheck = round(dx / (dt * vStar) / 20.);
    size_t step = 0;

    double T = T0;

    bool meshIsNotExhausted = true;

    std::vector<double> rhsCE;
    std::vector<double> rhsCEs;
    std::vector<double> pressure;

    std::vector< std::vector<double> > velocities = zeroMatrixXY;
    std::vector< std::vector<double> > openingAtTheStep = zeroMatrixXY;

	std::vector<double> concentration = zeroVectorXY;			//Вектор концентраций проппанта Света
	std::vector<double>  concentration_temp = zeroVectorXY;		//Вектор новых концентраций проппанта Света
	std::vector< std::vector<double> > concentrationAtTheStep = zeroMatrixXY;

    // \todo: указывать число потоков
    std::cout << std::endl << "Running threads... OK. Size: "
        << getNumberOfThreads() << "." << std::endl;

    std::cout << std::endl << "Starting calculations..." << std::endl;

    if(runningFromGUI){
        std::cout << "Progress: 0.0" << std::endl;
    }else{
        ai::showProgressBar(0.);
    }

    auto startTime = ai::time();

    while(T1 >= T && meshIsNotExhausted){
        // calculatePressure(pressure, index, activeElements,
        //     partialInfluenceMatrix, WkNew, dSigmaCol, Wk.size());

            calculatePressureF(pressure, InfluenceVector1, InfluenceVector2, index, activeElements,
             WkNew, dSigmaCol, Wk.size());
//         		calculateOpeningSpeedProp(Wk, rhsCE, rhsCEs, pressure, mesh, index,
// 			activeElements, concentration, concentration_temp, T);

        calculateOpeningSpeed(Wk, rhsCE, rhsCEs, pressure, mesh, index,
            activeElements, T);

        if(0 == regime){
            for(size_t i = 0; i < ribbons.size(); ++i){
                const size_t iRibbon = ribbons[i].i;
                const size_t jRibbon = ribbons[i].j;

                distances[iRibbon][jRibbon] = 0.25 *  sqrt(0.5 * M_PI) * E
                    * Wk[index[iRibbon][jRibbon]] / Kic;
                distances[iRibbon][jRibbon] *= distances[iRibbon][jRibbon];
                if(epsilon > distances[iRibbon][jRibbon]){
                    distances[iRibbon][jRibbon] = epsilon;
                }
            }
        }else{
            calculateVelocity(velocities, index, ribbons, distances, Wk);

            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    distances[i][j] += velocities[i][j] * dt;
                }
            }
        }

        for(size_t i = 0; i < rhsCE.size(); ++i){
            Wk[i] = ai::max(Wk[i] + rhsCE[i] * dt, 0.);
        }

        for(size_t i = 0; i < rhsCEs.size(); ++i){
            WkNew[i] = ai::max(WkNew[i] + rhsCEs[i] * dt, 0.);
        }

        const size_t savedSize = activeElements.size();

        if(step == stepToCheck){
            if(runningFromGUI){
                std::cout << "Progress: " << (T - T0) / (T1 - T0) << std::endl;
            }else{
                ai::showProgressBar((T - T0) / (T1 - T0));
            }

            step = 0;
            if(0 != regime){
                stepToCheck = std::round(dx / (dt * ai::max(velocities)) / 20.);
            }

            std::vector<Ribbon> oldRibbons = ribbons;

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(1 == j || ySize - 1 == j || xSize - 1 == i){
                    meshIsNotExhausted = false;
                    break;
                }

                const double d = distances[i][j];

                if(0 < i){
                    findNewRibbons(i - 1, j, d, dMin1, dCenter1, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j]]);

                    findNewRibbons(i - 1, j - 1, d, dMin2, dCenter2, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j - 1]]);

                    findNewRibbons(i - 1, j + 1, d, dMin2, dCenter2, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j + 1]]);
                }

                findNewRibbons(i + 1, j, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j]]);

                findNewRibbons(i, j - 1, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i][j - 1]]);

                findNewRibbons(i, j + 1, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i][j + 1]]);

                findNewRibbons(i + 1, j - 1, d, dMin2, dCenter2, Wk, mesh,
                    index, activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j - 1]]);

                findNewRibbons(i + 1, j + 1, d, dMin2, dCenter2, Wk, mesh,
                    index, activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j + 1]]);

                if(dMin2 < d){
                    mesh[i][j].type = CHANNEL;
                }
            }

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(1 != mesh[i][j].type){
                    ribbons.erase(ribbons.begin() + k);
                    distances[i][j] = 0;
                }
            }

            if(savedSize != activeElements.size() && meshIsNotExhausted){
                setDerivatives(mesh, activeElements);

                buildPartialInfluenceMatrix(newInfluenceMatrix, activeElements,
                    Wk, sigmaCol, WkNew, partialInfluenceMatrix, dSigmaCol,
                    index
                );
            }

            height = 0.;
            length = 0.;
            for(size_t i = 0; i < xSize; ++i){
                if(1 == mesh[i][j00].type){
                    length = distances[i][j00]
                        + std::sqrt(std::pow(mesh[i][j00].x, 2)
                        + std::pow(mesh[i][j00].y, 2));

                    break;
                }
            }
            for(size_t j = 0; j < ySize; ++j){
                if(1 == mesh[i00][j].type){
                    height = distances[i00][j]
                        + std::sqrt(std::pow(mesh[i00][j].x, 2)
                        + std::pow(mesh[i00][j].y, 2));

                    break;
                }
            }
            fracture.push_back(
                std::vector<double>{
                    T, pressure[index[i00][j00]], length, height
                }
            );


			////////////////////////////////////////////////////////////////////////
			////
			////	Костыль!!! убираем несимметрию в матрице раскрытия и давления
			////
			////////////////////////////////////////////////////////////////////////
			//for (size_t k = 0; k < activeElements.size(); ++k) {
			//	const size_t i = activeElements[k][0];
			//	const size_t j = activeElements[k][1];

			//	if (LEFT != mesh[i][j].xDerivative && RIGHT != mesh[i][j].xDerivative	&& BOTTOM != mesh[i][j].yDerivative && TOP != mesh[i][j].yDerivative && ALONE != mesh[i][j].xDerivative)

			//		if (j == j00)
			//		{
			//			pressure[index[i][j]] = pressure[index[0][i + j00]];
			//			Wk[index[i][j]] = Wk[index[0][i + j00]];
			//		}

			//}

            if(saveSteps){
                for(size_t i = 0; i < xSize; ++i){
                    for(size_t j = 0; j < ySize; ++j){
                        openingAtTheStep[i][j] = 1000 * wn * Wk[index[i][j]];
                        concentrationAtTheStep[i][j] = concentration[index[i][j]];
                    }
                }
                ai::saveMatrix(ai::string("./Results/Opening/time=")
                    + ai::string(T), openingAtTheStep);
                ai::saveMatrix(ai::string("./Results/Concentration/time=")
                    + ai::string(T), concentrationAtTheStep);
            }

            if(0 == regime){
                const size_t lastIndex = fracture.size() - 1;

                const double maxDelta = ai::max(fracture[lastIndex][2]
                    - fracture[lastIndex - 1][2], fracture[lastIndex][3]
                    - fracture[lastIndex - 1][3]) / (fracture[lastIndex][0]
                    - fracture[lastIndex - 1][0]);

                if(epsilon < maxDelta){
                    stepToCheck = ai::min(stepToCheck,
                        (size_t) std::round(dx / (dt * maxDelta) / 20.));
                }
            }
        }

        T += dt;
        ++step;

        if(T > timeToChangeInjection){
            Q0 = injection[injectionIndex][1];
            ++injectionIndex;

            if(injection.size() > injectionIndex){
                timeToChangeInjection = injection[injectionIndex][0];
            }else{
                timeToChangeInjection = T1 + 10 * dt;
            }

            std::string line("Time = ");
            line += ai::string(injection[injectionIndex - 1][0])
                + ai::string(". Injection changed.");
            ai::printLine(line);
        }
    }

    auto finishTime = ai::time();

    if(meshIsNotExhausted){
        if(runningFromGUI){
            std::cout << "Progress: 1.0" << std::endl;
        }else{
            ai::showProgressBar(1.);
        }

        std::cout << std::endl;
    }else{
        if(runningFromGUI){
            std::cout << "Progress: " << (T - T0) / (T1 - T0) << std::endl;
        }else{
            ai::showProgressBar((T - T0) / (T1 - T0));
        }

        std::cout << std::endl;
        std::cout << "Attention. Mesh is exhausted! "
            << "Cannot continue calculation." << std::endl;
        std::cout << "Finished at time = " << T << "." << std::endl;
    }

    std::cout << "Time used: " << ai::duration(startTime, finishTime, "s")
        << "s" << std::endl;

    std::vector< std::vector<double> > openingAtTheEnd = zeroMatrixXY;
    std::vector< std::vector<double> > concentrationAtTheEnd = zeroMatrixXY;
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            openingAtTheEnd[i][j] = 1000 * wn * Wk[index[i][j]];
            concentrationAtTheEnd[i][j] = concentration[index[i][j]];
        }
    }

    std::cout << "Saving results... ";
    ai::saveMatrix("./Results/opening", openingAtTheEnd);
    ai::saveMatrix("./Results/concentration", concentrationAtTheEnd);

    height = 0.;
    length = 0.;
    for(size_t i = 0; i < xSize; ++i){
        if(1 == mesh[i][j00].type){
            length = distances[i][j00]
                + std::sqrt(std::pow(mesh[i][j00].x, 2)
                + std::pow(mesh[i][j00].y, 2));

            break;
        }
    }
    for(size_t j = 0; j < ySize; ++j){
        ///TODO: исправить типы
        if(1 == mesh[i00][j].type){
            height = distances[i00][j]
                + std::sqrt(std::pow(mesh[i00][j].x, 2)
                + std::pow(mesh[i00][j].y, 2));

            break;
        }
    }
    fracture.push_back(
        std::vector<double>{T, pressure[index[i00][j00]], length, height}
    );
    std::stringstream comment;
    comment << std::right << std::setw(14) << "Time\t"
        << std::right << std::setw(14) << "Pressure\t"
        << std::right << std::setw(14) << "Length\t"
        << std::right << std::setw(14) << "Height\t";
    ai::saveMatrix("./Results/fracture", fracture);
    comment.str(std::string());

    std::cout << "OK." << std::endl;




    std::vector< std::vector<double> > fluid;
    fluid.push_back(std::vector<double>{
        0.,
        T1,
        Q0 * wn / ts,
        bn,
        mu
    });
    std::vector< std::vector<double> > proppant;
    proppant.push_back(std::vector<double>{
        0.,
        T1,
        0.
    });
    std::vector< std::vector<double> > layers;
    for(size_t i = 0; i < barriers.size(); ++i){
        layers.push_back(std::vector<double>{
            -1. * barriers[i][1],
            -1. * barriers[i][0],
            E,
            C,
            barriers[i][2] * wn * E
        });
        if(i < barriers.size() - 1){
            if(0 > barriers[i][0] * barriers[i + 1][1]){
                layers.push_back(std::vector<double>{
                    -1. * barriers[i][0],
                    -1. * barriers[i + 1][1],
                    E,
                    C,
                    0.
                });
            }
        }
    }

    saveInitialDataJSON(
        "./Results/parameters",
        T1,
        ySize,
        xSize,
        dy,
        dx,
        fluid,
        proppant,
        layers
    );

    return 0;
}
