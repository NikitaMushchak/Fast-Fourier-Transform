/*#ifdef OPENMP
    #include <omp.h>
#endif*/

#include "ai.hh"
#include "planar3D.hh"
 #include "Circulant_tools.hh"
/*!
\detailed __multiply__ функция умножает матрицу на вектор и возвращает результат по ссылке

\param[in] matrix - исходная матрица (N*N)
\param[in] vector - исходный вектор (1*N)
\param[out] result - результирующий вектор (1*N)
*/
inline void multiply(
    std::vector< std::vector<double> > &matrix,
    std::vector<double> &vector,
    std::vector<double> &result
){
    result.resize(matrix.size());

    // #pragma omp parallel for
    for(size_t i = 0; i < matrix.size(); ++i){
        size_t j = 0;

        for(; j <= vector.size() - 16; j += 16){
            result[i] += matrix[i][j] * vector[j]
                + matrix[i][j + 1] * vector[j + 1]
                + matrix[i][j + 2] * vector[j + 2]
                + matrix[i][j + 3] * vector[j + 3]
                + matrix[i][j + 4] * vector[j + 4]
                + matrix[i][j + 5] * vector[j + 5]
                + matrix[i][j + 6] * vector[j + 6]
                + matrix[i][j + 7] * vector[j + 7]
                + matrix[i][j + 8] * vector[j + 8]
                + matrix[i][j + 9] * vector[j + 9]
                + matrix[i][j + 10] * vector[j + 10]
                + matrix[i][j + 11] * vector[j + 11]
                + matrix[i][j + 12] * vector[j + 12]
                + matrix[i][j + 13] * vector[j + 13]
                + matrix[i][j + 14] * vector[j + 14]
                + matrix[i][j + 15] * vector[j + 15];
        }

        for(; j < vector.size(); ++j){
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

/*!
\detailed __calculatePressure__ - функция, рассчитывающая давление в активных элементах с помощью матрицы коэффициентов влияния, раскрытия и заданного контраста напряжений

\param[in] pressure
\param[in] index
\param[in] activeElements
\param[in] influenceMatrix
\param[in] opening
\param[in] sigma
\param[in] size
*/
void calculatePressure(
    std::vector<double> &pressure,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector<double> &opening,
    std::vector<double> &sigma,
    const size_t size
){
    pressure.resize(size);
    std::fill(pressure.begin(), pressure.end(), 0);

    std::vector<double> pressureS;
    multiply(influenceMatrix, opening, pressureS);

    // #pragma omp parallel for
    for(size_t i = 0; i < activeElements.size(); ++i){
        const size_t ai = activeElements[i][0];
        const size_t aj = activeElements[i][1];
        pressure[index[ai][aj]] = pressureS[i] / dx + sigma[i];
    }
}


void calculatePressureF(
    std::vector<double> &pressure,
     std::vector<double> &InfluenceVector1,
     std::vector<double> &InfluenceVector2,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<double> &opening,
    std::vector<double> &sigma,
    const size_t size
){
    // std::cout<<"size = "<<size<<std::endl;
    pressure.resize(InfluenceVector1.size());
    // std::cout<<"....";
    std::fill(pressure.begin(), pressure.end(), 0);
// std::cout<<"calc pressure";
    std::vector<double> pressureS1;
    std::vector<double> pressureS2;
    pressureS1.resize(InfluenceVector1.size());
    pressureS2.resize(InfluenceVector1.size());


    for(size_t i = 0; i < InfluenceVector1.size(); ++i){

        pressureS1[i] = 0.;
        pressureS2[i] = 0.;
    }
    // std::cout<<"pressureS size =" << pressureS.size()<<std::endl;
    std::vector<double> product1;
    std::vector<double> product2;
        product1.resize(pressureS1.size());
        product2.resize(pressureS1.size());
        // for(size_t i = 0 ;i<pressureS.size();++i)
        // {
        //     product[i].resize(2);
        // }

        for(size_t i = 0; i < activeElements.size(); ++i){
        const size_t ai = activeElements[i][0];
        const size_t aj = activeElements[i][1];
            // pressureS[index[ai][aj]][0] = opening[i];
            pressureS1[index[ai][aj]] = opening[i];
    }
    fft2_RealSymm(pressureS1,pressureS2);
    getWiseElement2(InfluenceVector1,InfluenceVector2, pressureS1,pressureS2, product1,product2);

    ifft2_RealSymm(product1, product2);
    // std::cout<<"..Solve..";
    // #pragma omp parallel for
    for(size_t i = 0; i < activeElements.size(); ++i){
        const size_t ai = activeElements[i][0];
        const size_t aj = activeElements[i][1];
       pressure[index[ai][aj]] = product1[index[ai][aj]] / dx + sigma[i];
      // pressure[index[ai][aj]] = a[index[ai][aj]] / dx + sigma[i];
    }
}



/*!
\detailed __calculateVelocity__ - функция, рассчитывающая скорость фронта по заданному режиму развития

\param[in] velocities
\param[in] index
\param[in] ribbons
\param[in] distances
\param[in] Wt
*/
void calculateVelocity(
    std::vector< std::vector<double> > &velocities,
    std::vector< std::vector<size_t> > &index,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &distances,
    std::vector<double> &Wt
){
    velocities.resize(distances.size());
    std::vector<double> zeroVector(distances[0].size(), 0);
    std::fill(velocities.begin(), velocities.end(), zeroVector);

    if(1 == regime){
        //viscosity dominated regime
        const double D3 = 1. / (1. - alpha);

        for(size_t i = 0; i < ribbons.size(); ++i){
            const size_t iRibbon = ribbons[i].i;
            const size_t jRibbon = ribbons[i].j;

            const double distance = distances[iRibbon][jRibbon];

            if(epsilon <= distance){
                velocities[iRibbon][jRibbon] =
                    std::pow(Wt[index[iRibbon][jRibbon]]
                        / (Amu * std::pow(distance, alpha)), D3);
            }else{
                velocities[iRibbon][jRibbon] = 0.;
            }
        }
    }else{
        if(-1 == regime){
            //leak-off dominated regime
            const double alphaL = (bn + 4.) / (4. * bn + 4.);
            const double D3 = 3. / (1. - alphaL);
            const double Ainf = std::pow(Amu, (1. + 2. * alphaL) / 3.);

            for(size_t i = 0; i < ribbons.size(); ++i){
                const size_t iRibbon = ribbons[i].i;
                const size_t jRibbon = ribbons[i].j;

                const double distance = distances[iRibbon][jRibbon];

                if(epsilon <= distance){
                    velocities[iRibbon][jRibbon] =
                        std::pow(Wt[index[iRibbon][jRibbon]]
                            / (Ainf * std::pow(distance, alphaL)), D3)
                        / std::pow(4. * C, 2);
                }else{
                    velocities[iRibbon][jRibbon] = 0.;
                }
            }
        }
    }
}

/*!
\detailed __calculateOpeningSpeed__ - функция, рассчитывающая производную раскрытия по времени с учётом давления и раскрытия в соседних элементах, а также утечки жидкости в пласт без проппанта

\param[in] Wt !!!!тут все нужно скопировать с твд схемы!!!!
\param[in] dWdt
\param[in] dWdts
\param[in] pressure
\param[in] mesh
\param[in] index
\param[in] activeElements
\param[in] currentTime
*/
void calculateOpeningSpeed(
	std::vector<double> &Wt,
	std::vector<double> &dWdt,
	std::vector<double> &dWdts,
	std::vector<double> &pressure,
	std::vector< std::vector<Cell> > &mesh,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<size_t> > &activeElements,
	double currentTime
) {
	dWdt.resize(Wt.size());
	std::fill(dWdt.begin(), dWdt.end(), 0.);

	dWdts.resize(activeElements.size());
	std::fill(dWdts.begin(), dWdts.end(), 0.);

	const double WPow = (2. * bn + 1.) / bn;
	const double PPow = 1. / bn;

	double flowIMinusHalf;
	double flowIPlusHalf;
	double flowJMinusHalf;
	double flowJPlusHalf;

	for (size_t k = 0; k < activeElements.size(); ++k) {
		const size_t i = activeElements[k][0];
		const size_t j = activeElements[k][1];

		const double Pij = pressure[index[i][j]];

		flowIMinusHalf = 0.;
		flowIPlusHalf = 0.;
		flowJMinusHalf = 0.;
		flowJPlusHalf = 0.;

		if (LEFT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIMinusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i - 1][j]]);
			flowIMinusHalf = -ai::sign(pressure[index[i - 1][j]] - Pij)
				* std::pow(openingIMinusHalf, WPow)
				* std::pow(std::abs(pressure[index[i - 1][j]] - Pij) / dx, PPow);
		}

		if (RIGHT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIPlusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i + 1][j]]);
			flowIPlusHalf = ai::sign(pressure[index[i + 1][j]] - Pij)
				* std::pow(openingIPlusHalf, WPow)
				* std::pow(std::abs(pressure[index[i + 1][j]] - Pij) / dx, PPow);
		}

		if (i00 == i) {
			flowIMinusHalf = -flowIPlusHalf;
		}

		if (BOTTOM != mesh[i][j].yDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingJMHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i][j - 1]]);
			flowJMinusHalf = -ai::sign(pressure[index[i][j - 1]] - Pij)
				* std::pow(openingJMHalf, WPow)
				* std::pow(std::abs(pressure[index[i][j - 1]] - Pij) / dx, PPow);
		}

		if (TOP != mesh[i][j].yDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingJPlusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i][j + 1]]);
			flowJPlusHalf = ai::sign(pressure[index[i][j + 1]] - Pij)
				* std::pow(openingJPlusHalf, WPow)
				* std::pow(std::abs(pressure[index[i][j + 1]] - Pij) / dx, PPow);
		}

		const double xDerivative = (flowIPlusHalf - flowIMinusHalf) / dx;
		const double yDerivative = (flowJPlusHalf - flowJMinusHalf) / dy;

		const double value = xDerivative + yDerivative
			+ 0.25 * Q0 * (double)(i == i00 && j == j00) / std::pow(dx, 2)
      + 0.25 * Q0 * (double)(i == i01 && j == j01) / std::pow(dx, 2)
			- C / std::sqrt(currentTime - mesh[i][j].activationTime);

		dWdt[index[i][j]] = value;
		dWdts[k] = value;
	}
}

/*!
\detailed __FluxLimiter__ - функция перерасчета постоков с учетом потоков вторго порядка.

\param[in] type - тип ограничителя потока !!!Дописать значения параметров!!!!
\param[in] r - значение величины потока второго порядка пропанта в ячейке
*/

double FluxLimiter(int type, double r)
{
	double F = 0.;
	switch (type)
	{
	case 0: //'Upwind'
		F = 0.;
		break;
	case 1: //'VanLeer'
		F = (r + std::fabs(r)) / (1 + r * r);
		break;
	case 2: //'minmod'
				F = std::fmax(0, std::fmin(1, r));
		break;
	case 3: //'MC'
				F = std::fmax(0, std::fmin(2 * r, std::fmin(0.5*(1 + r), 2)));
		break;
	case 4: //'Superbee'
				F = std::fmax(0, std::fmax(std::fmin(2 * r, 1), std::fmin(r, 2)));
		break;

	default:
		F = 0.;
		break;
	}
	return F;
}



/*!
\detailed __calculateOpeningSpeedProp_TVD__ - функция, рассчитывающая производную раскрытия по времени с учётом давления и раскрытия в соседних элементах, а также утечки жидкости в пласт, а также переноса проппанта. Расчет TVD схемы ведется явным методом.

\param[in] Wt
\param[in] dWdt
\param[in] dWdts
\param[in] pressure
\param[in] mesh
\param[in] index
\param[in] activeElements
\param[in] concentration -вектор концентраций проппанта
\param[in] concentration_temp - вектор новых концентраций проппанта
\param[in] currentTime
*/
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
){


	dWdt.resize(Wt.size());
    std::fill(dWdt.begin(), dWdt.end(), 0.);

    dWdts.resize(activeElements.size());
    std::fill(dWdts.begin(), dWdts.end(), 0.);

    const double WPow = (2. * bn + 1.) / bn;
    const double PPow = 1. / bn;

    double flowIMinusHalf;		//Егоровские потоки??!!
    double flowIPlusHalf;
    double flowJMinusHalf;
    double flowJPlusHalf;

	////Данные от Светы
	double p_velocityIMinusHalf;						//расчет скорости пропанта в ячейке
	double p_velocityIPlusHalf;
	double p_velocityJMinusHalf;
	double p_velocityJPlusHalf;
	double p_flowIMinusHalf;							//расчет потоков пропанта в ячейке Егор помни, это не твои потоки!!!
	double p_flowIPlusHalf;
	double p_flowJMinusHalf;
	double p_flowJPlusHalf;
	double p_flowIMinusHalf_sdOrder;							//расчет потоков второго порядка пропанта в ячейке  Егор помни, это ваще не твои потоки!!!
	double p_flowIPlusHalf_sdOrder;
	double p_flowJMinusHalf_sdOrder;
	double p_flowJPlusHalf_sdOrder;
	int sflowLeft;									//Знак направления
	int sflowRight;
	int sflowTop;
	int sflowBottom;

	double Rop = 2.65 * 1000;											//плотность проппанта
	double Rof = 1.1*1000;												//плотность жидкости
	double g = 9.81;
	double d = 2 * 0.001; 												//диаметр частицы проппанта
	double betta = -2.5*bn;  											//степень для вязкости
	double C_max = 0.6;  												//максимальная концентрация проппанта
	double alfa = 5;  													//степень для корректирующего коэффициента в уравнении для скорости
	double vs = g * (Rop - Rof)*std::pow(d, (bn + 1.))/(std::pow(3.,(bn - 1))*18.);	// скорость оседания пропанта
	double cp = -2.5;				// у Светы это расчетная величина:	cp = betta / bn;
	int type = 0;														//Тип ограничителя потока, используется в функции FluxLimiter
																		//0: 'Upwind'
																		//1: 'VanLeer'
																		//2: 'minmod'
																		//3: 'MC'
																		//4: 'Superbee'

	double injFunc;														//Скорость закачки
	double gradient_concentration;										//Градиент изменения концентрации. Меняется при наличии закачки в центре терщины

	injFunc = 0.5;

	int minI = 10;
	for (size_t k = 0; k < activeElements.size(); ++k) {
		const size_t i = activeElements[k][0];
		const size_t j = activeElements[k][1];

		const double Pij = pressure[index[i][j]];

		flowIMinusHalf = 0.;
		flowIPlusHalf = 0.;
		flowJMinusHalf = 0.;
		flowJPlusHalf = 0.;

		p_flowIMinusHalf = 0.;
		p_flowIPlusHalf = 0.;
		p_flowJMinusHalf = 0.;
		p_flowJPlusHalf = 0.;

		p_velocityIMinusHalf = 0.;
		p_velocityIPlusHalf = 0.;
		p_velocityJMinusHalf = 0.;
		p_velocityJPlusHalf = 0.;

		gradient_concentration = 0.;

		p_flowIMinusHalf_sdOrder = 0.;
		p_flowIPlusHalf_sdOrder = 0.;
		p_flowJMinusHalf_sdOrder = 0.;
		p_flowJPlusHalf_sdOrder = 0.;

		sflowLeft = 0;
		sflowRight = 0;
		sflowTop = 0;
		sflowBottom = 0;

		if (LEFT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIMinusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i - 1][j]]);
			const double concentrationIMinusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i - 1][j]]);
			const double pressure_gradient = pressure[index[i - 1][j]] - Pij;

			flowIMinusHalf = -ai::sign(pressure_gradient)
				* std::pow(openingIMinusHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationIMinusHalf / C_max), cp);

			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i - 1][j]] >= C_max && pressure_gradient < 0))
				//Нужно уточнить знаки больше или меньше
			{
				p_velocityIMinusHalf = 0;
			}
			else
			{
				p_velocityIMinusHalf = flowIMinusHalf / openingIMinusHalf;
				//расчет новых значений потоков в ячейке в зависимости от направления скорости перетоков между ячейками

				if (p_velocityIMinusHalf < 0)
				{
					p_flowIMinusHalf = p_velocityIMinusHalf * concentration[index[i - 1][j]] * Wt[index[i - 1][j]];
				}
				else
				{
					p_flowIMinusHalf = p_velocityIMinusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}
			//Реализация TVD схемы
			//Левая граница
			//if (concentration[index[i - 1][j]] >epsilon && concentration[index[i][j]] >epsilon)
				if (concentration[index[i - 1][j]] != 0 && concentration[index[i][j]] != 0)
				{
				sflowLeft = ai::sign(p_flowIMinusHalf);
				p_flowIMinusHalf_sdOrder = (Wt[index[i + sflowLeft][j]] * concentration[index[i + sflowLeft][j]] - Wt[index[i - 1 + sflowLeft][j]] * concentration[index[i - 1 + sflowLeft][j]])
					/ (Wt[index[i][j]] * concentration[index[i][j]] - Wt[index[i - 1][j]] * concentration[index[i - 1][j]]);
			}
			//Проводим перерасчет постоков с учетом потоков вторго порядка (функция FluxLimiter)
			p_flowIMinusHalf = p_flowIMinusHalf + 0.5*FluxLimiter(type, p_flowIMinusHalf_sdOrder)*((double)sflowLeft- p_velocityIMinusHalf* dt_step/dx)*p_velocityIMinusHalf
								*(Wt[index[i][j]] * concentration[index[i][j]] - Wt[index[i - 1][j]] * concentration[index[i - 1][j]]);
		}



		if (RIGHT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIPlusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i + 1][j]]);
			const double concentrationIPlusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i + 1][j]]);
			const double pressure_gradient = pressure[index[i + 1][j]] - Pij;

			flowIPlusHalf = ai::sign(pressure_gradient)
				* std::pow(openingIPlusHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationIPlusHalf / C_max), cp);

			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i + 1][j]] >= C_max && pressure_gradient < 0))
			{
				p_velocityIPlusHalf = 0;
			}
			else
			{
				p_velocityIPlusHalf = flowIPlusHalf / openingIPlusHalf;

				//расчет новых значений потоков в ячейке в зависимости от направления скорости перетоков между ячейками
				if (p_velocityIPlusHalf > 0)
				{
					p_flowIPlusHalf = p_velocityIPlusHalf * concentration[index[i + 1][j]] * Wt[index[i + 1][j]];
				}
				else
				{
					p_flowIPlusHalf = p_velocityIPlusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}
			//Реализация TVD схемы
			//Правая граница
			//if (concentration[index[i + 1][j]] >epsilon && concentration[index[i][j]] >epsilon)
				if (concentration[index[i + 1][j]] != 0 && concentration[index[i][j]] != 0)
				{
				sflowRight = ai::sign(p_flowIPlusHalf);
				p_flowIPlusHalf_sdOrder = (Wt[index[i + 1 + sflowRight][j]] * concentration[index[i + 1 + sflowRight][j]] - Wt[index[i + sflowRight][j]] * concentration[index[i + sflowRight][j]])
					/ (Wt[index[i + 1][j]] * concentration[index[i + 1][j]] - Wt[index[i][j]] * concentration[index[i][j]]);
			}
			//Проводим перерасчет постоков с учетом потоков вторго порядка (функция FluxLimiter)
			p_flowIPlusHalf = p_flowIPlusHalf + 0.5*FluxLimiter(type, p_flowIPlusHalf_sdOrder)*((double)sflowRight - p_velocityIPlusHalf * dt_step / dx)*p_velocityIPlusHalf
								*(Wt[index[i + 1][j]] * concentration[index[i + 1][j]] - Wt[index[i][j]] * concentration[index[i][j]]);
		}

		/////////////////////////////////////////////////////////
		//Для нулевого слоя "отражаем значения потоков для симметрии
			if(i00 == i){
		    flowIMinusHalf = -flowIPlusHalf;
			p_flowIMinusHalf = -p_flowIPlusHalf;
			p_velocityIMinusHalf = -p_velocityIPlusHalf;
			}


			if (BOTTOM != mesh[i][j].yDerivative
				&& ALONE != mesh[i][j].xDerivative
				) {
				const double openingJMHalf = 0.5
					* (Wt[index[i][j]] + Wt[index[i][j - 1]]);
				const double concentrationJMinusHalf = 0.5
					* (concentration[index[i][j]] + concentration[index[i][j - 1]]);
				const double pressure_gradient = pressure[index[i][j - 1]] - Pij;

				flowJMinusHalf = -ai::sign(pressure_gradient)
					* std::pow(openingJMHalf, WPow)
					* std::pow(std::abs(pressure_gradient) / dx, PPow)
					/ std::pow((1. - concentrationJMinusHalf / C_max), cp);

				if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i][j - 1]] >= C_max && pressure_gradient < 0))
					//Нужно уточнить знаки больше или меньше
				{
					p_velocityJMinusHalf = 0;

				}
				else
				{
					p_velocityJMinusHalf = flowJMinusHalf / openingJMHalf;

					//Оседание проппанта
			//		p_velocityJMinusHalf = p_velocityJMinusHalf - vs * std::pow((1. - concentrationJMinusHalf / C_max), alfa);
					if (p_velocityJMinusHalf < 0)
					{
						p_flowJMinusHalf = p_velocityJMinusHalf * concentration[index[i][j - 1]] * Wt[index[i][j - 1]];
					}
					else
					{
						p_flowJMinusHalf = p_velocityJMinusHalf * concentration[index[i][j]] * Wt[index[i][j]];
					}
				}
				//Реализация TVD схемы
				//Нижняя граница
				//if (concentration[index[i][j - 1]] >epsilon && concentration[index[i][j]] >epsilon)
					if (concentration[index[i][j - 1]] != 0 && concentration[index[i][j]] != 0)
					{
					sflowBottom = ai::sign(p_flowJMinusHalf);
					p_flowJMinusHalf_sdOrder = (Wt[index[i][j + sflowBottom]] * concentration[index[i][j + sflowBottom]] - Wt[index[i][j - 1 + sflowBottom]] * concentration[index[i][j - 1 + sflowBottom]])
						/ (Wt[index[i][j]] * concentration[index[i][j]] - Wt[index[i][j - 1]] * concentration[index[i][j - 1]]);
				}
				//Проводим перерасчет постоков с учетом потоков вторго порядка (функция FluxLimiter)
				p_flowJMinusHalf = p_flowJMinusHalf + 0.5*FluxLimiter(type, p_flowJMinusHalf_sdOrder)*((double)sflowBottom - p_velocityJMinusHalf * dt_step / dy)*p_velocityJMinusHalf
								*(Wt[index[i][j]] * concentration[index[i][j]] - Wt[index[i][j - 1]] * concentration[index[i][j - 1]]);
			}

			if (TOP != mesh[i][j].yDerivative
				&& ALONE != mesh[i][j].xDerivative
				) {
				const double openingJPlusHalf = 0.5
					* (Wt[index[i][j]] + Wt[index[i][j + 1]]);
				const double concentrationJPlusHalf = 0.5
					* (concentration[index[i][j]] + concentration[index[i][j + 1]]);
				const double pressure_gradient = pressure[index[i][j + 1]] - Pij;

				flowJPlusHalf = ai::sign(pressure_gradient)
					* std::pow(openingJPlusHalf, WPow)
					* std::pow(std::abs(pressure_gradient) / dx, PPow)
					/ std::pow((1. - concentrationJPlusHalf / C_max), cp);

				if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i][j + 1]] >= C_max && pressure_gradient < 0))
				{
					p_velocityJPlusHalf = 0;
				}
				else
				{
					p_velocityJPlusHalf = flowJPlusHalf / openingJPlusHalf;

					//Оседание проппанта
			//		p_velocityJPlusHalf = p_velocityJPlusHalf - vs * std::pow((1. - concentrationJPlusHalf / C_max), alfa);
					if (p_velocityJPlusHalf > 0)
					{
						p_flowJPlusHalf = p_velocityJPlusHalf * concentration[index[i][j + 1]] * Wt[index[i][j + 1]];
					}
					else
					{
						p_flowJPlusHalf = p_velocityJPlusHalf * concentration[index[i][j]] * Wt[index[i][j]];
					}
				}
				//Реализация TVD схемы
				//Верхняя граница
				//if (concentration[index[i][j + 1]] >epsilon && concentration[index[i][j]] >epsilon)
					if (concentration[index[i][j + 1]] != 0 && concentration[index[i][j]] != 0)
					{
					sflowTop = ai::sign(p_flowJPlusHalf);
					p_flowJPlusHalf_sdOrder = (Wt[index[i][j + 1 + sflowTop]] * concentration[index[i][j + 1 + sflowTop]] - Wt[index[i][j + sflowTop]] * concentration[index[i][j + sflowTop]])
						/ (Wt[index[i][j + 1]] * concentration[index[i][j + 1]] - Wt[index[i][j]] * concentration[index[i][j]]);
				}
				//Проводим перерасчет постоков с учетом потоков вторго порядка (функция FluxLimiter)
				p_flowJPlusHalf = p_flowJPlusHalf + 0.5*FluxLimiter(type, p_flowJPlusHalf_sdOrder)*((double)sflowTop - p_velocityJPlusHalf * dt_step / dy)*p_velocityJPlusHalf
									*(Wt[index[i][j + 1]] * concentration[index[i][j + 1]] - Wt[index[i][j]] * concentration[index[i][j]]);

			}

        const double xDerivative = (flowIPlusHalf - flowIMinusHalf) / dx;
        const double yDerivative = (flowJPlusHalf - flowJMinusHalf) / dy;

        const double value = xDerivative + yDerivative
            + Q0 * (double)(i == i00 && j == j00) / std::pow(dx, 2)
            - C / std::sqrt(currentTime-mesh[i][j].activationTime);

		dWdt[index[i][j]] = value;
		dWdts[k] = value;



		//Костыль!!!!!!!!!!!!!!!!!!
		//Проверяем непопадание на границу трещины
		if (value > epsilon && Wt[index[i][j]] > epsilon)
		{
			if (i == i00 && j == j00)
			{
				gradient_concentration = ((((p_flowIPlusHalf - p_flowIMinusHalf) / dx - (p_flowJMinusHalf - p_flowJPlusHalf) / dy) + injFunc) *dt_step);
			}
			else
			{
				gradient_concentration = (((p_flowIPlusHalf - p_flowIMinusHalf) / dx - (p_flowJMinusHalf - p_flowJPlusHalf) / dy) *dt_step);
			}
			concentration_temp[index[i][j]] = (gradient_concentration
				+ concentration[index[i][j]] * Wt[index[i][j]]) / (Wt[index[i][j]] + value * dt_step);
		}
	}
	concentration = concentration_temp;
	return;
}



/*!
\detailed __calculateOpeningSpeedProp__ - функция, рассчитывающая производную раскрытия по времени с учётом давления и раскрытия в соседних элементах,
а также утечки жидкости в пласт, а также переноса проппанта Явная противопоточная разностная схема
<br>Уравнение переноса проппанта:
<br>\f$[\frac{\partial \left( {{}_{p}}w \right)}{\partial t}+\nabla \cdot \left( {{C}_{p}}w{{\mathbf{v}}_{\mathbf{}}} \right)=0,]\f$
<br>где \f${{C}_{p}}\f$– концентрация проппанта, \f$[{{\mathbf{v}}_{\mathbf{}}}]\f$– скорость проппанта.
<br>\f$\frac{u_{i,j}^{k+1}-u_{i,j}^{k}}{\vartriangle t}+\frac{H_{i+0.5,j}^{k}-H_{i-0.5,j}^{k}}{\vartriangle x}+\frac{Q_{i,j+0.5}^{k}-Q_{i,j-0.5}^{k}}{\vartriangle y}=0,\f$
<br>где \f$u_{i,j}^{k}=\left( {{C}_{p}}w \right)_{i,j}^{k}\f$
<br>\f$H_{i+0.5,j}^{k}=h_{i+0.5,j}^{k}+\frac{1}{2}\delta \left( r \right)\left[ sgn \left( \left( {{v}_{p}} \right)_{i+0.5,j}^{k} \right)-\lambda \left( {{v}_{p}} \right)_{i+0.5,j}^{k} \right]\left( {{v}_{p}} \right)_{i+0.5,j}^{k}\left( u_{i+1,j}^{k}-u_{i,j}^{k} \right)\f$
<br>\f$[h_{i+0.5,j}^{k}=\left\{ \begin{align}
  & u_{i,j}^{k}\cdot \left( {{v}_{p}} \right)_{i+0.5,j}^{k},\,\,\,\,\,\left( {{v}_{p}} \right)_{i+0.5,j}^{k}\ge 0 \\
 & u_{i+1,j}^{k}\cdot \left( {{v}_{p}} \right)_{i+0.5,j}^{k},\,\,\left( {{v}_{p}} \right)_{i+0.5,j}^{k}\le 0 \\
\end{align} \right.]\f$
<br>\f$[\delta \left( r \right)]\f$ – ограничитель потока (Flux limiter)
<br>
\f$ [\begin{align}
& \left. 1 \right)\delta \left( r \right)=\min \bmod \left( 1,r \right), \\
& \left. 2 \right)\delta \left( r \right)=\frac{r+\left| r \right|}{1+r{}^{2}}, \\
& \left. 3 \right)\delta \left( r \right)=\max \left[ 0,\min \left( 2r,1 \right),\min \left( 2,r \right) \right] \\
\end{align}] \f$
<br>\f$r=\frac{u_{i+1+\sigma ,j}^{k}-u_{i+\sigma ,j}^{k}}{u_{i+1,j}^{k}-u_{i,j}^{k}},\f$
<br>\f$ [\sigma =sgn \left[ \left( {{v}_{p}} \right)_{i+0.5,j}^{k} \right].] \f$
<br>\f$\lambda =\frac{\vartriangle t}{\vartriangle x}\f$
<br>
<br>



\param[in] Wt !!!!тут все нужно скопировать с твд схемы!!!!
\param[in] dWdt
\param[in] dWdts
\param[in] pressure
\param[in] mesh
\param[in] index
\param[in] activeElements
\param[in] concentration
\param[in] concentration_temp
\param[in] currentTime

*/

void calculateOpeningSpeedProp(
	std::vector<double> &Wt,
	std::vector<double> &dWdt,
	std::vector<double> &dWdts,
	std::vector<double> &pressure,
	std::vector< std::vector<Cell> > &mesh,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<size_t> > &activeElements,
	std::vector<double> &concentration,				//Вектор концентраций проппанта Света
	std::vector<double> &concentration_temp,		//Вектор новых концентраций проппанта Света
	double currentTime
) {
	dWdt.resize(Wt.size());
	std::fill(dWdt.begin(), dWdt.end(), 0.);

	dWdts.resize(activeElements.size());
	std::fill(dWdts.begin(), dWdts.end(), 0.);

	const double WPow = (2. * bn + 1.) / bn;
	const double PPow = 1. / bn;

	double flowIMinusHalf;		//Егоровские потоки??!!
	double flowIPlusHalf;
	double flowJMinusHalf;
	double flowJPlusHalf;

	////Данные от Светы
	double p_velocityIMinusHalf;						//расчет скорости пропанта в ячейке
	double p_velocityIPlusHalf;
	double p_velocityJMinusHalf;
	double p_velocityJPlusHalf;
	double p_flowIMinusHalf;							//расчет потоков пропанта в ячейке Егор помни, это не твои потоки!!!
	double p_flowIPlusHalf;
	double p_flowJMinusHalf;
	double p_flowJPlusHalf;


	double Rop = 2.65 * 1000;											//плотность проппанта
	double Rof = 1.1 * 1000;												//плотность жидкости
	double g = 9.81;
	double d = 2 * 0.001; 												//диаметр частицы проппанта
	double betta = -2.5*bn;  											//степень для вязкости
	double C_max = 0.6;  												//максимальная концентрация проппанта
	double alfa = 5;  													//степень для корректирующего коэффициента в уравнении для скорости
	double vs = g * (Rop - Rof)*std::pow(d, (bn + 1.)) / (std::pow(3., (bn - 1))*18.);	// скорость оседания пропанта
	double cp = -2.5;				// у Светы это расчетная величина:	cp = betta / bn;

	double injFunc;														//Скорость закачки
	double gradient_concentration;										//Градиент изменения концентрации. Меняется при наличии закачки в центре терщины

		injFunc = 0.5;


	int minI = 10;
	for (size_t k = 0; k < activeElements.size(); ++k) {
		const size_t i = activeElements[k][0];
		const size_t j = activeElements[k][1];

		const double Pij = pressure[index[i][j]];

		flowIMinusHalf = 0.;
		flowIPlusHalf = 0.;
		flowJMinusHalf = 0.;
		flowJPlusHalf = 0.;

		p_flowIMinusHalf = 0.;
		p_flowIPlusHalf = 0.;
		p_flowJMinusHalf = 0.;
		p_flowJPlusHalf = 0.;

		p_velocityIMinusHalf = 0.;
		p_velocityIPlusHalf = 0.;
		p_velocityJMinusHalf = 0.;
		p_velocityJPlusHalf = 0.;

		gradient_concentration = 0.;

		if (LEFT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIMinusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i - 1][j]]);
			const double concentrationIMinusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i - 1][j]]);
			const double pressure_gradient = pressure[index[i - 1][j]] - Pij;

			flowIMinusHalf = -ai::sign(pressure_gradient)
				* std::pow(openingIMinusHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationIMinusHalf / C_max), cp);

			///////////////////////////////////////////////////////
			if ((1. - concentrationIMinusHalf / C_max)<0)
			{
				std::cout << "left " << (1. - concentrationIMinusHalf / C_max) << " ";
			}
			///////////////////////////////////////////////////////////



			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i - 1][j]] >= C_max && pressure_gradient < 0))
				//Нужно уточнить знаки больше или меньше
			{
				p_velocityIMinusHalf = 0;
			}
			else
			{
				p_velocityIMinusHalf = flowIMinusHalf / openingIMinusHalf;
				//расчет новых значений потоков в ячейке в зависимости от направления скорости перетоков между ячейками

				///////////////////////////////////////////////////////
				if ((p_velocityIMinusHalf*pressure_gradient)>0)
				{
					std::cout << " LEFT!!! " << "endl";
				}
				/////////////////////////////////////////////////////////


				if (p_velocityIMinusHalf < 0)
				{
					p_flowIMinusHalf = p_velocityIMinusHalf * concentration[index[i - 1][j]] * Wt[index[i - 1][j]];
				}
				else
				{
					p_flowIMinusHalf = p_velocityIMinusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}
		}



		if (RIGHT != mesh[i][j].xDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingIPlusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i + 1][j]]);
			const double concentrationIPlusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i + 1][j]]);
			const double pressure_gradient = pressure[index[i + 1][j]] - Pij;

			flowIPlusHalf = ai::sign(pressure_gradient)
				* std::pow(openingIPlusHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationIPlusHalf / C_max), cp);

			/////////////////////////////////////////////////////////
			if ((1. - concentrationIPlusHalf / C_max)<0)
			{
				std::cout << "right " << (1. - concentrationIPlusHalf / C_max) << " ";
			}
			///////////////////////////////////////////////////////////

			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i + 1][j]] >= C_max && pressure_gradient < 0))
			{
				p_velocityIPlusHalf = 0;
			}
			else
			{
				p_velocityIPlusHalf = flowIPlusHalf / openingIPlusHalf;

				///////////////////////////////////////////////////////
				if ((p_velocityIPlusHalf*pressure_gradient)<0)
				{
					std::cout << " RIGHT!!! " << "endl";
				}
				/////////////////////////////////////////////////////////




				//расчет новых значений потоков в ячейке в зависимости от направления скорости перетоков между ячейками
				if (p_velocityIPlusHalf > 0)
				{
					p_flowIPlusHalf = p_velocityIPlusHalf * concentration[index[i + 1][j]] * Wt[index[i + 1][j]];
				}
				else
				{
					p_flowIPlusHalf = p_velocityIPlusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}

		}


		//Для нулевого слоя "отражаем значения потоков для симметрии
		if (i00 == i) {
			flowIMinusHalf = -flowIPlusHalf;
			//p_flowIMinusHalf = -p_flowIPlusHalf;
			//p_velocityIMinusHalf = -p_velocityIPlusHalf;
			p_flowIMinusHalf = 0.;
			p_velocityIMinusHalf = 0.;
		}


		if (BOTTOM != mesh[i][j].yDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingJMHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i][j - 1]]);
			const double concentrationJMinusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i][j - 1]]);
			const double pressure_gradient = pressure[index[i][j - 1]] - Pij;

			flowJMinusHalf = -ai::sign(pressure_gradient)
				* std::pow(openingJMHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationJMinusHalf / C_max), cp);

			/////////////////////////////////////////////////////////
			if ((1. - concentrationJMinusHalf / C_max)<0)
			{
				std::cout << "bottom " << (1. - concentrationJMinusHalf / C_max) << " ";
			}
			///////////////////////////////////////////////////////////



			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i][j - 1]] >= C_max && pressure_gradient < 0))
				//Нужно уточнить знаки больше или меньше
			{
				p_velocityJMinusHalf = 0;

			}
			else
			{
				p_velocityJMinusHalf = flowJMinusHalf / openingJMHalf;


				///////////////////////////////////////////////////////
				if ((p_velocityJMinusHalf*pressure_gradient)>0)
				{
					std::cout << " BOTTOM!!! " << "endl";
				}
				/////////////////////////////////////////////////////////


				//Оседание проппанта
				//p_velocityJMinusHalf = p_velocityJMinusHalf - vs * std::pow((1. - concentrationJMinusHalf / C_max), alfa);
				if (p_velocityJMinusHalf < 0)
				{
					p_flowJMinusHalf = p_velocityJMinusHalf * concentration[index[i][j - 1]] * Wt[index[i][j - 1]];
				}
				else
				{
					p_flowJMinusHalf = p_velocityJMinusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}
		}

		if (TOP != mesh[i][j].yDerivative
			&& ALONE != mesh[i][j].xDerivative
			) {
			const double openingJPlusHalf = 0.5
				* (Wt[index[i][j]] + Wt[index[i][j + 1]]);
			const double concentrationJPlusHalf = 0.5
				* (concentration[index[i][j]] + concentration[index[i][j + 1]]);
			const double pressure_gradient = pressure[index[i][j + 1]] - Pij;

			flowJPlusHalf = ai::sign(pressure_gradient)
				* std::pow(openingJPlusHalf, WPow)
				* std::pow(std::abs(pressure_gradient) / dx, PPow)
				/ std::pow((1. - concentrationJPlusHalf / C_max), cp);


			/////////////////////////////////////////////////////////
			if ((1. - concentrationJPlusHalf / C_max)<0)
			{
				std::cout << "top " << (1. - concentrationJPlusHalf / C_max) << " ";
			}
			///////////////////////////////////////////////////////////



			if ((concentration[index[i][j]] >= C_max && pressure_gradient > 0) || (concentration[index[i][j + 1]] >= C_max && pressure_gradient < 0))
			{
				p_velocityJPlusHalf = 0;
			}
			else
			{
				p_velocityJPlusHalf = flowJPlusHalf / openingJPlusHalf;

				///////////////////////////////////////////////////////
				if ((p_velocityJPlusHalf*pressure_gradient)<0)
				{
					std::cout << " TOP!!! " << "endl";
				}
				/////////////////////////////////////////////////////////



				//Оседание проппанта
				//p_velocityJPlusHalf = p_velocityJPlusHalf - vs * std::pow((1. - concentrationJPlusHalf / C_max), alfa);
				if (p_velocityJPlusHalf > 0)
				{
					p_flowJPlusHalf = p_velocityJPlusHalf * concentration[index[i][j + 1]] * Wt[index[i][j + 1]];
				}
				else
				{
					p_flowJPlusHalf = p_velocityJPlusHalf * concentration[index[i][j]] * Wt[index[i][j]];
				}
			}
		}

		const double xDerivative = (flowIPlusHalf - flowIMinusHalf) / dx;
		const double yDerivative = (flowJPlusHalf - flowJMinusHalf) / dy;

		const double value = xDerivative + yDerivative
			+ Q0 * (double)(i == i00 && j == j00) / std::pow(dx, 2)
			- C / std::sqrt(currentTime - mesh[i][j].activationTime);

		dWdt[index[i][j]] = value;
		dWdts[k] = value;



		//Костыль!!!!!!!!!!!!!!!!!!
		//Проверяем непопадание на границу трещины
		if (value > epsilon && Wt[index[i][j]] > epsilon)
			{
			if (i == i00 && j == j00)
			{
				gradient_concentration = ((  ((p_flowIPlusHalf - p_flowIMinusHalf) / dx - (p_flowJMinusHalf - p_flowJPlusHalf) / dy ) + injFunc) *dt_step);
			}
			else
			{
				gradient_concentration = ( ((p_flowIPlusHalf - p_flowIMinusHalf) / dx - (p_flowJMinusHalf - p_flowJPlusHalf) / dy) *dt_step);
			}
			concentration_temp[index[i][j]] = (gradient_concentration
				+ concentration[index[i][j]] * Wt[index[i][j]]) / (Wt[index[i][j]] + value * dt_step);
		}
	}
	concentration = concentration_temp;
	return;
}
