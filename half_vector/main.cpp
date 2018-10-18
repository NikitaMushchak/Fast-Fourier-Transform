//﻿#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <random>
#include <vector>
#include <time.h>
#include "Circulant_tools.h"
#include "Circulant_tools2.hh"
#include "influence_vector.h"
#include "Circulant_tools3.hh"

int main()
{
	size_t N = 4096;//length of vector 2^n!!!!!!!!!!!!
    size_t counter = 1;
	std::cout << "Dimension of matrix : " << N << std::endl;;
	//if ((size_t)log2(N) % 2 != 0)
	if (std::pow(2, (size_t)log2(N)) != N)
	{
		std::cout << "Dimension of the matrix will be 2^n!!!!" ;

	}

	// generate vector
    std::vector<double> vectDirect0, vectDirect1, vectDirect2, vectDirect3, influence, influence1, influencef;
    vectDirect1.resize(N);
    vectDirect2.resize(N);
    vectDirect3.resize(N);
	influence.resize(N);
	influencef.resize(N);
	for (size_t i = 0 ; i < N;++i)
	{
		influencef[i].resize(2);
	}
	// vector with complex data
	//std::cout << "vector of expantion " << std::endl;

    // быстрое преобразование Фурье для вектора и циркулянта

    // long t = clock();
    double Eps = 0, EpsRe = 0, EpsIm = 0;

    // for (int i = 0; i < counter; ++i)
    // {
        generateVector(influence, N);
        generateSymm(vectDirect0, N); //symmetric opening
        vectDirect2 = vectDirect0; // opening vector copy
		// std::cout<<"influence"<<std::endl;
		// ai::printVector(influence);

			double Norm = 0;

        for (size_t i = 0; i < N; i++)
        {
            vectDirect1[i] = vectDirect3[i] = 0;
			influence1.push_back(0.);
			influencef[i][0]=influence[i]);
            double n1 = vectDirect0[i];
            Norm += n1 * n1;
        }
		 // ai::printVector(vectDirect1);
        // Norm /= N;
        // Norm = sqrt(Norm);


        fft2(vectDirect2, vectDirect3);


		fft2_RealSymm(influence, influence1);
        // double eps_r = 0, eps_i = 0;
        // for (size_t i = 0; i < N; i++)
        // {
        //     double e1 = vectDirect0[i] - vectDirect2[i];
        //     eps_r += e1 * e1;
        //     //printf("%f\n", e1);
        //     //system("pause");
		//
        //     double e2 = vectDirect1[i] - vectDirect3[i];
        //     eps_i += e2 * e2;
        // }
        // eps_i /= N;
        // eps_i = sqrt(eps_i);
        // EpsIm += eps_i;
        // eps_r /= N;
        // eps_r = sqrt(eps_r);
        // EpsRe += eps_r;

		auto time1 = ai::time();
		fft2_RealSymm(vectDirect0, vectDirect1);
		std::vector<double> product1, product2;
		 product1.resize(N);
		 product2.resize(N);
		getWiseElement2(vectDirect0, vectDirect1,influence, influence1,product1, product2);
						//real       imagine

		ifft2_RealSymm(product1, product2);
		auto time2 = ai::time();

        // ifft2(vectDirect2, vectDirect3);

        // double eps = 0;
        // for (size_t i = 0; i < N; i++)
        // {
        //     double e1 = vectDirect0[i] - vectDirect2[i];
        //     eps += e1 * e1;
        // }
        // eps /= N;
        // eps = sqrt(eps) / Norm;
        // Eps += eps;
    // }

    // Eps /= counter;
    // EpsIm /= counter;
    // EpsRe /= counter;
    // t = clock() - t;
	//
    // printf("time= %.3f sec\n", t / 1000.);
    // printf("Eps= %e\n", Eps);
    // printf("EpsIm= %e\n", EpsIm);
    // printf("EpsRe= %e\n", EpsRe);
    // system("pause");

	return 0;
}
