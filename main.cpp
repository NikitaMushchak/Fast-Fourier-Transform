
#include <iostream>
#include <random>
#include <vector>
#include "ai.hh"
#include "Circulant_tools.h"
#include "influence_vector.h"

int main()
{
	size_t N = 128;//length of vector 2^n!!!!!!!!!!!!
	std::cout << "Dimention of matrix : " << N << std::endl;;
	if ((size_t)log2(N) % 2 != 0)
	{
		//std::cout << "Dimention of the matrix will be 2^n!!!!" ;
		
	}
	/// Complex circulant
	CMatrix A; 
	CArray y; // vector with Complex data
	CArray x(N); // first coloumn of circulant
	//CArray x = {1.0 , 12.0, 11.0, 10.0};
	//CArray x = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
	// для матрицы
	std::vector<double> vectDirect2; // vector
	std::vector<double> vectDirect1; // generated vector 
	std::vector< std::vector<double> > loaded; // created circulant
	/*for (size_t i = 0; i < C.size(); i++)
	{
		C[i].resize(N);
	}*/
	// generate vector
	generateVector(vectDirect1, N);
	// generate circulant
	createCirculantMatrix(loaded, vectDirect1);
	
	//ai::parseFileInMatrix("C:\\Temp\\circ16.txt", '\t', loaded);
	/*std::cout<<"Martrix!!!"<<std::endl;
	ai::printMatrix(loaded);*/
	
	// circulant with Complex data
	A.resize(loaded.size());
	for (size_t i = 0; i < loaded.size(); i++)
	{
		A[i].resize(loaded[i].size());
		for (size_t j = 0; j < loaded[i].size(); j++)
		{
			A[i][j] = loaded[i][j];
		}
	}
	// vector with complex data
	y.resize(vectDirect1.size());
	for (size_t i = 0; i < vectDirect1.size(); i++)
	{
		y[i] = vectDirect1[i];
		
	}

	//std::cout << " first coloumn of circulant ";
	for (size_t i = 0; i<N; i++)
	{
		for (size_t j = 0; j < N; j++)
			if (j < 1)  x[i]=A[i][j];// для матрицы
			//std::cout << x[i] << " ";
		//std:: cout << std::endl;
	}
    // CArray y = { 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0,1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0 };
	 //generate vector
	 ///std::vector<double> vectEgor1 = { 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0,1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0 };
//	 std::vector<double> vectDirect2;
	 auto timeEgorStart = ai::time();
	 for (int i = 0; i < 1000; ++i)
	 {
		 multiply(loaded, vectDirect1, vectDirect2);
	 }
	 auto timeEgorFinish = ai::time();
	 std::cout << "Direct multiplication time:  "<<ai::duration(timeEgorStart, timeEgorFinish, "us") << "us" << std::endl;
	/* std::cout << "Direct product" << std::endl;
	 for (size_t i = 0; i < N; ++i)
	 {
		 std::cout << vectDirect2[i] << std::endl;
	 }*/

	CArray g(N);
	CArray f(N);
	CArray h(N);
	CArray z(N);
	
	/*
	    1. вычисляем для вектора f = FFT(y),
		2. вычисляем для первого вектора циркулянта g = FFT(x),
		3. вычиляем the element wise vector - vector product h = f∗ g,
		4. вычиляем z = IFFT(h) 
			/**/
	
	// быстрое преобразование Фурье для вектора и циркулянта
	
	auto start = ai::time();
	for (int i = 0; i < 1000; ++i)
	{
		//от вектора
		fft(y, f, N);
		// от циркулянта
		fft(x, g, N);
		getWiseElement(f, g, h, N);

		ifft(h, z, N);
	}
	auto finish = ai::time();
	//auto finish = ai::time();
	
	/*std::cout << "product" << std::endl;
	for (size_t i = 0; i < N; ++i)
	{
		std::cout << z[i] << std::endl;
	}*/
	
	std::cout << "Fourier: " << ai::duration(start, finish, "us") << "us" << std::endl;

	std::cout << "Direct to Fourier ratio"<< ai::duration(timeEgorStart, timeEgorFinish) /
		ai::duration(start, finish, "us") << "us" << std::endl;
	/*
	std::cout << "FFT of vector" << std::endl;
	for (size_t i = 0; i < N; ++i)
	{
		std::cout << prody[i] << std::endl;
		//std::cout << " sup"<<f << std::endl;
	}
	/**/
	
	/*
	std::cout << "wise element" << std::endl;
	for (size_t i = 0; i < N; ++i)
	{
		std::cout << prodw[i] << std::endl;
	}
	std::cout << std::endl << "IFFT" << std::endl;
	/**/
	
	//_getch();

	return 0;
}