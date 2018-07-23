
#include <iostream>
#include <random>
#include <vector>
#include "ai.hh"
#include "Circulant_tools.h"
#include "influence_vector.h"

int main()
{
	size_t N = 8192;//length of vector 2^n!!!!!!!!!!!!
    size_t counter = 10;
	std::cout << "Dimention of matrix : " << N << std::endl;;
	if ((size_t)log2(N) % 2 != 0)
	{
		//std::cout << "Dimention of the matrix will be 2^n!!!!" ;

	}
	/// Complex circulant
	//CMatrix A;//  now not actual 20.07
	CArray y(N); // vector with Complex data
	CArray x(N); // first coloumn of circulant complex data
				 //CArray x = {1.0 , 12.0, 11.0, 10.0};
				 //CArray x = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
				 // для матрицы
	std::vector<  double> vectDirect2; // vector
	std::vector< double> vectDirect1; // generated vector 
	std::vector< std::vector<double> > loaded; // created circulant
	for (size_t i = 0; i < y.size(); i++)
	{
		y[i].resize(2);
	}
	//std::cout<<"y resized"<<" "<<std::endl;
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i].resize(2);
	}
	//std::cout << "x resized" << " " << std::endl;
	// generate vector
	generateVector(vectDirect1, N);
	// generate circulant
	createCirculantMatrix(loaded, vectDirect1);

	//ai::parseFileInMatrix("C:\\Temp\\circ16.txt", '\t', loaded);
	/*
	std::cout<<"Circulant martrix"<<std::endl;
	ai::printMatrix(loaded);
	/**/
	// circulant with Complex data
	/*A.resize(loaded.size());
	for (size_t i = 0; i < loaded.size(); i++)
	{
	A[i].resize(loaded[i].size());
	for (size_t j = 0; j < loaded[i].size(); j++)
	{
	A[i][j] = loaded[i][j];
	}
	}*/

	// vector with complex data
	//std::cout << "vector of expantion " << std::endl;
	y.resize(vectDirect1.size());
	for (size_t i = 0; i < vectDirect1.size(); i++)
	{
		y[i][0] = vectDirect1[i];
		y[i][1] = 0.;
		//std::cout << "(" << y[i][0] << "," << y[i][1] << ")" << std::endl;

	}

	//std::cout << " first coloumn of circulant " << std::endl;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
			if (j < 1)
			{
				x[i][0] = loaded[i][j];// для матрицы
				x[i][1] = 0.;
				//std::cout << "(" << x[i][0] << "," << x[i][1] << ")" << std::endl;
				//std:: cout << std::endl;

			}
	}

	// redefine first row vector as Complex


	// CArray y = { 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0,1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0 };
	//generate vector
	///std::vector<double> vectEgor1 = { 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0,1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 7.0, 8.0 };
	//	 std::vector<double> vectDirect2;
	auto timeEgorStart = ai::time();
	for (int i = 0; i < counter; ++i)
	{
	multiply(loaded, vectDirect1, vectDirect2);
	}
	auto timeEgorFinish = ai::time();
	std::cout << "Direct multiplication time:  " << ai::duration(timeEgorStart, timeEgorFinish, "ms") << "ms" << std::endl;
	/*std::cout << "Direct product" << std::endl;
	for (size_t i = 0; i < N; ++i)
	{
		std::cout << vectDirect2[i] << std::endl;
	}*/

	CArray h(N);
	for (size_t i = 0; i < h.size(); i++)
	{
		h[i].resize(2);
	}
	/*
	1. вычисляем для вектора f = FFT(y),
	2. вычисляем для первого вектора циркулянта g = FFT(x),
	3. вычиляем the element wise vector - vector product h = f ∗ g,
	4. вычиляем z = IFFT(h)
	/**/
	
	// быстрое преобразование Фурье для вектора и циркулянта
	/**/
	auto start = ai::time();
	for (int i = 0; i < counter; ++i)
	{
		//от вектора
		fft(y, N);
		//ifft(f, fou, N);


	//auto finish = ai::time();

	// от циркулянта
		fft(x, N);
		getWiseElement(x, y, h, N);

		ifft(h, N);

	}
	auto finish = ai::time();
	std::cout << "Fourier: " << ai::duration(start, finish, "ms") << "ms" << std::endl;
	std::cout << "Coeff FFT/Direct: " << ai::duration(start, finish, "ms") / ai::duration(timeEgorStart, timeEgorFinish, "ms") << std::endl;
	/**/
	
	//auto finish = ai::time();
	//std::cout << "FFT vector" << std::endl;
	//for (size_t i = 0; i < N; ++i)
	//{
	//	std::cout << "(" << f[i][0] << "," << f[i][1] << ")" << std::endl;
	//}
	
	//std::cout << "product of circulant" << std::endl;
	//for (size_t i = 0; i < N; ++i)
	//{
	//	std::cout << "(" << g[i][0] << "," << g[i][1] << ")" << std::endl;
	//}

	//std::cout << "wise element" << std::endl;
	//for (size_t i = 0; i < N; ++i)
	//{
	//	std::cout << "(" << h[i][0] << "," << h[i][1] << ")" << std::endl;
	//}
	///**/
	/*std::cout << "product Fourier" << std::endl;
	for (size_t i = 0; i < N; ++i)
	{
		std::cout << "(" << z[i][0] << "," << z[i][1] << ")" << std::endl;
	}
	std::cout << "Fourier: " << ai::duration(start, finish, "us") << "us" << std::endl;*/

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
