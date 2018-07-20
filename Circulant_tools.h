#pragma once
#include "ai.hh"

//generate circulant
void createCirculantMatrix(
	std::vector< std::vector<double> > &matrix,
	std::vector<double> &vector
) {
	size_t length = vector.size();

	size_t displacement = 0;

	matrix.resize(length);

	for (size_t i = 0; i < length; ++i) {
		matrix[i].resize(length);

		size_t k = length - 1;
		for (size_t j = displacement; j < length; ++j) {
			matrix[i][k] = vector[j];

			--k;
		}
		for (size_t j = 0; j < displacement; ++j) {
			matrix[i][k] = vector[j];

			--k;
		}
		++displacement;
	}
}
///wise element miltiplication
void getWiseElement(CArray &a, CArray &b, CArray &w, size_t m)
{
	for (size_t i = 0; i < m; ++i)
	{
		//w[i] = a[i] * b[i];
	}
}
void conjugate(CArray &x, CArray & x1, size_t N)
{
	for (size_t i = 0; i < N; i++)
	{
		x1[i][1] = -x[i][1];
	}
	//for (auto& x1 : x)
	//	x1 = std::conj(x1);
}
// Cooley-Tukey FFT 
void multiplication(CArray& x, CArray y, )
{

}
void fft(CArray &x, CArray &y, size_t N)
{
	// DFT
	//unsigned int N = x.size(), k = N, n;
	size_t k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT[0] = phiT[0] * phiT[0];
		phiT[1] = phiT[1] * phiT[1];

		Complex T = Complex(1.0L , 0);

		for (size_t l = 0; l < k; l++)
		{
			for (size_t a = l; a < N; a += n)
			{
				size_t b = a + k;
				Complex t = Complex (x[a][0] - x[b][0], x[a][1] - x[b][1]);
				x[a][0] += x[b][0];
				x[a][1] += x[b][0];
				x[b] = t * T; //write vector/vector multiplication

			}
			T *= phiT;
		}
	}
	// Decimate
	size_t m = (size_t)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		size_t b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
			//y[b] = x[b];
		}

		y = x;

	}
}
// inverse fft 
/**/
void ifft(CArray& x, CArray& v, size_t N)
{

	CArray y(N);
	// conjugate the complex numbers
	conjugate(x, v, N);
	// fft
	fft(v, y, N);
	// conjugate the complex numbers again
	conjugate(y, v, N);
	// scale the numbers
	for (size_t i = 0; i < N; i++)
	{
		v[i] /= N;

	}

}
