#pragma once
#include "ai.hh"

///wise element miltiplication vector-vector
inline void getWiseElement(std::vector< std::vector<double> > &a, std::vector< std::vector<double> > &b, std::vector< std::vector<double> > &c , size_t m)
{
	for (size_t i = 0; i < m; ++i)
	{
		c[i][0] = a[i][0] * b[i][0] - a[i][1] * b[i][1];
		c[i][1] = a[i][0] * b[i][1] + a[i][1] * b[i][0];
	}
}

//cojugate vector with std::vector<double> data
inline void conjugate(std::vector< std::vector<double> > &x, size_t N)
{
	size_t i = 0;

	for (; i <= N - 4; i += 4)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];

		x[i+1][0] =  x[i+1][0];
		x[i+1][1] = -x[i+1][1];

		x[i + 2][0] = x[i + 2][0];
		x[i + 2][1] = -x[i + 2][1];

		x[i + 3][0] = x[i + 3][0];
		x[i + 3][1] = -x[i + 3][1];
	}
	for (; i < N; ++i)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];
	}
}
// Cooley-Tukey FFT

/**/
inline void fft(std::vector< std::vector<double> > &x, size_t N)
{
	// DFT
	size_t k = N, n;
	double thetaT = 3.14159265358979 / ((double) N);

	std::vector<double> phiT = { cos(thetaT), -sin(thetaT) };
	std::vector<double> T;
	while (k > 1)
	{
		n = k;
		k >>= 1;

		phiT = { phiT[0] * phiT[0] - phiT[1] * phiT[1], 2 * phiT[0] * phiT[1] };
		std::vector<double> T { 1. , 0. };

		for (size_t l = 0; l < k; l++)
		{
			for (size_t a = l; a < N; a += n)
			{
				size_t b = a + k;
				std::vector<double> t = { x[a][0] - x[b][0], x[a][1] - x[b][1] };
				x[a][0] += x[b][0];
				x[a][1] += x[b][1];
				x[b][0] = t[0] * T[0] - t[1] * T[1];
				x[b][1] = t[0] * T[1] + t[1] * T[0];

			}
			double tempT = T[0];
			T[0] = T[0]* phiT[0] - T[1] * phiT[1];
			T[1] = tempT* phiT[1] + T[1] * phiT[0];

		}
	}
	// Decimate
	size_t m = (size_t)log2(N);
	for (size_t a = 0; a < N; ++a)
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
			std::vector<double> t { x[a][0] , x[a][1] };
			x[a][0] = x[b][0];
			x[a][1] = x[b][1];
			x[b][0] = t[0];
			x[b][1] = t[1];

		}
	}
}
// inverse fft
/**/
/**/
inline void ifft(std::vector< std::vector<double> >& x, size_t N)
{
	// conjugate the std::vector<double> numbers
	conjugate(x, N);
	// fft
	fft(x, N);
	// scale the numbers
	for (size_t i = 0; i < N; i++)
	{
		x[i][0] /= N;
		//x[i][1].erase();

	}
}
/**/
