#pragma once
#include <algorithm>
// Cooley-Tukey FFT

/**/
// void fft(std::vector<std::vector<double>> &x, size_t N)
// {
// 	// DFT
// 	size_t k = N, n;
// 	double thetaT = 3.14159265358979 / ((double) N);
//
// 	std::vector<double> phiT = { cos(thetaT), -sin(thetaT) };
// 	while (k > 1)
// 	{
// 		n = k;
// 		k >>= 1;
//
// 		phiT = { phiT[0] * phiT[0] - phiT[1] * phiT[1], 2 * phiT[0] * phiT[1] };
// 		std::vector<double> T { 1. , 0. };
//
// 		for (size_t l = 0; l < k; l++)
// 		{
// 			for (size_t a = l; a < N; a += n)
// 			{
// 				size_t b = a + k;
// 				std::vector<double> t = { x[a][0] - x[b][0], x[a][1] - x[b][1] };
// 				x[a][0] += x[b][0];
// 				x[a][1] += x[b][1];
// 				x[b][0] = t[0] * T[0] - t[1] * T[1];
// 				x[b][1] = t[0] * T[1] + t[1] * T[0];
//
// 			}
// 			double tempT = T[0];
// 			T[0] = T[0]* phiT[0] - T[1] * phiT[1];
// 			T[1] = tempT* phiT[1] + T[1] * phiT[0];
//
// 		}
// 	}
// 	// Decimate
// 	size_t m = (size_t)log2(N);
// 	for (size_t a = 0; a < N; ++a)
// 	{
// 		size_t b = a;
// 		// Reverse bits
// 		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
// 		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
// 		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
// 		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
// 		b = ((b >> 16) | (b << 16)) >> (32 - m);
// 		if (b > a)
// 		{
// 			std::vector<double> t= { x[a][0] , x[a][1] };
// 			x[a][0] = x[b][0];
// 			x[a][1] = x[b][1];
// 			x[b][0] = t[0];
// 			x[b][1] = t[1];
//
// 		}
// 	}
// }

// void fft(std::vector<double> &x, std::vector<double> &y, size_t N)
// {
//     // DFT
//     size_t k = N, n;
//     double thetaT = 3.14159265358979 / ((double)N);
//
//     double phiT0 = cos(thetaT);
//     double phiT1 = -sin(thetaT);
//     while (k > 1)
//     {
//         n = k;
//         k >>= 1;
//
//         double phiT00 = phiT0;
//         phiT0 = phiT0 * phiT0 - phiT1 * phiT1;
//         phiT1 = 2 * phiT00 * phiT1;
//         double T0 = 1.;
//         double T1 = 0.;
//
//         for (size_t l = 0; l < k; l++)
//         {
//             for (size_t a = l; a < N; a += n)
//             {
//                 size_t b = a + k;
//                 double t0 = x[a] - x[b];
//                 double t1 = y[a] - y[b];
//                 x[a] += x[b];
//                 y[a] += y[b];
//                 x[b] = t0 * T0 - t1 * T1;
//                 y[b] = t0 * T1 + t1 * T0;
//             }
//             double tempT = T0;
//             T0 = T0 * phiT0 - T1 * phiT1;
//             T1 = tempT * phiT1 + T1 * phiT0;
//         }
//     }
//
//     // Decimate
//     size_t m = (size_t)log2(N);
//     for (size_t a = 0; a < N; ++a)
//     {
//         size_t b = a;
//         // Reverse bits
//         b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
//         b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
//         b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
//         b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
//         b = ((b >> 16) | (b << 16)) >> (32 - m);
//         if (b > a)
//         {
//             double t0 = x[a];
//             double t1 = y[a];
//             x[a] = x[b];
//             y[a] = y[b];
//             x[b] = t0;
//             y[b] = t1;
//         }
//     }
// }

inline void fft2(std::vector<double> &vector0, std::vector<double> &vector1) {
//inline void fft2(std::vector<double> &Vector0, std::vector<double> &Vector1) {
//    double* vector0 = Vector0.data();
//    double* vector1 = Vector1.data();

    std::size_t k = vector0.size();
    std::size_t j = 0;
    std::size_t n = 0;

    const double length = (double)k;

    double thetaT = 3.14159265358979 / length;
    double swap0 = 0.;
    double swap1 = 0.;
    double T0 = 1.;
    double T1 = 0.;
    double phiT0 = cos(thetaT);
    double phiT1 = -sin(thetaT);

    while (k > 1) {
        n = k;

        k >>= 1;

        swap0 = phiT0;
        swap1 = phiT1;

        phiT0 = swap0 * swap0 - swap1 * swap1;
        phiT1 = 2. * swap0 * swap1;

        T0 = 1.;
        T1 = 0.;

        for (std::size_t l = 0; l < k; ++l)
        {
            for (std::size_t i = l; i < length; i += n)
            {
                j = i + k;

                swap0 = vector0[i] - vector0[j];
                swap1 = vector1[i] - vector1[j];

                vector0[i] += vector0[j];
                vector1[i] += vector1[j];

                vector0[j] = swap0 * T0 - swap1 * T1;
                vector1[j] = swap0 * T1 + swap1 * T0;
            }

            swap0 = T0;

            T0 = swap0 * phiT0 - T1 * phiT1;
            T1 = swap0 * phiT1 + T1 * phiT0;
        }
    }

    std::size_t m = (std::size_t) log2(length);

    for (std::size_t i = 0; i < length; ++i) {
        j = i;

        j = (((j & 0xaaaaaaaa) >> 1) | ((j & 0x55555555) << 1));
        j = (((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2));
        j = (((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4));
        j = (((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8));
        j = ((j >> 16) | (j << 16)) >> (32 - m);

        if (j > i) {
            swap0 = vector0[i];
            swap1 = vector1[i];

            vector0[i] = vector0[j];
            vector1[i] = vector1[j];

            vector0[j] = swap0;
            vector1[j] = swap1;
        }
    }
}

inline void conjugate2(std::vector<double> &vector0, std::vector<double> &vector1) {
    std::size_t length = vector0.size();

    for (std::size_t i = 0; i < length; ++i) {
        vector1[i] = -vector1[i];
    }
}


inline void ifft2(std::vector<double> &vector0, std::vector<double> &vector1) {
    const double length = (double)vector0.size();

    conjugate2(vector0, vector1);

    fft2(vector0, vector1);

    for (std::size_t i = 0; i < length; ++i) {
        vector0[i] /= length;
    }
}
 inline void getWiseElement2(
    std::vector<double> &a0, std::vector<double> &a1,
    std::vector<double> &b0, std::vector<double> &b1,
    std::vector<double> &c0, std::vector<double> &c1){
    std::size_t length = a1.size();

    for(std::size_t i = 0; i < length; ++i){
        double  c0i   = a0[i] * b0[i] - a1[i] * b1[i];
                c1[i] = a0[i] * b1[i] + a1[i] * b0[i];
                c0[i] = c0i;
    }
}
//_______________________________________________________________________________________________
// ������� ������� � ��������� �������� �������������� ����� ��� ������������� �������,
// ������������� ������������ �������� ��������, ������ (-1/2) � (N/2 - 1/2),
//             �.�. ReVector[i] == ReVector[N - 1 - i],     ImVector[i] == 0
//             - �� ����� fft2_RealSymm, � �� ������ ifft2_RealSymm

inline void fft2_RealSymm(std::vector<double> &ReVector, std::vector<double> &ImVector) {

    std::size_t N = ReVector.size();
    std::size_t N_2 = N/2;

    double thetaT = 3.14159265358979 / N;
    double swap0 = 0.;
    //double swap1 = 0.;
    double T0 = 1.;
    double T1 = 0.;
    double phiT0 = cos(thetaT);
    double phiT1 = -sin(thetaT);

    std::vector<double> vector0, vector1;
    vector0.resize(N_2);
    vector1.resize(N_2);

    std::size_t i;
    for (i = 0; i < N_2; ++i) { vector0[i] = ReVector[2*i]; vector1[i] = 0; }

    fft2(vector0, vector1);

    for (i = 0; i < N; ++i) {

        double Re_X_ = T0 * vector0[i%N_2] - T1 * vector1[i%N_2];

        ReVector[i] =  2 * T0 * Re_X_;
        ImVector[i] = -2 * T1 * Re_X_;

        swap0 = T0;

        T0 = swap0 * phiT0 - T1 * phiT1;
        T1 = swap0 * phiT1 + T1 * phiT0;
    }
}

inline void ifft2_RealSymm(std::vector<double> &ReVector, std::vector<double> &ImVector) {

    std::size_t N = ReVector.size();
    std::size_t N_2 = N / 2;

    std::vector<double> vector0, vector1;
    vector0.resize(N_2);
    vector1.resize(N_2);

    std::size_t i;
    for (i = 0; i < N_2; ++i) {
        vector0[i] = (ReVector[i + N_2] + ReVector[i]) / 2;
        vector1[i] = (ImVector[i + N_2] + ImVector[i]) / 2;
    }

    ifft2(vector0, vector1);

    for (i = 0; i < N_2; ++i) {

        ReVector[2 * i] = ReVector[N - 1 - 2 * i] = vector0[i];
        ImVector[2 * i] = ImVector[N - 1 - 2 * i] = 0;
    }
}
