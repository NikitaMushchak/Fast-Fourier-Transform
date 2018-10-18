// #pragma once
// #include <algorithm>
//
// #include "ai.hh"
//
// inline void getWiseElement2(
//     std::vector< std::vector<double> > &a,
//     std::vector< std::vector<double> > &b,
//     std::vector< std::vector<double> > &c
// ){
//     std::size_t length = a.size();
//
//     for(std::size_t i = 0; i < length; ++i){
//         c[i][0] = a[i][0] * b[i][0] - a[i][1] * b[i][1];
//         c[i][1] = a[i][0] * b[i][1] + a[i][1] * b[i][0];
//     }
// }
//
// inline void conjugate2(std::vector< std::vector<double> > &vector){
//     std::size_t length = vector.size();
//
//     for(std::size_t i = 0; i < length; ++i){
//         vector[i][1] = -vector[i][1];
//     }
// }
//
// inline void fft2(std::vector< std::vector<double> > &vector){
//     std::size_t k = vector.size();
//     std::size_t j = 0;
//     std::size_t n = 0;
//
//     const double length = (double) k;
//
//     double thetaT = 3.14159265358979 / length;
//     double swap0 = 0.;
//     double swap1 = 0.;
//     double T0 = 1.;
//     double T1 = 0.;
//     double phiT0 = cos(thetaT);
//     double phiT1 = -sin(thetaT);
//
//     while(k > 1){
//         n = k;
//
//         k >>= 1;
//
//         swap0 = phiT0;
//         swap1 = phiT1;
//
//         phiT0 = swap0 * swap0 - swap1 * swap1;
//         phiT1 = 2. * swap0 * swap1;
//
//         T0 = 1.;
//         T1 = 0.;
//
//         for(std::size_t l = 0; l < k; ++l)
//         {
//             for(std::size_t i = l; i < length; i += n)
//             {
//                 j = i + k;
//
//                 swap0 = vector[i][0] - vector[j][0];
//                 swap1 = vector[i][1] - vector[j][1];
//
//                 vector[i][0] += vector[j][0];
//                 vector[i][1] += vector[j][1];
//
//                 vector[j][0] = swap0 * T0 - swap1 * T1;
//                 vector[j][1] = swap0 * T1 + swap1 * T0;
//             }
//
//             swap0 = T0;
//
//             T0 = swap0 * phiT0 - T1 * phiT1;
//             T1 = swap0 * phiT1 + T1 * phiT0;
//         }
//     }
//
//     std::size_t m = (std::size_t) log2(length);
//
//     for(std::size_t i = 0; i < length; ++i){
//         j = i;
//
//         j = (((j & 0xaaaaaaaa) >> 1) | ((j & 0x55555555) << 1));
//         j = (((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2));
//         j = (((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4));
//         j = (((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8));
//         j = ((j >> 16) | (j << 16)) >> (32 - m);
//
//         if(j > i){
//             swap0 = vector[i][0];
//             swap1 = vector[i][1];
//
//             vector[i][0] = vector[j][0];
//             vector[i][1] = vector[j][1];
//
//             vector[j][0] = swap0;
//             vector[j][1] = swap1;
//         }
//     }
// }
//
// inline void ifft2(std::vector< std::vector<double> > &vector){
//     const double length = (double) vector.size();
//
//     conjugate2(vector);
//
//     fft2(vector);
//
//     for(std::size_t i = 0; i < length; ++i){
//         vector[i][0] /= length;
//     }
// }
