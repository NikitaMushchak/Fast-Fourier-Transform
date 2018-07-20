//#include <iostream>
//#include <complex>
//#include <vector>
//
//#include "ai.hh"
//
//int main(){
//    const size_t N = std::pow(10, 9);
//    
//    std::complex<double> a(1., 1.);
//    std::complex<double> b(3., 2.);
//    std::complex<double> c(0., 0.);
//    
//    std::vector<double> A = {1., 1.};
//    std::vector<double> B = {3., 2.};
//    std::vector<double> C = {0., 0.};
//    
//    auto start1 = ai::time();
//    
//    for(size_t i = 0; i < N; ++i){
//        c = a * b;
//    }
//    
//    auto finish1 = ai::time();
//    
//    std::cout << c << std::endl;
//    
//    std::cout << "Time complex: " << ai::duration(start1, finish1, "ms") 
//        << "ms" << std::endl;
//    
//    auto start2 = ai::time();
//    
//    for(size_t i = 0; i < N; ++i){
//        C[i][0] = A[0] * B[0] - A[1] * B[1];
//        C[i][1] = A[0] * B[1] + A[1] * B[0];
//    }
//    
//    auto finish2 = ai::time();
//    
//    std::cout << "Time vector: " << ai::duration(start2, finish2, "ms") 
//        << "ms" << std::endl;
//    
//    std::cout << "(" << C[0] << "," << C[1] << ")" << std::endl;
//    
//    return 0;
//}