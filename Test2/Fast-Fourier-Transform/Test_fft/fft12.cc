#include <iostream>
#include <complex>
#include <string>
#include <vector>

//#include <ai>

#include "Circulant_tools.hh"
#include "Circulant_tools2.hh"
#include "ai.hh"

void compareFFT12(
    const std::size_t power = 8,
    const std::string timeUnit = "ms"
){
    const std::size_t N = std::pow(10, power);

    std::cout << "Running " << N << " (10^" << power << ") test(s)"
        << std::endl;


    std::cout << std::endl;

    const std::size_t size = std::pow(2, 12);

    std::vector<double> source1;
    ai::generateRandomVector(source1, size, -10., 10.);
    std::vector<double> source2;
    ai::generateRandomVector(source2, size, -10., 10.);

    std::vector< std::vector<double> > a;
    for(std::size_t i = 0; i < size; ++i){
        std::vector<double> f;
        f[0] = source1[i];
        f[2] = 0.;
        a.push_back(f);
    }
    std::vector< std::vector<double> > b;
    for(std::size_t i = 0; i < size; ++i){
        std::vector <double> h ;
        h[0]= source2[i];
        h[1]= 0.;
        b.push_back(h);
    }
    std::vector< std::vector<double> > c;
    for(std::size_t i = 0; i < size; ++i){
        std::vector<double>  n ;
        n[0]=0.;
        n[1]= 0.;
        c.push_back(n);
    }

    std::vector< std::vector<double> > A = a;
    std::vector< std::vector<double> > B = b;
    std::vector< std::vector<double> > C = c;

    fft(a, a.size());
    fft(b, b.size());

    auto startFFT1 = ai::time();

    for(std::size_t i = 0; i < N; ++i){
        getWiseElement(a, b, c, a.size());
        ifft(c, c.size());
    }

    auto finishFFT1 = ai::time();

    double fft1Time = ai::duration(startFFT1, finishFFT1, timeUnit);

    std::cout << "Time for FFTv1: " << fft1Time << timeUnit
        << std::endl;

    std::cout << "First final value is " << c[0][0] << std::endl;

    std::cout << std::endl;

    fft2(A);
    fft2(B);

    auto startFFT2 = ai::time();

    for(std::size_t i = 0; i < N; ++i){
        getWiseElement2(A, B, C);
        ifft2(C);
    }

    auto finishFFT2 = ai::time();

    double fft2Time = ai::duration(startFFT2, finishFFT2, timeUnit);

    std::cout << "Time for FFTv2: " << fft2Time << timeUnit
        << std::endl;

    std::cout << "First final value is " << C[0][0] << std::endl;
}

int main(const int argc, const char *argv[]){
    std::size_t power = 3;
    std::string timeUnit = "ms";

    for(int i = 1; i < argc; ++i){
        if("-h" == std::string(argv[i]) || "--help" == std::string(argv[i])){
            std::cout << "usage: fft12 [options]"
                << std::endl
                << "    -h  --help            print this usage and exit"
                << std::endl << std::endl
                << "    --power=<value>       specify number of tests in tens "
                << "[integer]"
                << std::endl
                << "    --time-unit=<value>   specify time unit [s, ms, us]"
                << "[string]";

            return 0;
        }

        if(
            ai::assignParameter(argv[i], "--power=", power)
            || ai::assignByCheckingParameter(
                argv[i], "--time-unit=s", timeUnit, std::string("s")
            )
            || ai::assignByCheckingParameter(
                argv[i], "--time-unit=ms", timeUnit, std::string("ms")
            )
            || ai::assignByCheckingParameter(
                argv[i], "--time-unit=us", timeUnit, std::string("us")
            )
        ){
            continue;
        }
    }

    compareFFT12(power, timeUnit);

    return 0;
}
