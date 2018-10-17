#include <vector>
#include <iostream>

#include "ai.hh"
#include "planar3D.hh"

/// \brief Main-функция
/// \details Первичная функция, считывающая параметры, переданные при запуске
/// программы, и запускающая расчёт распространения трещины по планарной модели
int main(const int argc, const char *argv[]){
    std::string compiler = "%compiler%";
    std::string target = "%target%";
    std::string timestamp = "%timestamp%";
    buildVersion = "2.1.1+%platform%";

    bool saveSteps = false;
    bool runningFromGUI = false;

    std::string pathToBarriersFile = std::string();
    std::string pathToInjectionFile = std::string();

    double cellSize = 5.;//6.;
    double meshSize = 6.2;//5.;
    double mu = 0.4;
    double T1 = 12.0;//25.0;
    double ts = 60.;

    E = 36 * pow(10, 9);
    epsilon = 0.0000001;//std::pow(10., -8);
    C = 50. * std::pow(10., -6);//0. * std::pow(10., -6);
    Kic = 0.;//0. * std::pow(10., 6);
    Q0 = 0.08;//2.2 / 60.;//4.4/60.
    bn = 1.;

    ///TODO: specify regime
    regime = 1;

    for(int i = 1; i < argc; ++i){
        if("-v" == std::string(argv[i]) || "--version" == std::string(argv[i])){
            std::cout << "  Build: "  << buildVersion<< "." << std::endl;
            std::cout << "  Compiler: "  << compiler<< "." << std::endl;
            std::cout << "  Target: "  << target<< "." << std::endl;
            std::cout << "  AiLibrary: " << ai::getVersion() << "."
                << std::endl;
            std::cout << "  Compilation timestamp: "  << timestamp << "."
                << std::endl;

            return 0;
        }

        if("-h" == std::string(argv[i]) || "--help" == std::string(argv[i])){
            std::cout << "usage: planar3D [options]"
                << std::endl
                << "    -h  --help            print this usage and exit"
                << std::endl
                << "    -v  --version         print build info and exit"
                << std::endl
                << "    --list-errors         print possible errors ans exit"
                << std::endl << std::endl

                << "  Liquid parameters" << std::endl
                << "    --Q=<value>           injection [double, m^3 / s]"
                << std::endl
                << "    --mu=<value>          viscosity [double, Pa * s]"
                << std::endl
                << "    --n=<value>           index of the liquid's rheology "
                << "[double, n/d]"
                << std::endl << std::endl

                << "  Shelf parameters" << std::endl
                << "    --E=<value>           flat Young's modulus "
                << "[double, GPa]"
                << std::endl
                << "    --C=<value>           Carter coefficient [double, "
                << "um / s^0.5]"
                << std::endl
                << "    --Kic=<value>         stress intensity factor [double, "
                << "MPa/m^0.5]"
                << std::endl << std::endl

                << "  Time parameters" << std::endl
                << "    --time=<value>        modeling time [double, min]"
                << std::endl
                << "    --ts=<value>          time scale [double, s]"
                << std::endl << std::endl

                << "  Mesh parameters" << std::endl
                << "    --mesh-size=<value>    how many initial cracks fit the "
                << "mesh [uint, n/d]"
                << std::endl
                << "    --cell-size=<value>    how many cells fit the initial "
                << "crack [uint, n/d]"
                << std::endl << std::endl

                << "  Fracture regime" << std::endl
                << "    --viscosity           viscosity dominated regime "
                << "{defualt}"
                << std::endl
                << "    --toughness           toughness dominated regime "
                << std::endl
                << "    --leak-off            leak-off dominated regime "
                << std::endl << std::endl

                << "  Compressive stress barriers" << std::endl
                << "    --barriers=<path>     path to a txt-file with a list "
                << "of barriers [string]"
                << std::endl
                << "    --barriers            default load from "
                << "./InitialConditions/barriers.txt"
                << std::endl << std::endl

                << "  Dynamic Injection" << std::endl
                << "    --injection=<path>     path to a txt-file with a list "
                << "of injections [string]"
                << std::endl
                << "    --injection            default load from "
                << "./InitialConditions/injection.txt"
                << std::endl << std::endl

                << "  Other flags" << std::endl
                << "    --save-steps           save pressure and opening at "
                << "every check"
                << std::endl
                << "    --env=gui             flag for planar3DGUI"
                << std::endl;

            return 0;
        }

        if("--list-errors" == std::string(argv[i])){
            std::cout << "  User input errors" << std::endl
                << "    Code 11. Cell size is not a positive integer."
                << std::endl
                << "    Code 12. Cell size is less than 5."
                << std::endl
                << "    Code 13. Mesh size is not a positive integer."
                << std::endl
                << "    Code 14. Carter coefficient isn't positive in a "
                << "leak-off dominated regime."
                << std::endl
                << "    Code 15. Stress intensity factor isn't positive in a "
                << "toughness dominated "
                << std::endl
                << "             regime."
                << std::endl
                << "    Code 16. Viscosity isn't positive in a "
                << "viscosity dominated regime."
                << std::endl << std::endl

                << "  Enviroment errors" << std::endl
                << "    Code 21. Cannot find %folderName%."
                << std::endl
                << "    Code 22. Cannot open %fileName%."
                << std::endl;

            return 0;
        }

        if(
            ai::assignAbsDoubleParameter(argv[i], "--Q=", Q0)
            || ai::assignAbsDoubleParameter(argv[i], "--n=", bn)
            || ai::assignAbsDoubleParameter(argv[i], "--mu=", mu)
            || ai::assignAbsDoubleParameter(argv[i], "--time=", T1)
            || ai::assignAbsDoubleParameter(argv[i], "--ts=", ts)
            || ai::assignAbsDoubleParameter(argv[i], "--cell-size=", cellSize)
            || ai::assignAbsDoubleParameter(argv[i], "--mesh-size=", meshSize)
            || ai::assignStringParameter(argv[i],
                "--barriers=", pathToBarriersFile
            )
            || ai::assignStringParameter(argv[i],
                "--injection=", pathToInjectionFile
            )
        ){
            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--E=", E)){
            E *= std::pow(10, 9);
            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--C=", C)){
            C *= std::pow(10., -6);
            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--Kic=", Kic)){
            Kic *= std::pow(10., 6);
            continue;
        }

        if("--viscosity" == std::string(argv[i])){
            regime = 1;
            continue;
        }

        if("--toughness" == std::string(argv[i])){
            regime = 0;
            continue;
        }

        if("--leak-off" == std::string(argv[i])){
            regime = -1;
            continue;
        }

        if("--barriers" == std::string(argv[i])){
            pathToBarriersFile = "./InitialConditions/barriers.txt";
            continue;
        }

        if("--injection" == std::string(argv[i])){
            pathToInjectionFile = "./InitialConditions/injection.txt";
            continue;
        }

        if("--save-steps" == std::string(argv[i])){
            saveSteps = true;
            continue;
        }

        if("--env=gui" == std::string(argv[i])){
            runningFromGUI = true;
            continue;
        }
    }

    // if(cellSize != (int) cellSize){
    //     std::cerr << "Cell size should be a positive integer." << std::endl;
    //     return 11;
    // }
    //
    // if(5. > cellSize){
    //     std::cerr << "Cell size should be at least 5." << std::endl;
    //     return 12;
    // }
    //
    // if(meshSize != (int) meshSize){
    //     std::cerr << "Mesh size should be a positive integer." << std::endl;
    //     return 13;
    // }

    if(-1 == regime && epsilon > C){
        std::cerr << "Carter coefficient should be positive in a leak-off "
            << "dominated regime."
            << std::endl;
        return 14;
    }

    if(0 == regime && epsilon > Kic){
        std::cerr << "Stress intensity factor should be positive in a "
            << "toughness dominated regime."
            << std::endl;
        return 15;
    }

    if(1 == regime && epsilon > mu){
        std::cerr << "Viscosity should be positive in a viscosity dominated"
            << "regime."
            << std::endl;
        return 16;
    }

    const double theta = 2. * std::pow(2. * (2. * bn + 1.) / bn, 1. / bn);
    alpha = 2. / (bn + 2.);
    const double BAlpha = 0.25 * alpha * tan(0.5 * M_PI - M_PI * (1 - alpha));
    Amu = std::pow(BAlpha * (1 - alpha), -0.5 * alpha);

    return planar3D(
        T1, mu, theta, ts, cellSize, meshSize,
        pathToBarriersFile, pathToInjectionFile, runningFromGUI, saveSteps
    );
}
