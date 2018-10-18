#include <vector>
#include <iostream>

#include "ai.hh"
#include "planar3D.hh"

/// \details Функция рассчитывает начальное раскрытие в точке по 
/// автомодельному решению для раскрытия и радиуса трещины
double getInitialOpening(
    const double x,
    const double y,
    const double radius,
    const std::vector<double> zP,
    const std::vector<double> W0
){
    const double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2)) / radius;
    
    double opening = 0;
    
    if(1 > distance){
        for(size_t i = 0; i < zP.size() - 1; ++i){
            if(distance >= zP[i] && distance < zP[i + 1]){
                opening = W0[i] * zP[i + 1] - W0[i + 1] * zP[i]
                    - distance * (W0[i] - W0[i + 1]);
                opening /= (zP[i + 1] - zP[i]);
            }
        }
    }
    
    return opening;
}

/// \details Функция подгружает данные о контрасте напряжений в пласте и 
/// автомодельное решение (радиус трещины, давление и раскрытие) для жидкости
/// заданной реологии из соответствующих текстовых файлов
bool setInitialData(
    const double bn,
    double &xStar0,
    std::vector< std::vector<double> > &A,
    const std::string pathToBarriersFile,
    std::vector< std::vector<double> > &barriers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
){
    double approximateBn = 1;
    
    const std::vector< std::vector<double> > bnIntervals = {
        {0, 0.2, 0.7296},
        {0.2, 0.4, 0.7236},
        {0.4, 0.6, 0.7162},
        {0.6, 1, 0.6975754},
        {1, 1, 0.6975754}
    };
    
    ai::assignFromVectorByIntervalCondition(approximateBn, xStar0, bn, 
        bnIntervals);

    std::cout << "Reading initial conditions from text files:" << std::endl;
    
    try{
        std::cout << "AS" + ai::string(approximateBn) + "... ";
        
        ai::parseFileInMatrix("./InitialConditions/A" 
            + ai::string(approximateBn) + ".txt", '\t', A);
        
        std::cout << "OK. Size: " << A.size() << "." << std::endl;
    }catch(const std::exception &error){
        std::cout << "Fail!" << std::endl;
        
        std::cerr << "Cannot open ./InitialConditions/A"
            << ai::string(approximateBn) + ".txt, provide file and restart." 
            << std::endl;
        
        std::cerr << error.what() << std::endl;
        
        return false;
    }
    
    if(std::string() != pathToBarriersFile){
        try{
            std::cout << "Barriers... ";
            
            ai::parseFileInMatrix(pathToBarriersFile, ' ', barriers);
            
            std::cout << "OK. Size: " << barriers.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;
            
            std::cerr << "Cannot open " << pathToBarriersFile << ", "
                << "provide file and restart." << std::endl;
            
            std::cerr << error.what() << std::endl;
            
            return false;
        }
    }else{
        std::cout << "Barriers... No file." << std::endl;
    }
    
    if(std::string() != pathToInjectionFile){
        try{
            std::cout << "Injection... ";
            
            ai::parseFileInMatrix(pathToInjectionFile, ' ', injection);
            
            std::cout << "OK. Size: " << injection.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;
            
            std::cerr << "Cannot open " << pathToInjectionFile << ", "
                << "provide file and restart." << std::endl;
            
            std::cerr << error.what() << std::endl;
            
            return false;
        }
    }else{
        std::cout << "Injection... No file." << std::endl;
    }
    
    return true;
}

/// \details Функция сохраняет начальные параметры в файл формата JSON по 
/// установленному шаблону
void saveInitialDataJSON(
    const std::string filename,
    const double modelingTime,
    const std::size_t meshHeight,
    const std::size_t meshLength,
    const double cellHeight,
    const double cellLength,
    std::vector< std::vector<double> > &fluid,
    std::vector< std::vector<double> > &proppant,
    std::vector< std::vector<double> > &layers,
    const std::string tab = std::string("    ")
){
    std::ofstream output(filename + ai::string(".json"));

    if(!output.good()){
        throw std::runtime_error(
            ai::string("Exception while saving the data into the file: ") 
            + filename
        );
    }
    
    std::size_t indentLevel = 0;
    
    auto indent = [&indentLevel, &tab]() -> std::string
    {
        std::string output;
        
        for(std::size_t i = 0; i < indentLevel; ++i){
            output += tab;
        }
        
        return output;
    };
    
    output << "{" << std::endl;
    ++indentLevel;

    output << indent() << "\"model\": \"planar3D\"," << std::endl;
    output << indent() << "\"time\": " << modelingTime << "," << std::endl;
    
    output << indent() << "\"mesh\": {" << std::endl;
    ++indentLevel;
        output << indent() << "\"height\": " << meshHeight << "," << std::endl;
        output << indent() << "\"lenght\": " << meshLength << "," << std::endl;
        output << indent() << "\"cell\": {" << std::endl;
        ++indentLevel;
            output << indent() << "\"height\": " << cellHeight << "," 
                << std::endl;
            output << indent() << "\"length\": " << cellLength
                << std::endl;
        --indentLevel;
        output << indent() << "}" << std::endl;
    --indentLevel;
    output << indent() << "}," << std::endl;
    
    output << indent() << "\"fluid\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < fluid.size(); ++i){
            if(fluid[i].size() != 5){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: fluid") 
                );
            }
            output << indent() << "\"" << ai::string(i) << "\": {" 
                << std::endl;
            ++indentLevel;
                output << indent() << "\"start time\": " << fluid[i][0]
                    << "," << std::endl;
                output << indent() << "\"stop time\": " << fluid[i][1]
                    << "," << std::endl;
                output << indent() << "\"injection\": " << fluid[i][2]
                    << "," << std::endl;
                output << indent() << "\"rheology index\": " << fluid[i][3]
                    << "," << std::endl;
                output << indent() << "\"viscosity\": " << fluid[i][4]
                    << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(fluid.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}," << std::endl;
    
    output << indent() << "\"proppant\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < proppant.size(); ++i){
            if(proppant[i].size() != 3){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: proppant") 
                );
            }
            output << indent() << "\"" << ai::string(i) << "\": {" 
                << std::endl;
            ++indentLevel;
                output << indent() << "\"start time\": " << proppant[i][0]
                    << "," << std::endl;
                output << indent() << "\"stop time\": " << proppant[i][1]
                    << "," << std::endl;
                output << indent() << "\"injection\": " << proppant[i][2] 
                    << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(fluid.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}," << std::endl;
    
    output << indent() << "\"layers\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < layers.size(); ++i){
            if(layers[i].size() != 5){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: layers") 
                );
            }
            output << indent() << "\"" << ai::string(i) << "\": {" 
                << std::endl;
            ++indentLevel;
                output << indent() << "\"y1\": " << layers[i][0] << ","
                    << std::endl;
                output << indent() << "\"y2\": " << layers[i][1] << ","
                    << std::endl;
                output << indent() << "\"plane strain modulus\": "
                    << layers[i][2] << "," << std::endl;
                output << indent() << "\"leak-off\": " << layers[i][3] << ","
                    << std::endl;
                output << indent() << "\"stress\": " << layers[i][4]
                    << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(layers.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}" << std::endl;
    
    output << "}";
}
