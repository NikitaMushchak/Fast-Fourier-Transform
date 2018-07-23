#pragma once

#include <array>
#include <cmath>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sys/stat.h>

// edited by N.D. Mushchak 19.07
////////////////////////////////
//defines
typedef std::vector<double> Complex;
typedef std::vector < std::vector <double >> CArray;
//typedef std::vector<CArray> CMatrix;






//#include <dirent.h>

/*
TODO: fix
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include "include/dirent.h"
#else
#include <dirent.h>
#endif
*/

namespace ai {
	/*/ **************************************************************** /*/

	/*/
	* функция, возвращающая версию библиотеки
	/*/
	inline std::string getVersion() {
		return "1.1.0";
	}


	/*/ **************************************************************** /*/

	/*/
	* функция, приводящее значение к типу строки
	/*/
	template<typename T>
	std::string string(const T value) {
		std::ostringstream stream;

		stream << value;

		return stream.str();
	}

	/*/ **************************************************************** /*/

	/*/
	* функция
	/*/
	inline std::size_t counter(std::size_t value = 0) {
		static std::size_t count = 0;

		if (0 != value) {
			count = value;
		}

		++count;

		return count;
	}

	/*/
	* функция
	/*/
	inline std::string marker(std::size_t value = 0) {
		static std::size_t count = 0;

		if (0 != value) {
			count = value;
		}

		++count;

		return "Marker #" + ai::string(count);
	}

	/*/
	* функция
	/*/
	inline void printMarker(std::size_t value = 0) {
		std::cout << marker(value) << std::endl;
	}

	/*/
	* функция, возвращающая знак переданного числа (-1, 0, 1)
	/*/
	template<typename T>
	T sign(const T a) {
		if (0 == a) {
			return (T)0;
		}

		return copysign(1, a);
	}

	/*/
	* функция, возвращающая максимальное значение в паре
	/*/
	template<typename T>
	T min(const T a, const T b) {
		if (a > b) {
			return b;
		}

		return a;
	}

	/*/
	* функция, возвращающая максимальное значение в паре
	/*/
	template<typename T>
	T max(const T a, const T b) {
		if (a < b) {
			return b;
		}

		return a;
	}

	/*/
	* функция, возвращающая максимальное значение в векторе
	/*/
	template<typename T>
	T max(const std::vector<T> input) {
		T maximum = input[0];

		for (std::size_t i = 1; i < input.size(); ++i) {
			maximum = ai::max(maximum, input[i]);
		}

		return maximum;
	}

	/*/
	* функция, возвращающая максимальное значение в матрице
	/*/
	template<typename T>
	T max(const std::vector< std::vector<T> > input) {
		T maximum = input[0][0];

		for (std::size_t i = 0; i < input.size(); ++i) {
			for (std::size_t j = 0; j < input[0].size(); ++j) {
				maximum = ai::max(maximum, input[i][j]);
			}
		}

		return maximum;
	}

	/*/ **************************************************************** /*/

	/*/
	* функция, возвращающая текущее время
	/*/
	inline std::chrono::high_resolution_clock::time_point time() {
		return std::chrono::high_resolution_clock::now();
	}

	/*/
	* функция, возвращающая разницу между двумя значениями времени в секундах
	/*/
	inline double duration(
		const std::chrono::high_resolution_clock::time_point start,
		const std::chrono::high_resolution_clock::time_point finish,
		const std::string scale = std::string("ms")
	) {
		if (std::string("s") == scale) {
			return (double)std::chrono::duration_cast
				<std::chrono::seconds> (finish - start).count();
		}

		if (std::string("us") == scale) {
			return (double)std::chrono::duration_cast
				<std::chrono::microseconds> (finish - start).count();
		}

		return (double)std::chrono::duration_cast
			<std::chrono::milliseconds> (finish - start).count();
	}

	/*/ **************************************************************** /*/

	/*/
	* функция, проверяющая существование директории с указанным адресом
	/*/
	inline bool folderExists(const std::string name) {
		struct stat buffer;

		return 0 == stat(name.c_str(), &buffer);
	}

	/*/ **************************************************************** /*/

	/*/
	* функция, проверяющая, начинается ли строка с указанной подстроки
	/*/
	inline bool hasPrefix(const std::string &str, const std::string &prefix) {
		return str.size() >= prefix.size()
			&& 0 == str.compare(0, prefix.size(), prefix);
	}

	/*/
	* функция, проверяющая, заканчивается ли строка указанной подстрокой
	/*/
	inline bool hasSuffix(const std::string &str, const std::string &suffix) {
		return str.size() >= suffix.size() &&
			str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

	/*/
	*
	/*/
	inline std::string replace(
		std::string text,
		const std::string &substring,
		const std::string &replacement
	) {
		std::size_t position = text.find(substring);
		std::size_t substringSize = substring.size();
		std::size_t replacementSize = replacement.size();

		while (std::string::npos != position) {
			text.replace(position, substringSize, replacement);
			position = text.find(substring, position + replacementSize);
		}
		return text;
	}

	/*/
	*
	/*/
	inline std::string parseParameter(
		const char *input,
		const std::string name
	) {
		std::string parameter(input);
		if (hasPrefix(parameter, name)) {
			parameter = parameter.substr(name.size());
			return parameter;
		}
		return NULL;
	}

	/*/
	* функция, присваивающая значение в зависимости от значения указанного
	* параметра
	/*/
	template<typename T>
	inline void assignFromVectorByIntervalCondition(
		T &value,
		const T parameter,
		const std::vector< std::vector<T> > intervals
	) {
		for (std::size_t i = 0; i < intervals.size(); ++i) {
			if (intervals[i][0] <= parameter && parameter <= intervals[i][1]) {
				value = intervals[i][2];

				return;
			}
		}
	}

	/*/
	* функция, присваивающая пару значений в зависимости от значения указанного
	* параметра
	/*/
	template<typename T>
	inline void assignFromVectorByIntervalCondition(
		T &firstValue,
		T &secondValue,
		const T parameter,
		const std::vector< std::vector<T> > intervals
	) {
		for (std::size_t i = 0; i < intervals.size(); ++i) {
			if (intervals[i][0] <= parameter && parameter <= intervals[i][1]) {
				firstValue = intervals[i][1];
				secondValue = intervals[i][2];

				return;
			}
		}
	}

	/*/
	* функция
	/*/
	inline bool assignBooleanParameter(
		const char *input,
		const std::string name,
		bool &value
	) {
		if (name == std::string(input)) {
			value = true;

			return true;
		}

		return false;
	}

	/*/
	* функция
	/*/
	inline bool assignCharParameter(
		const char *input,
		const std::string name,
		char &value
	) {
		std::string parameter = std::string(input);

		if (ai::hasPrefix(parameter, name)) {
			parameter = parameter.substr(name.size());

			if (std::string() != parameter) {
				value = (char)parameter[0];
			}
			else {
				value = ' ';
			}

			return true;
		}

		return false;
	}

	/*/
	* функция, присваивающая строковое значение в результате парсинга
	* указанного текстового параметра
	/*/
	inline bool assignStringParameter(
		const char *input,
		const std::string name,
		std::string &value
	) {
		//TODO: parseParameter
		std::string parameter = std::string(input);

		if (ai::hasPrefix(parameter, name)) {
			value = parameter.substr(name.size());

			return true;
		}

		return false;
	}

	/*/
	* TODO
	/*/
	template<typename T>
	inline bool assignParameter(
		const char *input,
		const std::string name,
		T &value
	) {
		std::string parameter = std::string(input);

		if (ai::hasPrefix(parameter, name)) {
			parameter = parameter.substr(name.size());

			if (std::istringstream(parameter) >> value) {
				return true;
			}
		}

		return false;
	}

	/*/
	* функция, присваивающая вещественное значение по модулю в результате
	* парсинга указанного текстового параметра
	/*/
	inline bool assignAbsDoubleParameter(
		const char *input,
		const std::string name,
		double &value
	) {
		std::string parameter = std::string(input);

		if (ai::hasPrefix(parameter, name)) {
			parameter = parameter.substr(name.size());
			value = std::abs(strtod(parameter.c_str(), nullptr));

			return true;
		}

		return false;
	}

	/*/
	* функция, выводящая текстовый индикатор прогресса шириной в 80 символов
	/*/
	inline void showProgressBar(double progress) {
		if (1 < progress) {
			progress = 1;
		}

		if (0.01 > progress) {
			progress = 0;
		}

		int width = progress * 73;

		std::cout << std::fixed;
		std::cout.precision(1);
		std::cout.flush();

		std::cout << "\r" << std::string(width, '=')
			<< std::string(73 - width, '-')
			<< " " << progress * 100. << "%";
	}

	/*/ **************************************************************** /*/

	/*/
	* функция, выводящая значения строки на экран
	/*/
	inline void printLine(const std::string line) {
		std::cout << "\r" << std::setw(80) << std::left << line << std::endl;
	}

	/*/ **************************************************************** /*/

	/*/
	* функция, считывающая значения из текстового файла в матрицу
	/*/
	template<typename T>
	inline void parseFileInMatrix(
		const std::string filename,
		const char separator,
		std::vector< std::vector<T> > &matrix
	) {
		std::ifstream input(filename);

		if (!input.good()) {
			throw std::runtime_error(
				ai::string("Exception while parsing the file into a matrix: ")
				+ filename
			);
		}

		std::string token;

		T value;

		for (std::string line; getline(input, line);) {
			if ('#' == line[0]) {
				continue;
			}

			std::istringstream stream(line);

			std::vector<double> row;

			while (std::getline(stream, token, separator)) {
				if (std::istringstream(token) >> value) {
					row.push_back(value);
				}
			}

			matrix.push_back(row);
		}
	}

	/*/
	* функция, считывающая значения из текстового файла в вектор
	/*/
	template<typename T>
	inline void parseFileInVector(
		const std::string filename,
		const char separator,
		std::vector<T> &vector
	) {
		std::ifstream input(filename);

		if (!input.good()) {
			throw std::runtime_error(
				ai::string("Exception while parsing the file into a vector: ")
				+ filename
			);
		}

		std::string token;

		T value;

		if ('\n' == separator) {
			for (std::string line; getline(input, line);) {
				if ('#' == line[0]) {
					continue;
				}

				if (std::istringstream(line) >> value) {
					vector.push_back(value);
				}
			}
		}
		else {
			std::string line;

			std::getline(input, line);

			std::istringstream stream(line);

			while (std::getline(stream, token, separator)) {
				if (std::istringstream(token) >> value) {
					vector.push_back(value);
				}
			}
		}
	}

	/*/
	* функция
	/*/
	template<typename T>
	void accumulateFileInMatrix(
		const std::string filename,
		const char separator,
		std::vector< std::vector<T> > &matrix
	) {
		std::vector< std::vector<T> > matrixToAdd;

		ai::parseFileInMatrix(filename, separator, matrixToAdd);

		for (std::size_t i = 0; i < matrix.size(); ++i) {
			std::transform(
				matrix[i].begin(),
				matrix[i].end(),
				matrixToAdd[i].begin(),
				matrix[i].begin(),
				std::plus<T>()
			);
		}
	}


	/*/
	* функция
	/*/
	template<typename T>
	void accumulateFileInVector(
		const std::string filename,
		const char separator,
		std::vector<T> &vector
	) {
		std::vector<T> vectorToAdd;

		ai::parseFileInVector(filename, separator, vectorToAdd);

		std::transform(
			vector.begin(),
			vector.end(),
			vectorToAdd.begin(),
			vector.begin(),
			std::plus<T>()
		);
	}

	/*/
	* функция, выводящая значения матрицы на экран
	/*/
	template<typename T>
	inline void printMatrix(
		const std::vector<std::vector <T> > matrix,
		const int precision = 5
	) {
		std::cout << std::scientific;
		std::cout.precision(precision);

		if (1 > matrix.size()) {
			throw std::runtime_error(
				ai::string("Exception while printing the matrix: ")
				+ ai::string("size should be at least 1.")
			);
		}

		std::cout << "Matrix[" << matrix.size()
			<< "x" << matrix[0].size() << "] = {" << std::endl;

		for (const std::vector<T> &vector : matrix) {
			std::size_t lastIndex = vector.size() - 1;

			for (std::size_t i = 0; i < lastIndex; ++i) {
				std::cout << vector[i] << ", ";
			}
			std::cout << vector[lastIndex] << std::endl;
		}

		std::cout << "}[" << matrix.size()
			<< "x" << matrix[0].size() << "]" << std::endl;
	}

	/*/
	* функция, выводящая значения вектора на экран
	/*/
	template<typename T>
	inline void printVector(
		const std::vector<T> vector,
		const int precision = 5
	) {
		std::size_t lastIndex = vector.size() - 1;

		std::cout << std::scientific;
		std::cout.precision(precision);

		std::cout << "Vector[" << vector.size() << "] = {" << std::endl;

		for (std::size_t i = 0; i < lastIndex; ++i) {
			std::cout << vector[i] << ", ";
		}
		std::cout << vector[lastIndex] << std::endl;

		std::cout << "}[" << vector.size() << "]" << std::endl;
	}

	/*/
	* функция, сохраняющая значения матрицы в текстовый файл
	/*/
	template<typename T>
	inline void saveMatrix(
		const std::string filename,
		const std::vector<std::vector <T> > matrix,
		std::string comment = std::string(),
		std::string type = std::string("text"),
		std::string delimiter = std::string(" ")
	) {
		std::string extension("_m.txt");
		std::string prefix("");
		std::string suffix("");

		if (std::string("wolfram") == type) {
			extension = std::string("_m.wm");
			prefix = std::string("{");
			delimiter = std::string(", ");
			suffix = std::string("}");
		}

		if (std::string("excel") == type) {
			extension = std::string("_m.csv");
			delimiter = std::string("; ");
		}

		if (std::string("data") == type) {
			extension = std::string("_m.dat");
			delimiter = std::string("\t");
		}

		std::ofstream output(filename + extension);

		if (!output.good()) {
			throw std::runtime_error(
				ai::string("Exception while saving the matrix into the file: ")
				+ filename
			);
		}

		if (std::string() != comment) {
			output << comment << std::endl;
		}

		output << prefix;

		//delimiter after each line
		//what about setw?
		for (const std::vector<T> &row : matrix) {
			output << prefix;

			const std::size_t lastIndex = row.size() - 1;

			for (std::size_t i = 0; i < lastIndex; ++i) {
				output << std::setw(14) << row[i] << delimiter;
			}

			output << std::setw(14) << row[lastIndex] << suffix << std::endl;
		}

		output << suffix;
	}

	/*/
	* функция, сохраняющая значения вектора в текстовый файл
	/*/
	template<typename T>
	inline void saveVector(
		const std::string filename,
		const std::vector <T> vector,
		std::string comment = std::string(),
		std::string type = std::string("text"),
		std::string delimiter = std::string("\n")
	) {
		std::string extension("_v.txt");
		std::string prefix("");
		std::string suffix("");

		if (std::string("wolfram") == type) {
			extension = std::string("_v.wm");
			prefix = std::string("{");
			delimiter = std::string(", ");
			suffix = std::string("}");
		}

		if (std::string("excel") == type) {
			extension = std::string("_v.csv");
			delimiter = std::string("; ");
		}

		if (std::string("data") == type) {
			extension = std::string("_v.dat");
			delimiter = std::string("\t");
		}

		std::ofstream output(filename + extension);

		if (!output.good()) {
			throw std::runtime_error(
				ai::string("Exception while saving the vector into the file: ")
				+ filename
			);
		}

		if (std::string() != comment) {
			output << comment << std::endl;
		}

		output << prefix;

		const std::size_t lastIndex = vector.size() - 1;

		for (std::size_t i = 0; i < lastIndex; ++i) {
			output << vector[i] << delimiter;
		}

		output << vector[lastIndex] << suffix;
	}

	/*/
	* функция, сохраняющая значения текстовой строки в текстовый файл
	/*/
	template<typename T>
	inline void saveLine(
		const std::string filename,
		const std::string line,
		std::string comment = std::string()
	) {
		std::ofstream output(filename + "_l.txt");

		if (!output.good()) {
			throw std::runtime_error(
				ai::string("Exception while saving the line into the file: ")
				+ filename
			);
		}

		if (std::string() != comment) {
			output << comment << "\n";
		}

		output << line;
	}

	/*/ **************************************************************** /*/

	/*/
	* функция
	/*/
	inline std::size_t countLinesInFile(
		const std::string filename,
		const std::string token = std::string()
	) {
		std::ifstream input(filename);

		if (!input.good()) {
			throw std::runtime_error(
				ai::string("Exception while counting lines in the file: ")
				+ filename
			);
		}

		std::size_t count = 0;

		if (std::string() == token) {
			for (std::string line; getline(input, line);) {
				++count;
			}
		}
		else {
			for (std::string line; getline(input, line);) {
				if (std::string::npos != line.find(token)) {
					++count;
				}
			}
		}

		return count;
	}

	///*/
	// * функция 
	///*/
	//inline std::vector<std::string> listFilesWithExtension(
	//    const std::string path,
	//    const std::string extension,
	//    const std::string prefix = std::string()
	//){
	//    DIR *dir;

	//    struct dirent *ent;

	//    dir = opendir(path.c_str());

	//    std::vector<std::string> files;

	//    if(NULL != dir){
	//        while(NULL != (ent = readdir (dir))){
	//            if(DT_REG == ent->d_type && hasSuffix(ent->d_name, extension)){
	//                files.push_back(prefix + ent->d_name);
	//            }
	//        }

	//        closedir(dir);
	//    }

	//	return files;
	//}

	///*/ **************************************************************** /*/
	//
	///*/
	// * функция 
	///*/
	//inline std::string execute(const std::string command){
	//    std::array<char, 128> buffer;

	//    std::string result;

	//    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);

	//    if(!pipe){
	//        throw std::runtime_error(
	//            ai::string("Exception while executing the command: ")
	//            + command
	//        );
	//    }

	//	while(!feof(pipe.get())){
	//        if(nullptr != fgets(buffer.data(), 128, pipe.get())){
	//            result += buffer.data();
	//        }
	//    }

	//    return result;
	//}
}