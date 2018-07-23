#pragma once
#include "ai.hh"

//generate vector
void generateVector(std::vector<double> &vector, size_t length)
{
	vector.resize(length);

	std::generate(vector.begin(), vector.end(), std::rand);

	for (size_t i = 0; i < length; ++i) {
		vector[i] /= std::pow(10, 1);
	}
}
//multiplication martix and vector directly
void multiply(
	std::vector< std::vector<double> > &matrix,
	std::vector<double> &vector,
	std::vector<double> &result)
{
	result.resize(matrix.size());

	for (size_t i = 0; i < matrix.size(); i++)
	{
		result[i] = 0.;

		/*for (size_t j = 0; j < matrix[i].size(); j++)
		{
			result[i] += matrix[i][j] * vector[j];
		}*/
        
        size_t j = 0;

        for(; j <= vector.size() - 4; j += 4){
            result[i] += matrix[i][j] * vector[j]
                + matrix[i][j + 1] * vector[j + 1]
                + matrix[i][j + 2] * vector[j + 2]
                + matrix[i][j + 3] * vector[j + 3];
        }

        for(; j < vector.size(); ++j){
            result[i] += matrix[i][j] * vector[j];
        }
	}
}