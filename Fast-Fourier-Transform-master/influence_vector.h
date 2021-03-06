#pragma once
#include "ai.hh"

//generate vector
void generateVector(std::vector<long double> &vector, size_t length) 
{
	vector.resize(length);

	std::generate(vector.begin(), vector.end(), std::rand);

	for (size_t i = 0; i < length; ++i) {
		vector[i] /= std::pow(10, 2);
	}
}
//multiplication martix and vector directly
void multiply(
	std::vector< std::vector<long double> > &matrix,
	std::vector<long double> &vector,
	std::vector<long double> &result)
{
	result.resize(matrix.size());

	for (size_t i = 0; i < matrix.size(); i++)
	{
		result[i] = 0.L;

		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			result[i] += matrix[i][j] * vector[j];
		}
	}
}