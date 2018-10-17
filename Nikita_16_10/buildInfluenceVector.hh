#pragma once
#include "planar3D.hh"


void buildInfluenceVector(std::vector<std::vector<double>> &InfluenceVector,
                                            const double axMax,
                                          std::vector<double>& zP,
                                          std::vector<double> &openingAtTheStart,
                                        std::vector< std::vector<double> > &barriers);
