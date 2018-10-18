#pragma once
#include "planar3D.hh"


void buildInfluenceVector(std::vector<double> &InfluenceVector1,
                            std::vector<double> &InfluenceVector2,
                                            const double axMax,
                                          std::vector<double>& zP,
                                          std::vector<double> &openingAtTheStart,
                                        std::vector< std::vector<double> > &barriers);
