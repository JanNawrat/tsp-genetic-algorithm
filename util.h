#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <limits>
#include <algorithm>
#include "/opt/homebrew/Cellar/libomp/17.0.6/include/omp.h"

void shuffle(std::vector< int >& P);

int getWeight(std::vector< int >& p);

int getWeightChange(std::vector< int >& p, int i, int j);

void invert(std::vector< int >& p, int i, int j);

void initMarsenneTwister();

void initAdjacencyMatrix(int number);

void printVector(std::vector<int> &p, int len);

int binarySearch(std::vector<int> &v, int element);