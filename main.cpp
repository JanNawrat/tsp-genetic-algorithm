#include "util.h"
#include "algorithm.h"

extern std::vector< std::vector< int > > M;

int main(int argc, char **argv) {
    if(argc != 2) {
        return 1;
    }
    initAdjacencyMatrix(atoi(argv[1]) - 1);
    int dim = M.size();
    initMarsenneTwister();
    //for (int i = 0; i < 20; i++) {
        int minimalWeight = populationBasedAlgorithm(5, 200, 20, dim);
        std::cout << minimalWeight << '\n';
    //}
}
