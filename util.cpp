#include "util.h"

std::string files[15] = {"xqf131", "xqg237", "pma343", "pka379", "bcl380", "pbl395", "pbk411", "pbn423", "pbm436", "xql662", "xit1083", "icw1483", "djc1785", "dcb2086", "pds2566"};

std::mt19937 mt;

std::vector< std::vector< int > > M;

void shuffle(std::vector< int >& P) {
    for(int i = P.size() - 1; i > 0; i--) {
        std::uniform_int_distribution<> dist(0, i);
        std::swap(P[dist(mt)], P[i]);
    }
}

int getWeight(std::vector< int >& p) {
    int weight = 0;
    for(int i = 1 ; i < p.size(); i++) {
        weight += M[p[i - 1]][p[i]];
    }
    weight += M[p[p.size()-1]][p[0]];
    return weight;
}

int getWeightChange(std::vector< int >& p, int i, int j) {
    int change = 0;
    int dim = p.size();
    if(i == 0) {
        change += M[p[j]][p[dim-1]] - M[p[i]][p[dim-1]];
        change += M[p[i]][p[j+1]] -  M[p[j]][p[j+1]];
    } else if(j == dim - 1) {
        change += M[p[0]][p[i]] - M[p[0]][p[j]];
        change += M[p[i-1]][p[j]] - M[p[i-1]][p[i]];
    } else {
        change += M[p[i-1]][p[j]] - M[p[i-1]][p[i]];
        change += M[p[i]][p[j+1]] -  M[p[j]][p[j+1]];
    }
    return change;
}

void invert(std::vector< int >& p, int i, int j) {
    for(int k = i; k <= (j + i - 1)/2; k++) {
        std::swap(p[k], p[j + i - k]);
    }
}

void initMarsenneTwister() {
    std::random_device rd;
    mt.seed(rd());
}

void initAdjacencyMatrix(int number) {
    std::fstream input;
    std::fstream output;
    std::string filename = files[number];
    std::cout << "Reading from: " << filename << '\n';
    input.open("data/"+filename+".tsp", std::ios::in);
    // output.open("data/annRes/"+filename+".txt", std::ios::out);
    if(!(input.good())||!(output.good())) {
        exit(1);
    }

    std::string data = "";
    while(data != "DIMENSION") {
        input >> data;
    }
    // dim - number of vertices
    // n - number of current vertex
    // x, y - coordinates of current vertex
    int dim, n, x, y;
    input >> data >> dim;
    std::vector< double > Vx(dim);
    std::vector< double > Vy(dim);
    while(data != "NODE_COORD_SECTION") {
        input >> data;
    }
    for(int i = 0; i < dim; i++) {
        input >> n >> x >> y;
        Vx[i] = x;
        Vy[i] = y;
    }

    // create adjacency matrix
    M = std::vector< std::vector<int> >(dim);
    for(int i = 0; i < dim; i++) {
        M[i] = std::vector<int>(dim);
    }

    // filling adjacency matrix
    for(int i = 0; i < dim; i++) {
        M[i][i] = 0;
        for(int j = i + 1; j < dim; j++) {
            M[i][j] = std::round( sqrt( pow(Vx[i] - Vx[j], 2) + pow(Vy[i] - Vy[j], 2) ) );
            M[j][i] = M[i][j];
        }
    }
}

void printVector(std::vector<int> &p, int len) {
    len = len > 0 ? len : p.size();
    for (int i = 0; i < len; i++) {
        std::cout << p[i] << "\t";
    }
    std::cout << '\n';
}

int binarySearch(std::vector<int> &v, int element) {
    int start = 0;
    int end = v.size() - 1;

    while (start <= end) {
        int midpoint = start + (end - start) / 2;
        if (v[midpoint] == element)
            return midpoint;
        if (v[midpoint] < element)
            start = midpoint + 1;
        else
            end = midpoint - 1;
    }
    return -1;
}