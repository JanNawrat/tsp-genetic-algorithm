#pragma once

#include "util.h"

extern std::mt19937 mt;

class Population {
    public:
    struct chromosome {
        std::vector< int > geneticInformation;
        int weight;
        double fitness;
    };
    int n;
    int dim;
    std::vector< chromosome > chromosomes;
    int minimalWeight;
    double totalFitness;
    std::uniform_real_distribution<> distribution;

    Population(int n, int dim) {
        this->n = n;
        this->dim = dim;
        chromosomes = std::vector< chromosome >(n);
        for (int i = 0; i < n; i++) {
            chromosomes[i] = {std::vector< int >(dim), 0};
            for (int j = 0; j < dim; j++) {
                chromosomes[i].geneticInformation[j] = j;
            }
            shuffle(chromosomes[i].geneticInformation);
        }
        computePartialWeights();
        computeFitness();
        findMinimalWeight();
    }

    void computePartialWeights() {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            chromosomes[i].weight = getWeight(chromosomes[i].geneticInformation);
        }
    }

    void computeFitness() {
        totalFitness = 0;
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            chromosomes[i].fitness = 1.0 / chromosomes[i].weight;
        }
        for (int i = 0; i < n; i++) {
            totalFitness += chromosomes[i].fitness;
        }
        distribution = std::uniform_real_distribution<>(0, totalFitness);
    }

    void findMinimalWeight() {
        minimalWeight = INT32_MAX;
        for (int i = 0; i < n; i++) {
            if (chromosomes[i].weight < minimalWeight) {
                minimalWeight = chromosomes[i].weight;
            }
        }
    }

    std::vector< int > rouletteSelection() {
        std::vector< double > positions(n);
        for (int i = 0; i < n; i++) {
            positions[i] = distribution(mt);
        }
        sort(positions.begin(), positions.end());

        std::vector< int > selectionResult(n);
        int positionsIterator = 0;
        int chromosomesIterator = 0;
        double fitnessSum = chromosomes[0].fitness;
        while (positionsIterator < n) {
            while (positions[positionsIterator] > fitnessSum) {
                chromosomesIterator += 1;
                fitnessSum += chromosomes[chromosomesIterator].fitness;
            }
            selectionResult[positionsIterator] = chromosomesIterator;
            positionsIterator += 1;
        }

        shuffle(selectionResult);
        return selectionResult;
    }

    void PMX(std::vector< int > &parents, double probability) {
        std::uniform_int_distribution<> mutdist(0,n-1);
        std::vector< chromosome > newChromosomes(n);
        std::uniform_int_distribution<> pmxdist(0, dim-1);
        #pragma omp parallel for
        for(int k = 0; k < n; k+=2) {
            if(mt() > probability) {
                newChromosomes[k].geneticInformation = chromosomes[parents[k]].geneticInformation;
                newChromosomes[k+1].geneticInformation = chromosomes[parents[k+1]].geneticInformation;
                continue;
            }

            // preparing child vectors
            newChromosomes[k] = {std::vector< int >(dim), 0};
            newChromosomes[k+1] = {std::vector< int >(dim), 0};

            // references
            std::vector< int > &child1 = newChromosomes[k].geneticInformation;
            std::vector< int > &child2 = newChromosomes[k+1].geneticInformation;
            std::vector< int > &parent1 = chromosomes[parents[k]].geneticInformation;
            std::vector< int > &parent2 = chromosomes[parents[k+1]].geneticInformation;

            for (int j = 0; j < dim; j++) {
                child1[j] = -1;
                child2[j] = -1;
            }

            // selecting crossover points
            int crossoverPoint1;
            int crossoverPoint2;
            do {
                crossoverPoint1 = pmxdist(mt);
                crossoverPoint2 = pmxdist(mt);
            } while (crossoverPoint1 == crossoverPoint2);
            if (crossoverPoint1 > crossoverPoint2) {
                std::swap(crossoverPoint1, crossoverPoint2);
            }

            // vectors for remembering copied genes
            std::vector<int> copiedGenes1(crossoverPoint2-crossoverPoint1+1);
            std::vector<int> copiedGenes2(crossoverPoint2-crossoverPoint1+1);

            // copy selected part from parent1 to child
            // seg start
            for (int i = crossoverPoint1; i <= crossoverPoint2; i++) {
                // std::cout << crossoverPoint1 << " " << crossoverPoint2 << " " << i << '\n';
                // std::cout << 1 << std::flush;
                child1[i] = parent1[i];
                // std::cout << 2 << std::flush;
                copiedGenes1[i - crossoverPoint1] = parent1[i];
                child2[i] = parent2[i];
                // std::cout << 3 << std::flush;
                copiedGenes2[i - crossoverPoint1] = parent2[i];
                // std::cout << 4 << std::endl;
            }
            // seg end
            sort(copiedGenes1.begin(), copiedGenes1.end());
            sort(copiedGenes2.begin(), copiedGenes2.end());

            // copy genes from second parent
            for (int i = crossoverPoint1; i <= crossoverPoint2; i++) {
                if (binarySearch(copiedGenes1, parent2[i]) == -1) {
                    int j;
                    int l = i;
                    do {
                        j = 0;
                        while(parent2[j] != child1[l]) {
                            j += 1;
                        }
                        l = j;
                    } while (child1[j] != -1);

                    child1[j] = parent2[i];
                }
                if (binarySearch(copiedGenes2, parent1[i]) == -1) {
                    int j;
                    int l = i;
                    do {
                        j = 0;
                        while(parent1[j] != child2[l]) {
                            j += 1;
                        }
                        l = j;
                    } while (child2[j] != -1);

                    child2[j] = parent1[i];
                }
            }

            // copy the rest of the genes
            for (int i = 0; i < crossoverPoint1; i++) {
                if (child1[i] == -1) {
                    child1[i] = parent2[i];
                }
                if (child2[i] == -1) {
                    child2[i] = parent1[i];
                }
            }

            for (int i = crossoverPoint2+1; i < dim; i++) {
                if (child1[i] == -1) {
                    child1[i] = parent2[i];
                }
                if (child2[i] == -1) {
                    child2[i] = parent1[i];
                }
            }
        }
        chromosomes = newChromosomes;
    }

    void mutation(double probability) {
        std::uniform_int_distribution<> mutdist(0,n-1);
        #pragma omp parallel for
        for (int k = 0; k < n; k++) {
            if(mt() > probability) {
                continue;
            }
            int i = mutdist(mt);
            int j = mutdist(mt);
            std::swap(chromosomes[k].geneticInformation[i],chromosomes[k].geneticInformation[j]);
        }
    }
    
    void memeticTabuSearch(int duration) {
        #pragma omp parallel for
        for (int k = 0; k < n; k++) {
            std::vector< int > currentChromosome = chromosomes[k].geneticInformation;
            int currentWeight = chromosomes[k].weight;
            std::vector< std::vector< int > > tabuList(duration);
            int minI, minJ;

            for (int t = 0; t < duration; t++) {
                int bestWeightImprovement = INT32_MAX;

                
                for(int l = 0; l < n; l++) {
                    std::uniform_int_distribution<> distI(0, dim-2);
                    int i = distI(mt);
                    std::uniform_int_distribution<> distJ(i + 1, dim-1);
                    int j = distJ(mt);

                    // inverted - same weight
                        if(j - i + 1 == dim) {
                            continue;
                        }
                        int weightChange = getWeightChange(currentChromosome, i, j);

                        if(weightChange < bestWeightImprovement) {
                            invert(currentChromosome, i, j);
                            if (std::find(tabuList.begin(), tabuList.end(), currentChromosome) == tabuList.end()) {
                                bestWeightImprovement = weightChange;
                                minI = i;
                                minJ = j;
                            }
                            invert(currentChromosome, i, j);
                        }
                }

                invert(currentChromosome, minI, minJ);
                currentWeight += bestWeightImprovement;

                // checking for improvement
                if(currentWeight < chromosomes[k].weight) {
                    chromosomes[k].weight = currentWeight;
                    chromosomes[k].geneticInformation = currentChromosome;
                }

                // adding new result to tabu list
                tabuList[t] = currentChromosome;
            }
        }
    }

    void sortByWeight() {
        std::sort(chromosomes.begin(), chromosomes.end(), [&](chromosome A, chromosome B) -> bool {return A.weight < B.weight;});
    }
};

void operations(Population *p) {
    for (int i = 0; i < 50; i++) {
        // selection
        p->computeFitness();
        std::vector< int > selection = p->rouletteSelection();
        // crossover
        p->PMX(selection, 0.7);
        // mutation
        p->mutation(0.3);
        p->computePartialWeights();
        // memetic
        p->memeticTabuSearch(20);
        // end
    }
}

int populationBasedAlgorithm(int numberOfIslands, int populationCount, int migrationCount, int dim) {
    std::vector< Population* > populations(numberOfIslands);
    std::vector< std::vector< Population::chromosome > > migrationVector(numberOfIslands);
    for (int i = 0; i < numberOfIslands; i++) {
        populations[i] = new Population(populationCount, dim);
        migrationVector[i] = std::vector< Population::chromosome >(migrationCount);
    }
    long minimalWeight = INT64_MAX;
    int iterations = 0;
    while (iterations < 4) {
        for(int i = 0; i < numberOfIslands; i++) {
            iterations += 1;
            operations(populations[i]);
            populations[i]->findMinimalWeight();
            if (populations[i]->minimalWeight < minimalWeight) {
                minimalWeight = populations[i]->minimalWeight;
                iterations = 0;
            }
            // migration
            for (int i = 0; i < numberOfIslands; i++) {
                populations[i]->sortByWeight();
                for (int j = 0; j < migrationCount; j++) {
                    migrationVector[i][j] = populations[i]->chromosomes[j];
                }
            }
            for (int i = 0; i < numberOfIslands; i++) {
                for (int j = 0; j < migrationCount; j++) {
                    int previous = i > 0 ? i - 1 : numberOfIslands - 1;
                    populations[i]->chromosomes[j] = migrationVector[previous][j];
                }
            }
        }
    }
    return minimalWeight;
}
