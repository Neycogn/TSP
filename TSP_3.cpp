#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <chrono>
#include <fstream>
#include <sstream>

using namespace std;

struct City {
    int index;
    double x, y;
};

double distance(const City &a, const City &b) {
    double xd = a.x - b.x;
    double yd = a.y - b.y;
    return round(sqrt((xd * xd + yd * yd) / 10.0));  // Chia cho 10 d? gi?m don v?
}

double calculateTotalDistance(const vector<City> &cities, const vector<int> &tour) {
    double totalDistance = 0.0;
    for (int i = 0; i < tour.size() - 1; ++i) {
        totalDistance += distance(cities[tour[i] - 1], cities[tour[i + 1] - 1]);
    }
    totalDistance += distance(cities[tour.back() - 1], cities[tour.front() - 1]);
    return totalDistance;
}

vector<int> createTour(int numCities) {
    vector<int> tour(numCities);
    iota(tour.begin(), tour.end(), 1);
    random_shuffle(tour.begin(), tour.end());
    return tour;
}

vector<vector<int>> initializePopulation(int populationSize, const vector<int> &initialTour) {
    vector<vector<int>> population;
    for (int i = 0; i < populationSize; ++i) {
        if (i == 0) {
            population.push_back(initialTour);  // Qu?n th? ban d?u
        } else {
            population.push_back(createTour(initialTour.size()));
        }
    }
    return population;
}

vector<int> selectParent(const vector<vector<int>> &population, const vector<double> &fitness) {
    double totalFitness = accumulate(fitness.begin(), fitness.end(), 0.0);
    double randomValue = (double) rand() / RAND_MAX * totalFitness;
    double runningSum = 0.0;

    for (int i = 0; i < population.size(); ++i) {
        runningSum += fitness[i];
        if (runningSum >= randomValue) {
            return population[i];
        }
    }
    return population.back();
}

vector<int> crossover(const vector<int> &parent1, const vector<int> &parent2) {
    int start = rand() % parent1.size();
    int end = start + rand() % (parent1.size() - start);
    
    vector<int> child(parent1.size(), -1);
    for (int i = start; i <= end; ++i) {
        child[i] = parent1[i];
    }

    int currentPos = 0;
    for (int i = 0; i < parent2.size(); ++i) {
        if (find(child.begin(), child.end(), parent2[i]) == child.end()) {
            while (child[currentPos] != -1) {
                currentPos++;
            }
            child[currentPos] = parent2[i];
        }
    }
    return child;
}

void mutate(vector<int> &tour, double mutationRate) {
    for (int i = 0; i < tour.size(); ++i) {
        if ((double) rand() / RAND_MAX < mutationRate) {
            int j = rand() % tour.size();
            swap(tour[i], tour[j]);
        }
    }
}

void saveToFile2(const string &filename, const vector<int> &bestTour, double bestFitness, double bestDistance, int numGenerations, 
                 const string &startTime, const string &endTime, double elapsedTime) {
    ofstream file(filename);
    if (file.is_open()) {
        file << "L?i gi?i t?t nh?t: ";
        for (int city : bestTour) {
            file << city << " ";
        }
        file << endl;
        file << "Giá tr? thích nghi c?a l?i gi?i: " << bestFitness << endl;
        file << "Giá tr? t?i uu c?a l?i gi?i: " << bestDistance << endl;
        file << "S? lu?ng th? h?: " << numGenerations << endl;
        file << "Th?i gian b?t d?u: " << startTime << endl;
        file << "Th?i gian k?t thúc: " << endTime << endl;
        file << "Th?i gian ch?y thu?t toán (s): " << elapsedTime << endl;
        file.close();
    }
}

void saveToFile3(const string &filename, const vector<vector<int>> &bestToursPerGeneration, 
                 const vector<double> &bestFitnessesPerGeneration, const vector<double> &bestDistancesPerGeneration) {
    ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < bestToursPerGeneration.size(); ++i) {
            file << "Th? h?: " << i + 1 << endl;
            file << "L?i gi?i t?t nh?t: ";
            for (int city : bestToursPerGeneration[i]) {
                file << city << " ";
            }
            file << endl;
            file << "Giá tr? thích nghi c?a l?i gi?i: " << bestFitnessesPerGeneration[i] << endl;
            file << "Giá tr? t?i uu c?a l?i gi?i: " << bestDistancesPerGeneration[i] << endl;
        }
        file.close();
    }
}

void loadCitiesAndInitialTour(const string &filename, vector<City> &cities, vector<int> &initialTour) {
    ifstream file(filename);
    int numCities;
    file >> numCities;

    cities.resize(numCities);
    for (int i = 0; i < numCities; ++i) {
        cities[i].index = i + 1;
        file >> cities[i].x >> cities[i].y;
    }

    // Ð?c qu?n th? ban d?u
    int city;
    while (file >> city && city != -1) {
        initialTour.push_back(city);
    }
}

int main() {
    srand(time(0));

    vector<City> cities;
    vector<int> initialTour;
    loadCitiesAndInitialTour("input1.txt", cities, initialTour);  // Ð?c d? li?u t? file input.txt

    const int POPULATION_SIZE = 100;
    const int NUM_GENERATIONS = 1000;
    const double MUTATION_RATE = 0.01;

    vector<vector<int>> population = initializePopulation(POPULATION_SIZE, initialTour);
    vector<vector<int>> bestToursPerGeneration;
    vector<double> bestFitnessesPerGeneration;
    vector<double> bestDistancesPerGeneration;

    auto start = chrono::high_resolution_clock::now();
    time_t start_time_t = chrono::system_clock::to_time_t(start);
    string startTime = ctime(&start_time_t);

    for (int generation = 0; generation < NUM_GENERATIONS; ++generation) {
        vector<double> fitness(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            fitness[i] = 1.0 / calculateTotalDistance(cities, population[i]);
        }

        vector<vector<int>> newPopulation;
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            vector<int> parent1 = selectParent(population, fitness);
            vector<int> parent2 = selectParent(population, fitness);
            vector<int> child = crossover(parent1, parent2);
            mutate(child, MUTATION_RATE);
            newPopulation.push_back(child);
        }
        population = newPopulation;

        vector<int> bestTour = population[0];
        double bestDistance = calculateTotalDistance(cities, bestTour);
        double bestFitness = 1.0 / bestDistance;
        for (int i = 1; i < POPULATION_SIZE; ++i) {
            double distance = calculateTotalDistance(cities, population[i]);
            double fitnessValue = 1.0 / distance;
            if (fitnessValue > bestFitness) {
                bestTour = population[i];
                bestDistance = distance;
                bestFitness = fitnessValue;
            }
        }
        bestToursPerGeneration.push_back(bestTour);
        bestFitnessesPerGeneration.push_back(bestFitness);
        bestDistancesPerGeneration.push_back(bestDistance);
    }

    auto end = chrono::high_resolution_clock::now();
    time_t end_time_t = chrono::system_clock::to_time_t(end);
    string endTime = ctime(&end_time_t);
    chrono::duration<double> elapsed = end - start;

    vector<int> bestTour = bestToursPerGeneration.back();
    double bestFitness = bestFitnessesPerGeneration.back();
    double bestDistance = bestDistancesPerGeneration.back();

    saveToFile2("file2.txt", bestTour, bestFitness, bestDistance, NUM_GENERATIONS, startTime, endTime, elapsed.count());
    saveToFile3("file3.txt", bestToursPerGeneration, bestFitnessesPerGeneration, bestDistancesPerGeneration);

    cout << "Best tour: ";
    for (int city : bestTour) {
        cout << city << " ";
    }
    cout << endl;
    cout << "Total distance: " << bestDistance << endl;
    cout << "Fitness: " << bestFitness << endl;
    cout << "Elapsed time (s): " << elapsed.count() << endl;

    return 0;
}

