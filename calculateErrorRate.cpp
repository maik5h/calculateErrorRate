// This code is used to calculate the lowest possible error rate in the following scenario:
// There are two discrete distributions of integer numbers. A set of n samples is drawn from
// either of the two, with probability 0.5. Based on these samples, how accurate can
// be predicted which one of the distributions was sampled from?
// 
// One of the distributions is an approximate distribution of feature numbers on experimental
// images, the other one was chosen to be equal with varying upper limits.
//
// It was decided to avoid vectors to have all data in static memory for much faster reading.
// Multithreading was not necessary.

#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <fstream>

// This defines the number of double values that are preserved in memory to hold the values of the discrete PDFs
// used for calculations. 27 is large enough to fit all examined distributions.
constexpr int reservedSize = 27;

// A class to store a finite discrete probability distribution. Every integer number in the range
// [0, size-1] is assigned a probability which is stored in the probabilities attribute.
template<int size>
class discreteProbabilityDistribution {
    double probabilities[size]; // Value at index i gives the probability for a random sample to have the value i.

    public:
    discreteProbabilityDistribution(double inProbabilities[size]) {
        for (int i=0; i<size; i++) {
            probabilities[i] = inProbabilities[i];
        }
    }

    // Returns the probability to draw the given samples from this distribution.
    // The total probability is the product of the probabilities for each sample.
    double getProbabilityForSamples(int *samples, int sampleSize) {

        double probability = 1.;

        for (int i=0; i<sampleSize; i++){
            probability *= probabilities[samples[i]];
        }

        return probability;
    }
};

// Experimental distribution, value at index i stands for probability that a random sample takes value i.
double EXPERIMENTAL[27] = {
    0.048, 0.024, 0.032, 0.052, 0.076, 0.088, 0.068, 0.052, 0.060,
    0.076, 0.080, 0.052, 0.036, 0.044, 0.032, 0.024, 0.020, 0.052,
    0.012, 0.020, 0.016, 0.008, 0.008, 0.008, 0.004, 0.000, 0.008,
};

// Returns the current time in ms since Unix epoch.
long getTime() {
    std::chrono::_V2::system_clock::time_point time = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
}

int getFactorial(int number){
    int result = 1;
    for (int i=0; i<number; i++){
        result *= i + 1;
    }
    return result;
}

// Appends the results of the calculation to a csv file.
//  - path: csv file path
//  - uL: equal distribution upper limit
//  - sN: number of samples
//  - eR: minimal error rate
//  - t: computation time
void writeResults(std::string path, int uL, int sN, double eR, double t) {
    std::ofstream logFile(path, std::ios_base::app);
    if (logFile.is_open()) {
        logFile << uL << ", " << sN << ", " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << eR << ", " << t << ";" << std::endl;
    }
}

// Returns the number of unique permutations of the input array, which has to be sorted in ascending order.
// It is calculated by dividing the factorial of the number of elements by the product of factorials of the numbers of equal samples. For {0, 0, 0, 1, 1, 1, 1, 1, 1, 1} this would be 10! / (3! * 7!).
int getCombinatoricFactor(int *samples, int numberSamples){
    int equals = 1; // number of elements that are equal to the current element
    int denominator = 1; // product of factorials of the number of which every possible element occurs in samples

    for (int i=0; i<numberSamples; i++){
        // since samples is sorted, check if next element is equal to previous and increase equals count
        if (samples[i] == samples[i+1]){
            equals ++;
        }
        // if a different value is found, update the denominator and reset the count of equal elements
        else{
            denominator *= getFactorial(equals);
            equals = 1;
        }
    }
    denominator *= getFactorial(equals);

    return getFactorial(numberSamples) / denominator;
}

// Updates the last sample in samples that has not reached max, such that samples is in ascending order. 
// Returns true while updating and false if input has reached the maximum e.g. all values are equal to max.
bool updateSamples(int *samples, int numberSamples, int max=26){
    // Return false if all samples are max on input.
    if (samples[0] == max && samples[numberSamples-1] == max){
        return false;
    }
    
    // Start iterating from back of array
    int i = numberSamples - 1;
    while (i >= 0){
        // If sample is lower than max, increase its value by one
        if (samples[i] < max){
            samples[i]++;
            break;
        }
        // If sample is max, procede with previous sample.
        else{
            i--;
        }
    }
    // Since ascending order is required, the value of last updated sample represents an upper bound to all
    // samples at a higher index. Reset all samples with higher index to this value.
    int lowerBound = samples[i];
    while (i < numberSamples){
        samples[i] = lowerBound;
        i++;
    }
    return true;
}

// Returns probability of samples beeing drawn from the less likely distribution, weighted by the number of unique permutations of the sample list.
double getPWrongGuess(int *samples, int numberSamples, discreteProbabilityDistribution<reservedSize> dist1, discreteProbabilityDistribution<reservedSize> dist2){
    double p1 = dist1.getProbabilityForSamples(samples, numberSamples);
    double p2 = dist2.getProbabilityForSamples(samples, numberSamples);
    int combinatoricFactor = getCombinatoricFactor(samples, numberSamples);
    return combinatoricFactor * ((p1 < p2) ? p1 : p2) / 2;
}

// Returns the minimum possible error rate when differentiating between the two distributions based on numberSamples samples.
double getErrorRate(discreteProbabilityDistribution<reservedSize> dist1, discreteProbabilityDistribution<reservedSize> dist2, int numberSamples){

    // Set up samples with every sample beeing zero
    int samples[numberSamples];
    for (int i=0; i<numberSamples; i++){
        samples[i] = 0;
    }

    // Start with initial samples
    double errorRate = getPWrongGuess(samples, numberSamples, dist1, dist2);

    // Iterate over all ordered combinations of samples
    while(updateSamples(samples, numberSamples)){
        errorRate += getPWrongGuess(samples, numberSamples, dist1, dist2);
    }

    return errorRate;
}

int main() {
    int sampleNumbers[2] = {1, 10}; // Number of samples to consider.
    int equalUpperLimits[6] = {7, 9, 11, 13, 15, 17}; // Defines the upper limit of the equal distribution (including this value).
    std::string outFile = "./results.csv";

    for (auto sampleNumber : sampleNumbers){
        for (auto equalUpperLimit : equalUpperLimits){
            // Set up an array for the equal distribution. Is filled up to index equalUpperLimit with equal values.
            double equal[reservedSize];
            for (int i=0; i<25; i++){
                equal[i] = (i <= equalUpperLimit) ? (double) 1 / (equalUpperLimit + 1) : 0.;
            }

            // Define two distributions:
            // First one is an approximate experimental distribution,
            // second one is an equal distribution over the interval [0, equalUpperLimit].
            discreteProbabilityDistribution<reservedSize> dist1(EXPERIMENTAL);
            discreteProbabilityDistribution<reservedSize> dist2(equal);

            // Keep track of computation time
            long startTime = getTime();
            std::cout << "running\n";
            
            double errorRate = getErrorRate(dist1, dist2, sampleNumber);
        
            // Log and write results
            long calcTime = getTime() - startTime;
            std::cout << "upper limit: " << equalUpperLimit << ", number of samples: " << sampleNumber << "\n";
            std::cout << "Calculation time: " << (float) calcTime / 1000 << "s.\n";
            std::cout << "error rate: " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << errorRate << "\n\n";

            writeResults(outFile, equalUpperLimit, sampleNumber, errorRate, (float) calcTime / 1000);
        }
    }
    return 0;
}