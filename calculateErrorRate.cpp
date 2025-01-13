#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

// experimental distribution, value at index i stands for probability that a random sample takes value i.
double EXPERIMENTAL[25] = {
    0.0483871,
    0.02419355,
    0.03225806,
    0.05241935,
    0.0766129,
    0.08870968,
    0.06854839,
    0.05241935,
    0.06048387,
    0.0766129,
    0.08064516,
    0.05241935,
    0.03629032,
    0.04435484,
    0.03225806,
    0.02419355,
    0.02016129,
    0.05241935,
    0.01209677,
    0.02016129,
    0.01612903,
    0.00806452,
    0.00806452,
    0.00806452,
    0.00403226
    };

int getFactorial(int number){
    int result = 1;
    for (int i=0; i<number; i++){
        result *= i + 1;
    }
    return result;
}

// Returns the number of unique permutations of the input List. It is calculated by dividing the factorial of the number of elements by the product of factorials of the numbers of equal samples. For {0, 0, 0, 1, 1, 1, 1, 1, 1, 1} this would be 10! / (3! * 7!).
int getCombinatoricFactor(std::vector<int> samples){
    // since samples is sorted, check if next element is equal to previous and increase equals count
    int equals = 1;
    int denominator = 1;
    for (int i=0; i<samples.size()-1; i++){
    if (samples[i] == samples[i+1]){
        equals ++;
    }
    else{
        // if a different value is found, update the denominator and reset the count of equal elements
        denominator *= getFactorial(equals);
        equals = 1;
    }
    }
    denominator *= getFactorial(equals);

    return getFactorial(samples.size()) / denominator;
}

// Returns the probability to draw the given set of samples from the given distribution.
double getSampleProbability(std::vector<int> samples, double distribution[25]){

    double probability = 1.;

    for (int i=0; i<samples.size(); i++){
        // distribution vector does only contain non-zero values, if sampled value is higher, it is assumed to be zero
        if (samples[i] >= 25){
            return 0.;
        }
        probability *= distribution[samples[i]];
    }

    return probability;
}

// updates samples by 1 such that it will always be sorted in ascending order. Returns true while updating and false if input has reached the maximum e.g. all values are equal to max
bool updateSamples(std::vector<int> &samples, int max=24){
    // return false if all samples are max on input
    if (samples.back() == max && samples.front() == max){
        return false;
    }
    // start iterating from back of array
    int i = samples.size() - 1;
    while (i >= 0){
        // if sample is lower than max, increase its value by one
        if (samples[i] < max){
            samples[i]++;
            break;
        }
        // else i will be index of largest sample that is < max
        else{
            i--;
        }
    }
    // since list has to be sorted, set all following samples to that value
    int j = i;
    while (i < samples.size()){
        samples[i] = samples[j];
        i++;
    }
    return true;
}

// returns probability of samples beeing drawn from the less likely distribution, weighted by the number of unique permutations of the sample list
double pWrongGuess(std::vector<int> samples, double distribution1[25], double distribution2[25]){
    double p1 = getSampleProbability(samples, distribution1);
    double p2 = getSampleProbability(samples, distribution2);
    int combinatoricFactor = getCombinatoricFactor(samples);
    return combinatoricFactor * ((p1 < p2) ? p1 : p2) / 2;
}

double getErrorRate(int numberSamples, double dist1[25], double dist2[25]){
    std::chrono::_V2::system_clock::time_point time = std::chrono::high_resolution_clock::now();
    long startTime = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    std::cout << "running\n";

    std::vector<int> samples = {};
    for (int i=0; i<numberSamples; i++){
        samples.push_back(0);
    }

    // start with initial samples, loop will change samples before calculating probability
    double errorRate = pWrongGuess(samples, dist1, dist2);
    long termsAdded = 1;

    while(updateSamples(samples)){
        errorRate += pWrongGuess(samples, dist1, dist2);
        termsAdded ++;
    }

    time = std::chrono::high_resolution_clock::now();
    long calcTime = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count() - startTime;
    std::cout << "Summed up " << termsAdded << " terms in " << (float) calcTime / 1000 << "s.\n";
    return errorRate;
}

int main() {
    int equalUpperLimit = 18;
    int numberSamples = 10;

    // static size is much faster so distributions are limited to size 25
    double equal[25];
    for (int i=0; i<25; i++){
        equal[i] = (i < equalUpperLimit) ? (double) 1 / equalUpperLimit : 0;
    }

    double errorRate = getErrorRate(numberSamples, EXPERIMENTAL, equal);
    std::cout << "error rate: " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << errorRate << "\n";

    return 0;
}