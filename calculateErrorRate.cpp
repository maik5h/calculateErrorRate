#include <cstdlib>
#include <iostream>
#include <ctime>
#include <chrono>
#include <map>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
using std::cout, std::map;

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

int getCombinatoricFactor(std::vector<int> samples){
    //Returns the number of unique permutations of the input List. It is calculated by dividing the factorial of the number of elements by the product of factorials over the numbers of equal samples. For {0, 0, 0, 1, 1, 1, 1, 1, 1, 1} this would be 10! / (3! * 7!).
    // since samples is sorted, check if next element is equal to previous and increase equals count
    int equals = 1;
    int denominator = 1;
    for (int i=0; i<samples.size()-1; i++){
    if (samples[i] == samples[i+1]){
        equals ++;
    }
    else{
        denominator *= getFactorial(equals);
        equals = 1;
    }
    }
    denominator *= getFactorial(equals);

    return getFactorial(samples.size()) / denominator;
}

double getSampleProbability(std::vector<int> samples, double distribution[25]){
    // Returns the probability to draw samples from distribution.
    double probability = 1.;

    for (int i=0; i<samples.size(); i++){
        // distribution vector does only contain non-zero values, if sampled value is higher, the probability must be zero
        if (samples[i] >= 25){
            return 0.;
        }
        probability *= distribution[samples[i]];
    }

    return probability;
}

bool updateSamples(std::vector<int> &samples, int max=24){
    // return false if all samples are max on input
    if (samples.back() == max && samples.front() == max){
        return false;
    }
    int i = samples.size() - 1;
    while (i >= 0){
        // if sample is lower than max, increase its value by one
        if (samples[i] < max){
            samples[i]++;
            break;
        }
        // i will be index of largest sample that is < max
        else{
            i--;
        }
    }
    // set all following samples to that value
    int j = i;
    while (i < samples.size()){
        samples[i] = samples[j];
        i++;
    }
    return true;
}

double pWrongGuess(std::vector<int> samples, double distribution1[25], double distribution2[25]){
    double p1 = getSampleProbability(samples, distribution1);
    double p2 = getSampleProbability(samples, distribution2);
    int combinatoricFactor = getCombinatoricFactor(samples);
    return combinatoricFactor*(((p1 < p2) ? p1 : p2) / 2);
}

double getErrorRate(int numberSamples, double dist1[25], double dist2[25]){
    std::chrono::_V2::system_clock::time_point time = std::chrono::high_resolution_clock::now();
    long startTime = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    cout << "running\n";

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
    cout << "Summed up " << termsAdded << " terms in " << (float) calcTime / 1000 << "s.\n";
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
    cout << "error rate: " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << errorRate << "\n";

    return 0;
}