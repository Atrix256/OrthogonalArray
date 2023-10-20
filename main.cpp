#include <stdio.h>
#include <vector>
#include <unordered_set>
#include <random>

// If true, will try options in a random order
#define RANDOMIZE_ORDER() true

// If true, will give the same random option order each run (useful for debugging)
#define DETERMINISTIC() true

// If deterministic, will use this value as the random seed (useful for debugging)
#define DETERMINISTIC_SEED() 0x7DA771D9

// An option to try for the orthogonal array
struct Option
{
    int index = 0;
    bool used = false;
};

std::mt19937 GetRNG()
{
#if DETERMINISTIC()
    std::mt19937 rng(DETERMINISTIC_SEED());
#else
    std::random_device rd;
    auto a = rd();
    printf("Seed = 0x%X\n", a);
    std::mt19937 rng(a);
#endif
    return rng;
}

// Gauss' formula
int Sum1ToN(int N)
{
    return (N * (N + 1)) / 2;
}

// Inverted Gauss' formula
// https://blog.demofox.org/2023/08/21/inverting-gauss-formula/
int WhatNForSum(int sum)
{
    return (int)std::floor(std::sqrt(2.0f * float(sum) + 0.25) - 0.5f);
}

// greatest common divisor
int CalculateGCD(int a, int b)
{
    // if a > b, GCD(a,b) is the same as GCD(a - b, b)
    // aka Euclids algorithm
    while (1)
    {
        if (a == b)
            return a;

        if (a < b)
            std::swap(a, b);

        a -= b;
    }
}

// Lowest common multiple
int CalculateLCM(int a, int b)
{
    return a * b / CalculateGCD(a, b);
}

// Convert a column pair index into the columns it represents.
// With 3 columns, that would give Sum1ToN(3-1)=3 column pairs, in this order:
// 0,1
// 0,2
// 1,2
// The basic idea is that for N columns, there are N-1 pairs, then N-2 pairs, then N-3 pairs, all the way down to 1 pair, then none.
// This is a sum of integers from 1 to N, but in reverse.
// If we reverse the index, we can ask what N we would need to sum 1 to N to not go over the reversed index as a sum.
// that gives us our reversed first column index.
// We can then use the remainder to get the second column index.
// For more information on this code read this: https://blog.demofox.org/2023/08/21/inverting-gauss-formula/
void ColumnPairIndexToColumns(int columnPairIndex, int& firstColumnIndex, int& secondColumnIndex, int NumColumns, int NumColumnPairs)
{
    int reversedColumnPairIndex = NumColumnPairs - columnPairIndex - 1;
    int reversedGroupIndex = WhatNForSum(reversedColumnPairIndex);
    firstColumnIndex = NumColumns - reversedGroupIndex - 1;

    int reversedNextGroupStartIndex = Sum1ToN(reversedGroupIndex + 1);
    int distanceToNextGroup = reversedNextGroupStartIndex - reversedColumnPairIndex;
    secondColumnIndex = firstColumnIndex + distanceToNextGroup;

    // convert from 1 based indices to 0 based
    firstColumnIndex--;
    secondColumnIndex--;
}

void OptionIndexToOptions(int optionIndex, const std::vector<int>& levels, std::vector<int>& option)
{
    option.resize(levels.size());
    for (int i = (int)levels.size() - 1; i >= 0; --i)
    {
        int valueCount = levels[i];
        option[i] = optionIndex % valueCount;
        optionIndex /= valueCount;
    }
}

void UnapplyOption(const std::vector<int>& levels, std::vector<std::vector<int>>& columnPairValueCounts, std::vector<int>& option)
{
    for (int columnPairIndex = 0; columnPairIndex < (int)columnPairValueCounts.size(); ++columnPairIndex)
    {
        std::vector<int>& valueCounts = columnPairValueCounts[columnPairIndex];

        int firstColumn = 0;
        int secondColumn = 0;
        ColumnPairIndexToColumns(columnPairIndex, firstColumn, secondColumn, (int)levels.size(), (int)columnPairValueCounts.size());

        int value = option[firstColumn] * levels[secondColumn] + option[secondColumn];

        valueCounts[value]--;
    }
}

bool ApplyOption(const std::vector<int>& levels, std::vector<std::vector<int>>& columnPairValueCounts, std::vector<int>& option, int countMultiplier)
{
    // First see if we can apply the option
    for (int columnPairIndex = 0; columnPairIndex < (int)columnPairValueCounts.size(); ++columnPairIndex)
    {
        std::vector<int>& valueCounts = columnPairValueCounts[columnPairIndex];

        int firstColumn = 0;
        int secondColumn = 0;
        ColumnPairIndexToColumns(columnPairIndex, firstColumn, secondColumn, (int)levels.size(), (int)columnPairValueCounts.size());

        int value = option[firstColumn] * levels[secondColumn] + option[secondColumn];

        if (valueCounts[value] >= countMultiplier)
            return false;
    }

    // If we got here, we can apply the option, so do so
    for (int columnPairIndex = 0; columnPairIndex < (int)columnPairValueCounts.size(); ++columnPairIndex)
    {
        std::vector<int>& valueCounts = columnPairValueCounts[columnPairIndex];

        int firstColumn = 0;
        int secondColumn = 0;
        ColumnPairIndexToColumns(columnPairIndex, firstColumn, secondColumn, (int)levels.size(), (int)columnPairValueCounts.size());

        int value = option[firstColumn] * levels[secondColumn] + option[secondColumn];

        valueCounts[value]++;
    }

    return true;
}

int main(int argc, char** argv)
{
    int countMultiplier = 0;
    std::vector<int> levels;

    // Read CountMultiplier and Levels in.
    {
        if (argc < 3)
        {
            printf(
                "Usage: OrthogonalArrays <CountMultiplier> <Levels0> ... <LevelsN>\n"
                "\n"
                "    CountMultiplier - how many times each column value combination appears.\n"
                "                      1 is fewest tests, higher valeus for more tests.\n"
                "\n"
                "    Levels0..N      - How many values (levels) each parameter (factor) has.\n"
                "\n"
            );
            return 1;
        }

        if (sscanf_s(argv[1], "%i", &countMultiplier) != 1)
        {
            printf("Error, could not read CountMultiplier");
            return 1;
        }

        levels.resize(argc - 2);
        for (int i = 2; i < argc; ++i)
        {
            if (sscanf_s(argv[i], "%i", &levels[i - 2]) != 1)
            {
                printf("Error, could not read level%i", i - 2);
                return 1;
            }
        }

        printf("Generating experiments\n  CountMultiplier=%i\n", countMultiplier);
        for (int i = 0; i < levels.size(); ++i)
            printf("  Level%i = %i\n", i, levels[i]);
    }

    // Get the unique number of combinations needed for each pair of columns
    std::vector<int> combinationCounts;
    {
        // gather them into a set
        std::unordered_set<int> combinationCountsSet;
        for (int i = 0; i < levels.size() - 1; ++i)
            for (int j = i + 1; j < levels.size(); ++j)
                combinationCountsSet.insert(levels[i] * levels[j]);

        // gather them into a vector
        for (int count : combinationCountsSet)
            combinationCounts.push_back(count);
    }

    // Get the least common multiple of column combinations, so we know
    // what the minimum number of experiments is.
    int LCM = combinationCounts[0];
    for (size_t i = 1; i < combinationCounts.size(); ++i)
        LCM = CalculateLCM(LCM, combinationCounts[i]);

    // calculate the number of experiments we are going to do
    int experimentCount = LCM * countMultiplier;

    // Calculate the total number of experiments we would need for full factorial
    int fullFactorialCount = 1;
    for (int level : levels)
        fullFactorialCount *= level;

    // report the experiment counts
    printf("\n%i experiments instead of %i with full factorial\n", experimentCount, fullFactorialCount);

    // Make a data structure to keep track of our solution possibilities.
    // We want to keep information about each column pair.
    // The information we want to keep is the count of each combination of column values.
    // There are Sum1ToN(columns-1) column pairs.
    // There are LevelsA * LevelsB possible values for the column pair AB.
    std::vector<std::vector<int>> columnPairValueCounts;
    columnPairValueCounts.resize(Sum1ToN((int)levels.size()-1));
    for (int columnPairIndex = 0; columnPairIndex < (int)columnPairValueCounts.size(); ++columnPairIndex)
    {
        int firstColumnIndex = 0;
        int secondColumnIndex = 0;
        ColumnPairIndexToColumns(columnPairIndex, firstColumnIndex, secondColumnIndex, (int)levels.size(), (int)columnPairValueCounts.size());
        columnPairValueCounts[columnPairIndex].resize(levels[firstColumnIndex] * levels[secondColumnIndex], 0);
    }

    // Make a list of all options
    std::vector<Option> options(fullFactorialCount);
    for (int optionIndex = 0; optionIndex < fullFactorialCount; ++optionIndex)
        options[optionIndex].index = optionIndex;

    // shuffle the options if we should
    #if RANDOMIZE_ORDER()
    std::shuffle(options.begin(), options.end(), GetRNG());
    #endif

    // Fill the first columnPairValueCounts in order, and let the other slots fill them in, in the process.
    // As an example, if column 0 and column 1 were 2 level factors (binary parameters),
    // we'd fill 00, then 01, then 10, then 11.  We need "countMultiplier" of each.
    // We end when stack is "experimentCount" in size
    std::vector<int> solutionStack;
    {
        std::vector<int> scanIndices(experimentCount, 0);

        int firstColumnValueCount = levels[0];
        int secondColumnValueCount = levels[1];
        int firstTwoColumnCombinationCount = firstColumnValueCount * secondColumnValueCount;
        std::vector<int> option;
        while (solutionStack.size() < experimentCount)
        {
            // get the values we want for the first and second column
            int stackSize = (int)solutionStack.size();
            int columnPairValueIndex = stackSize / countMultiplier;
            int desiredSecondColumnValue = columnPairValueIndex % secondColumnValueCount;
            int desiredFirstColumnValue = columnPairValueIndex / secondColumnValueCount;

            int scanIndex = scanIndices[stackSize];
            if (scanIndex >= options.size())
            {
                printf("Could not find a solution (1)!\n\n");
                return 2;
            }

            while (scanIndex < options.size())
            {
                if (options[scanIndex].used)
                {
                    scanIndex++;
                    continue;
                }

                OptionIndexToOptions(options[scanIndex].index, levels, option);
                if (option[0] != desiredFirstColumnValue || option[1] != desiredSecondColumnValue)
                {
                    scanIndex++;
                    continue;
                }

                if (ApplyOption(levels, columnPairValueCounts, option, countMultiplier))
                {
                    options[scanIndex].used = true;
                    solutionStack.push_back(scanIndex);
                    scanIndices[stackSize] = scanIndex + 1;

                    if (solutionStack.size() < options.size() && ((stackSize + 1) % countMultiplier != 0))
                        scanIndices[stackSize + 1] = scanIndex + 1;

                    break;
                }

                scanIndex++;
            }

            // If we hit this, we hit a dead end, so need to back track
            if (scanIndex >= options.size())
            {
                if (solutionStack.size() == 0)
                {
                    printf("Could not find a solution (2)!\n\n");
                    return 3;
                }

                int lastSolutionScanIndex = *solutionStack.rbegin();
                solutionStack.pop_back();
                OptionIndexToOptions(options[lastSolutionScanIndex].index, levels, option);
                UnapplyOption(levels, columnPairValueCounts, option);
                options[lastSolutionScanIndex].used = false;
                scanIndices[solutionStack.size()] = lastSolutionScanIndex + 1;
                scanIndices[solutionStack.size() + 1] = 0;
            }
        }

        // sort the solution stack
        std::sort(solutionStack.begin(), solutionStack.end(),
            [&](int scanIndexA, int scanIndexB)
            {
                return options[scanIndexA].index < options[scanIndexB].index;
            }
        );

        // Print out the results!
        for (int scanIndex : solutionStack)
        {
            int solution = options[scanIndex].index;
            printf("[%.3i]", solution);
            OptionIndexToOptions(solution, levels, option);
            bool first = true;
            for (int value : option)
            {
                printf("%s%i", first ? " " : ", ", value);
                first = false;
            }
            printf("\n");
        }
    }

    return 0;
}
/*
TODO: when randomized, it can fail to find a solution. investigate!
 
NOTES:
* just strength 2 arrays. you can modify the code to have strength 3 or higher arrays
* could run probably faster as algorithm x with dancing links
* the rabbit hole goes far deeper

*/