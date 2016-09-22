#include "Spin.h"

namespace yap {

//-------------------------
const SpinProjectionVector projections(unsigned two_j)
{
    SpinProjectionVector spv;
    spv.reserve(two_j + 1);
    for (int two_m = -two_j; two_m <= (int)two_j; two_m += 2)
        spv.push_back(two_m);
    return spv;
}

//-------------------------
const std::vector<SpinProjectionVector> projections(const SpinVector& two_J)
{
    // initialize vector of spin projections to -two_j
    std::vector<int> two_M;
    two_M.reserve(two_J.size());
    std::transform(two_J.begin(), two_J.end(), std::back_inserter(two_M),
    [](const SpinVector::value_type & two_j) {return -two_j;});

    std::vector<SpinProjectionVector> SPV;
    // fill SPV with "odometer"-style looping
    while (two_M.back() <= (int)two_J.back()) {
        SPV.push_back(two_M);
        two_M[0] += 2;
        for (size_t i = 0; (i < two_M.size() - 1) and (two_M[i] > (int)two_J[i]); ++i) {
            two_M[i] = -two_J[i];
            two_M[i + 1] += 2;
        }
    }
    return SPV;
}

//-------------------------
const SpinVector triangle(unsigned two_A, unsigned two_B)
{
    unsigned two_s0 = abs(two_A - two_B);
    unsigned two_s1 = two_A + two_B;
    SpinVector two_S;
    two_S.reserve(two_s1 - two_s0 + 1);
    for (unsigned two_s = two_s0; two_s <= two_s1; two_s += 2)
        two_S.push_back(two_s);
    return two_S;
}
    
}
