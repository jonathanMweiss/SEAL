#pragma once
#include "fractured_polynomial.h"

namespace seal::fractures
{
    typedef PolynomialFracture EvaluatedPoint;

    class PolynomialEvaluator
    {
    public:
        // TODO: instead of receiving a value and breaking it, receive some random values already in the correct mod.
        explicit PolynomialEvaluator(seal::fractures::Essence e);
        EvaluatedPoint Evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &value);
        EvaluatedPoint Evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &&value);

    private:
        Essence essence;
    };
} // namespace seal::fractures
