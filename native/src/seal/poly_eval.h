#pragma once
#include "fractured_ciphertext.h"
#include "fractured_polynomial.h"

namespace seal::fractures
{
    typedef PolynomialFracture EvaluatedPoint;
    typedef CiphertextFracture EvaluatedCipherPoint;
    class PolynomialEvaluator
    {
    public:
        // TODO: instead of receiving a value and breaking it, receive some random values already in the correct mod.
        explicit PolynomialEvaluator(seal::fractures::Essence e);
        EvaluatedPoint evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &value) const;
        EvaluatedPoint evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &&value) const;

        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, std::vector<std::uint64_t> &&value) const;
        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, std::vector<std::uint64_t> &value) const;

    private:
        Essence essence;
        EvaluatedPoint evalute(seal::util::ConstRNSIter rns_iter, const std::vector<std::uint64_t> &value) const;

        EvaluatedPoint evalute(
            seal::util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus,
            const std::vector<std::uint64_t> &value) const;

        static void validate_value_to_evaluate(
            const std::vector<std::uint64_t> &value, const std::vector<Modulus> &coeff_modulus);
    };
} // namespace seal::fractures
