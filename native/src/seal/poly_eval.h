#pragma once
#include "fractured_ciphertext.h"
#include "fractured_polynomial.h"

namespace seal::fractures
{
    /**
     * The evaluated point is consisting of number of scalars (since it is in CRT form) that represent the polynomial.
     */
    typedef PolynomialFracture EvaluatedPoint;

    /**
     * The inner structure of these evaluated point is more complex. in essence it is a HyperPolynomial with
     * deg 0 polynomial as its coefficients (tuple of scalars). This is important, since multiplying two such points
     * follow regular polynomial multiplication, and not point-2-point multiplication.
     */
    typedef CiphertextFracture EvaluatedCipherPoint;
    class PolynomialEvaluator
    {
    public:
        PolynomialEvaluator(const SEALContext &ctx) : context(ctx){};

        EvaluatedPoint evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &&value) const
        {
            return evaluate(p, value);
        }
        EvaluatedPoint evaluate(seal::Plaintext &p, const std::vector<std::uint64_t> &value) const;

        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, std::vector<std::uint64_t> &&value) const
        {
            return evaluate(ctx, value);
        }

        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, const std::vector<std::uint64_t> &value) const;

    private:
        EvaluatedPoint evaluate_singe_RNS_polynomial(
            seal::util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus,
            const std::vector<std::uint64_t> &value) const;

        SEALContext context;
        void validate_value_to_evaluate(
            const std::vector<std::uint64_t> &value, const std::vector<Modulus> &coeff_modulus) const;
    };
} // namespace seal::fractures
