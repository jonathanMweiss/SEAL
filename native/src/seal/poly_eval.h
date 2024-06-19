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
        PolynomialEvaluator(const SEALContext &ctx) : context(ctx), ev(ctx){};

        EvaluatedPoint evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &&value) const
        {
            return evaluate(p, value);
        }

        EvaluatedPoint evaluate(const seal::Plaintext &p, const std::vector<std::uint64_t> &value) const
        {
            seal::Plaintext tmp_cpy(p);
            return evaluate(tmp_cpy, value);
        }

        EvaluatedPoint evaluate(seal::Plaintext &p, const std::vector<std::uint64_t> &value) const;

        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, std::vector<std::uint64_t> &&value) const
        {
            return evaluate(ctx, value);
        }

        EvaluatedCipherPoint evaluate(const seal::Ciphertext &ctx, const std::vector<std::uint64_t> &value) const;

        /**
         * Generate a root of unity for a given modulus and n.
         * @param m The modulus.
         * @param n The n-th root of unity.
         * @return The root of unity.
         */
        static std::uint64_t gen_root_of_unity(const seal::Modulus &m, std::uint64_t n);

        /**
         * generate possible evaluation values for polynomials.
         *
         * To evaluate polynomials over Z_q[X]/(X^n + 1),
         * one should use values which can be used to create an evaluation map $f_a(P(X))$ for a polynomial P(X),
         * such that $f_a$ is a ring homomorphism. That is, $f_a(P(X) + Q(X)) = f_a(P(X)) + f_a(Q(X))$ and
         * $f_a(P(X) * Q(X)) = f_a(P(X)) * f_a(Q(X))$ for every $P(X), Q(X) \in Z_q[X]/(X^n + 1)$.
         *
         * $f_a$ is a ring homomorphism if and only if $a$ is a root of $X^n+1$.
         *
         * It is hard to find such elements. Luckily, we know that $(X^2n-1)=(X^n+1)(X^n-1)$ (thus if $X^n+1$ has roots
         * they must be roots for $X^2n-1$, too). In addition, we know that the $2n-th$ roots of unity are the roots for
         * $X^2n-1$. (since every $(w_i)^2n = 1$) and since all of them can't be the roots of $X^n-1$ (i.e., $X^n-1$ has
         * at most $n$ solutions)., the other $n$ elements must be roots for $X^n+1$.
         * @return
         */
        static std::vector<std::uint64_t> generate_first_polyeval_ring_homomorphism_value(
            const seal::SEALContext &cntx, seal::parms_id_type pid);

        static void multiply_scalar(
            const std::vector<std::uint64_t> &a, std::vector<std::uint64_t> &dest,
            const std::vector<seal::Modulus> &modulus);

        static std::vector<std::vector<std::uint64_t>> generate_possible_polyeval_ring_homomorphism_values(
            const seal::SEALContext &cntx, seal::parms_id_type pid);

    private:
        EvaluatedPoint evaluate_singe_RNS_polynomial(
            seal::util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus,
            const std::vector<std::uint64_t> &value) const;

        SEALContext context;
        seal::Evaluator ev;
        void validate_value_to_evaluate(
            const std::vector<std::uint64_t> &value, const std::vector<Modulus> &coeff_modulus) const;
    };
} // namespace seal::fractures
