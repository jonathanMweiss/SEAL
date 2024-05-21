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

        /**
         * Generate a root of unity for a given modulus and n.
         * @param m The modulus.
         * @param n The n-th root of unity.
         * @return The root of unity.
         */
        static std::uint64_t gen_root_of_unity(const seal::Modulus &m, std::uint64_t n)
        {
            // $a$ is primitive element (root) if each element that isn't $0$ can be expressed as
            // $a^k$ for some $k$.
            std::uint64_t a;
            seal::util::try_minimal_primitive_root(n, m, a);

            // a^(q-1)/n is an n-th root of unity.
            // proof: we know that $a^(q-1)$ is 1 (fermat's little theorem),
            // as a result, $a^k\ne 1$ for all $k < q-1$. otherwise there is a cycle that doesn't go through every
            // element other than 0 in the field, in contradiction to $a$ being a generator/ primitive element.
            // Thus, $w=a^(q-1)/n$ is an n-th root of unity, because w^n = a^(q-1) = 1 and w^k != 1 for all k < n.
            return seal::util::exponentiate_uint_mod(a, (m.value() - 1 / n), m);
        }

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
            const seal::SEALContext &cntx, seal::parms_id_type pid)
        {
            auto parms = cntx.get_context_data(pid)->parms();

            std::vector<std::uint64_t> r;
            for (std::uint64_t i = 0; i < parms.coeff_modulus().size(); ++i)
            {
                auto m = parms.coeff_modulus()[i];
                r.emplace_back(gen_root_of_unity(m, parms.poly_modulus_degree() * 2));
            }
            return r;
        }

        static void multiply_scalar(
            const std::vector<std::uint64_t> &a, std::vector<std::uint64_t> &dest,
            const std::vector<seal::Modulus> &modulus)
        {
            if (a.size() != dest.size())
            {
                throw std::invalid_argument("multiply_scalar: a and dest must have the same size");
            }
            if (a.size() > modulus.size())
            {
                throw std::invalid_argument("multiply_scalar: a and modulus must have the same size");
            }

            for (std::uint64_t i = 0; i < a.size(); ++i)
            {
                dest[i] = util::multiply_uint_mod(a[i], dest[i], modulus[i]);
            }
        }

        static std::vector<std::vector<std::uint64_t>> generate_possible_polyeval_ring_homomorphism_values(
            const seal::SEALContext &cntx, seal::parms_id_type pid)
        {
            auto parms = cntx.get_context_data(pid)->parms();
            auto r = generate_first_polyeval_ring_homomorphism_value(cntx, pid);

            std::vector<std::vector<std::uint64_t>> roots;
            std::vector<std::uint64_t> cur(r);
            for (std::uint64_t i = 0; i < parms.poly_modulus_degree() * 2; ++i)
            {
                if (0 == (i & 1))
                {
                    roots.push_back({ cur });
                }
                multiply_scalar(r, cur, parms.coeff_modulus());
            }

            return roots;
        }

    private:
        EvaluatedPoint evaluate_singe_RNS_polynomial(
            seal::util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus,
            const std::vector<std::uint64_t> &value) const;

        SEALContext context;
        void validate_value_to_evaluate(
            const std::vector<std::uint64_t> &value, const std::vector<Modulus> &coeff_modulus) const;
    };
} // namespace seal::fractures
