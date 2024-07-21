#include "poly_eval.h"
#include "util/uintarith.h"

namespace seal::fractures
{

    EvaluatedPoint seal::fractures::PolynomialEvaluator::evaluate(
        seal::Plaintext &p, const std::vector<std::uint64_t> &value) const
    {
        if (p.is_ntt_form())
        {
            throw std::invalid_argument("Plaintext is already in NTT form, cannot evaluate it.");
        }

        auto ctx_data = context.first_context_data();
        // otherwise, we need to uplift the ptx first, and i didn't want to implement it/ extract it out of the
        // evaluator.
        if (!ctx_data->qualifiers().using_fast_plain_lift)
        {
            throw std::invalid_argument(
                "haven't implemented polynomial evaluation without 'using_fast_plain_lift' set to true!.");
        }

        ev.plain_to_coeff_space(p, context.first_parms_id());

        auto &parms = ctx_data->parms();
        return evaluate_singe_RNS_polynomial(
            seal::util::ConstRNSIter(p.data(), parms.poly_modulus_degree()), parms.coeff_modulus(), value,
            context.first_parms_id());
    }

    EvaluatedCipherPoint PolynomialEvaluator::evaluate(
        const seal::Ciphertext &ctx, const std::vector<std::uint64_t> &value) const
    {
        if (ctx.is_ntt_form())
        {
            throw std::invalid_argument("Ciphertext is already in NTT form, cannot evaluate it.");
        }

        auto ctx_data = context.get_context_data(ctx.parms_id());
        auto &parms = ctx_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();

        validate_value_to_evaluate(value, coeff_modulus);

        // creating an empty one to reserve space.
        auto result = seal::fractures::EvaluatedCipherPoint::Empty(ctx.size(), 0, 0, coeff_modulus);

        std::uint64_t i = 0;
        SEAL_ITERATE(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter iter) {
            result.poly_fracs[i] = evaluate_singe_RNS_polynomial(iter, coeff_modulus, value, ctx.parms_id());
            i++;
        });

        return result;
    }
    void PolynomialEvaluator::validate_value_to_evaluate(
        const std::vector<std::uint64_t> &value, const std::vector<Modulus> &coeff_modulus) const
    {
        if (value.size() != coeff_modulus.size())
        {
            throw std::invalid_argument("value size must be equal to the number of moduli.");
        }
        for (size_t i = 0; i < coeff_modulus.size(); ++i)
        {
            if (value[i] >= coeff_modulus[i].value())
            {
                throw std::invalid_argument("value is not in the correct range.");
            }
        }
    }

    EvaluatedPoint PolynomialEvaluator::evaluate_singe_RNS_polynomial(
        util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus, const std::vector<std::uint64_t> &value,
        const parms_id_type p_id) const
    {
        // TODO: must use relevant context, for instance - padded ctx have different contexr!
        auto poly_degree = context.get_context_data(p_id)->parms().poly_modulus_degree();
        std::vector<std::uint64_t> result_vec(modulus.size());

        uint64_t i = 0;
        SEAL_ITERATE(seal::util::iter(rns_iter), modulus.size(), [&](auto coef_iter) {
            auto mod = modulus[i];

            uint64_t x = 1;
            uint64_t sum = 0;
            for (uint64_t j = 0; j < poly_degree; ++j)
            {
                // sum += coef*(value^i)
                sum = util::multiply_add_uint_mod(*coef_iter, x, sum, mod);

                // x = value^(i+1) Advance: after using the power.
                x = util::multiply_uint_mod(x, value[i], mod);

                coef_iter++;
            }

            result_vec[i] = sum;

            // advance to the next modulus on the RNS iterator.
            i++;
        });

        EvaluatedPoint point = EvaluatedPoint(0, 1, modulus.size());
        point.rns_coefficients.data = std::move(result_vec);
        return point;
    }
    void PolynomialEvaluator::multiply_scalar(
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
    std::vector<std::vector<std::uint64_t>> PolynomialEvaluator::generate_possible_polyeval_ring_homomorphism_values(
        const SEALContext &cntx, seal::parms_id_type pid)
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
    std::vector<std::uint64_t> PolynomialEvaluator::generate_first_polyeval_ring_homomorphism_value(
        const SEALContext &cntx, seal::parms_id_type pid)
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
    std::uint64_t PolynomialEvaluator::gen_root_of_unity(const Modulus &m, std::uint64_t n)
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
        return a;
        return seal::util::exponentiate_uint_mod(a, (m.value() - 1 / n), m);
    }
} // namespace seal::fractures
