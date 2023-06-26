#include "poly_eval.h"
#include "util/uintarith.h"

namespace seal::fractures
{
    PolynomialEvaluator::PolynomialEvaluator(seal::fractures::Essence e) : essence(std::move(e))
    {}
    EvaluatedPoint seal::fractures::PolynomialEvaluator::evaluate(
        seal::Plaintext &p, std::vector<std::uint64_t> &&value) const
    {
        auto tmp = value;
        return evaluate(p, tmp);
    }

    EvaluatedPoint seal::fractures::PolynomialEvaluator::evaluate(
        seal::Plaintext &p, std::vector<std::uint64_t> &value) const
    {
        if (p.is_ntt_form())
        {
            throw std::invalid_argument("Plaintext is already in NTT form");
        }
        auto ctx_data = essence.ctx.get_context_data(essence.ctx.first_parms_id());
        auto &parms = ctx_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();

        if (p.dyn_array().size() < coeff_count * coeff_modulus_size)
        {
            throw std::invalid_argument("Plaintext is too small. put in RNS representation.");
        }

        auto rns_iter = seal::util::ConstRNSIter(p.data(), essence.coeff_count);
        return evalute(rns_iter, value);
    }

    // TODO something is wrong here...
    EvaluatedPoint seal::fractures::PolynomialEvaluator::evalute(
        seal::util::ConstRNSIter rns_iter, const std::vector<std::uint64_t> &value) const
    {
        auto ctx_data = essence.ctx.get_context_data(essence.ctx.first_parms_id());
        auto &parms = ctx_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();

        std::vector<std::uint64_t> result_vec(coeff_modulus_size);

        uint64_t i = 0;

        rns_iter++;
        std::for_each_n(seal::util::iter(rns_iter), coeff_modulus_size, [&](auto coef_iter) {
            auto mod = coeff_modulus[i];
            // need to take the current powah and start going through it.
            uint64_t x = 1;
            uint64_t sum = 0;
            for (uint64_t j = 0; j < essence.coeff_count; ++j)
            {
                // sum += coef*(value^i)
                sum = util::multiply_add_uint_mod(*coef_iter, x, sum, mod);

                // x = value^i+1 Advance: after using the power.
                x = util::multiply_uint_mod(x, value[i], mod);
                coef_iter++;
            }

            // advance the iteration.
            result_vec[i] = sum;
            i++;
        });

        EvaluatedPoint point = EvaluatedPoint(0, 1, coeff_modulus_size);
        point.rns_coefficients.data = std::move(result_vec);
        return point;
    }

    EvaluatedCipherPoint PolynomialEvaluator::evaluate(const Ciphertext &ctx, std::vector<std::uint64_t> &value) const
    {
        if (ctx.is_ntt_form())
        {
            throw std::invalid_argument("Ciphertext is already in NTT form");
        }

        auto ctx_data = essence.ctx.get_context_data(essence.ctx.first_parms_id());
        auto &parms = ctx_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto result = seal::fractures::EvaluatedCipherPoint::Empty(ctx.size(), 0, 0, coeff_modulus);

        std::uint64_t i = 0;
        std::for_each_n(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter iter) {
            result.poly_fracs[i++] = evalute(iter, value);
        });

        return result;
    }

    EvaluatedCipherPoint PolynomialEvaluator::evaluate(const Ciphertext &ctx, std::vector<std::uint64_t> &&value) const
    {
        auto tmp = value;
        return evaluate(ctx, tmp);
    }
} // namespace seal::fractures
