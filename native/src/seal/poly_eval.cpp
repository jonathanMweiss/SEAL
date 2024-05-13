#include "poly_eval.h"
#include "util/uintarith.h"

namespace seal::fractures
{
    // TODO: When evaluating a ptx, we might need to transfer it first to the CTX modulus domain, and ONLY then
    //  evaluate. need to understand how to do it correctly.

    //    EvaluatedPoint seal::fractures::PolynomialEvaluator::evaluate(
    //        seal::Plaintext &p, std::vector<std::uint64_t> &value) const
    //    {
    //        if (p.is_ntt_form())
    //        {
    //            throw std::invalid_argument("Plaintext is already in NTT form");
    //        }
    //
    //        auto ctx_data = context.get_context_data(p.parms_id());
    //        auto &parms = ctx_data->parms();
    //        std::vector<seal::Modulus> plain_modulus{ parms.plain_modulus() };
    //
    //        validate_value_to_evaluate(value, plain_modulus);
    //
    //        size_t coeff_count = parms.poly_modulus_degree();
    //        size_t poly_size = p.dyn_array().size();
    //
    //        if (poly_size < coeff_count)
    //        {
    //            throw std::invalid_argument("Plaintext is too small. put in RNS representation.");
    //        }
    //
    //        return evalute(seal::util::ConstRNSIter(p.data(), essence.coeff_count), plain_modulus, value);
    //    }

    EvaluatedCipherPoint PolynomialEvaluator::evaluate(const Ciphertext &ctx, std::vector<std::uint64_t> &value) const
    {
        if (ctx.is_ntt_form())
        {
            throw std::invalid_argument("Ciphertext is already in NTT form");
        }

        auto ctx_data = context.get_context_data(ctx.parms_id());
        auto &parms = ctx_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();

        validate_value_to_evaluate(value, coeff_modulus);

        // creating an empty one to reserve space.
        auto result = seal::fractures::EvaluatedCipherPoint::Empty(ctx.size(), 0, 0, coeff_modulus);

        std::uint64_t i = 0;
        SEAL_ITERATE(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter iter) {
            result.poly_fracs[i] = evaluate_singe_RNS_polynomial(iter, coeff_modulus, value);
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
        seal::util::ConstRNSIter rns_iter, const std::vector<seal::Modulus> &modulus,
        const std::vector<std::uint64_t> &value) const
    {
        auto poly_degree = context.first_context_data()->parms().poly_modulus_degree();
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
} // namespace seal::fractures
