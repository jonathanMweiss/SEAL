#include "poly_eval.h"
#include "util/uintarith.h"

namespace seal::fractures
{
    PolynomialEvaluator::PolynomialEvaluator(seal::fractures::Essence e) : essence(std::move(e))
    {}
    EvaluatedPoint seal::fractures::PolynomialEvaluator::Evaluate(
        seal::Plaintext &p, std::vector<std::uint64_t> &&value)
    {
        auto tmp = value;
        return Evaluate(p, tmp);
    }

    EvaluatedPoint seal::fractures::PolynomialEvaluator::Evaluate(seal::Plaintext &p, std::vector<std::uint64_t> &value)
    {
        // TODO: Lift polynomial into RNS space and to be in the correct moduli of the coeff.
        //   need to learn how to do it. Otherwise i cannot work with it as a fracture/ mult it with the polys of a ctx.
        if (p.is_ntt_form())
        {
            throw std::invalid_argument("Plaintext is already in NTT form");
        }

        std::vector<std::uint64_t> result_vec(essence.coeff_modulus.size());

        auto rns_iter = seal::util::ConstRNSIter(p.data(), essence.coeff_count);
        std::uint64_t i = 0;

        rns_iter++;
        std::for_each_n(seal::util::iter(rns_iter), essence.coeff_modulus.size() - 1, [&](auto coef_iter) {
            auto mod = essence.coeff_modulus[i];
            // need to take the current powah and start going through it.
            std::uint64_t x = 1;
            std::uint64_t sum = 0;
            for (std::uint64_t j = 0; j < essence.coeff_count; ++j)
            {
                // sum += coef*(value^i)
                sum = seal::util::multiply_add_uint_mod(*coef_iter, x, sum, mod);

                // x = value^i+1 Advance: after using the power.
                x = seal::util::multiply_uint_mod(x, value[i], mod);
                coef_iter++;
            }

            // advance the iteration.
            result_vec[i] = sum;
            i++;
        });

        EvaluatedPoint point = seal::fractures::EvaluatedPoint(0, 1, essence.coeff_modulus_size);
        point.rns_coefficients.data = std::move(result_vec);
        return point;
    }
} // namespace seal::fractures
