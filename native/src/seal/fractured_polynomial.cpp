#include "fractured_polynomial.h"

namespace seal::fractures
{

    PolynomialFracture Polynomial::compute_fracture(
        const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
        uint64_t num_fractures, uint64_t index)
    {
        PolynomialFracture fracture(index, num_coeffs / num_fractures, modulus_size);

        auto rns_num = -1;
        std::for_each_n(seal::util::iter(rns_iter), modulus_size, [&](auto coef_iter) {
            rns_num++;

            coef_iter += index * num_coeffs / num_fractures;
            auto frac_coef_iter = fracture.rns_poly_iter(rns_num);
            for (uint64_t j = 0; j < num_coeffs / num_fractures; ++j)
            {
                *frac_coef_iter = *coef_iter;
                frac_coef_iter++;
                coef_iter++;
            }
        });

        return fracture;
    }

    Polynomial::Polynomial(const seal::util::ConstRNSIter &p, Essence e, std::uint64_t num_fractures) noexcept
        : poly_data(e.coeff_count, e.coeff_modulus_size), num_fractures(num_fractures), essence(std::move(e))
    {
        fractures.reserve(num_fractures);

        for (std::uint64_t i = 0; i < num_fractures; ++i)
        {
            fractures.push_back(compute_fracture(p, essence.coeff_modulus_size, essence.coeff_count, num_fractures, i));
        }
    }

    const PolynomialFracture &Polynomial::get_fracture(std::uint64_t index)
    {
        if (index >= fractures.size())
        {
            throw std::invalid_argument("index out of bounds");
        }

        return fractures[index];
    }
} // namespace seal::fractures