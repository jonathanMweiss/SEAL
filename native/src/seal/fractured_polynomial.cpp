#include "fractured_polynomial.h"

namespace seal::fractures {

    PolynomialFracture
    PolynomialShredder::compute_fracture(const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size,
                                         uint64_t num_coeffs, uint64_t num_fractures,
                                         uint64_t index) {

        PolynomialFracture fracture(index, num_coeffs / num_fractures, modulus_size);

        auto rns_num = -1;
        std::for_each_n(seal::util::iter(rns_iter), modulus_size, [&](auto coef_iter) {
            rns_num++;

            coef_iter += index * num_coeffs / num_fractures;
            auto frac_coef_iter = fracture.rns_poly_iter(rns_num);
            for (uint64_t j = 0; j < num_coeffs / num_fractures; ++j) {
                *frac_coef_iter = *coef_iter;
                frac_coef_iter++;
                coef_iter++;
            }
        });


        return fracture;
    }

    PolynomialShredder::PolynomialShredder(const seal::util::ConstRNSIter &p, std::uint64_t modulus_size,
                                           std::uint64_t num_coefficients, std::uint64_t num_fractures) noexcept {
        fractures.reserve(num_fractures);

        for (std::uint64_t i = 0; i < num_fractures; ++i) {
            fractures.push_back(compute_fracture(p, modulus_size, num_coefficients, num_fractures, i));
        }
    }

    const PolynomialFracture &PolynomialShredder::get_fracture(std::uint64_t index) {
        if (index >= fractures.size()) {
            throw std::invalid_argument("index out of bounds");
        }

        return fractures[index];
    }
}