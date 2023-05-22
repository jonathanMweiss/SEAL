

#include "rns_fracture.h"

// i have coef count/ num fractures
// i have num-cols as the
//seal::fractures::RnsFracture::RnsFracture(
//    const matrix<std::uint64_t> &rns_coefficients, std::uint64_t num_fractures, std::uint64_t frac_index,
//    std::uint64_t modulus_size)
//{
//    PolynomialFracture fracture(index, num_coeffs / num_fractures, modulus_size);
//
//    auto rns_num = -1;
//    std::for_each_n(seal::util::iter(rns_iter), modulus_size, [&](auto coef_iter) {
//        rns_num++;
//
//        coef_iter += index * num_coeffs / num_fractures;
//        auto frac_coef_iter = fracture.rns_poly_iter(rns_num);
//        for (uint64_t j = 0; j < num_coeffs / num_fractures; ++j)
//        {
//            *frac_coef_iter = *coef_iter;
//            frac_coef_iter++;
//            coef_iter++;
//        }
//    });
//};
