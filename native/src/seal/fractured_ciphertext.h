#pragma once

#include "fractured_polynomial.h"
#include <seal/seal.h>
#include <seal/util/polyarithsmallmod.h>

namespace seal::fractures {
    struct CiphertextFracture {
        std::vector<PolynomialFracture> poly_fracs;
        std::uint64_t index;
        std::vector<seal::Modulus> coeff_modulus;

        const CiphertextFracture &operator*=(const PolynomialFracture &poly);

        // changes self? no, creates a copy and resturns it.
        CiphertextFracture operator*(const PolynomialFracture &poly) {
            CiphertextFracture res(*this);
            return res *= poly;
        }


    };

/**
 * This class is responsible for fracturing ciphertext into multiple parts.
 * Each fractured part should be able to be multiplied with a fractured plainetext.
 * Using fractured ciphertexts and plaintexts, we can perform multiplication in parallel on multiple machines.
 */
    class CiphertextShredder {
    public:
        explicit CiphertextShredder(seal::Ciphertext &ctx,
                                    std::uint64_t modulus_size,
                                    std::uint64_t num_coefficients,
                                    const std::vector<seal::Modulus> &coeff_modulus,
                                    std::uint64_t num_fractures) noexcept;

        CiphertextFracture get_fracture(std::uint64_t index);

        CiphertextFracture &operator[](std::uint64_t index) {
            return ctx_parts[index];
        }

        //todo: func that turns this shredder back to full ciphertext.
        void into_ciphertext(seal::Ciphertext &ciphertext);

    private:
        std::uint64_t num_coefficients;
        std::vector<CiphertextFracture> ctx_parts;
    };
}