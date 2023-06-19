#pragma once

#include <seal/seal.h>
#include <seal/util/polyarithsmallmod.h>
#include "fractured_polynomial.h"

namespace seal::fractures
{
    struct CiphertextFracture
    {
        static CiphertextFracture Empty(
            std::uint64_t num_polys, std::uint64_t index, std::uint64_t num_coefs, std::vector<seal::Modulus> coef_mod);

        static CiphertextFracture Empty(const CiphertextFracture &ctxf);

        std::vector<PolynomialFracture> poly_fracs;
        std::uint64_t index;
        std::vector<seal::Modulus> coeff_modulus;

        const CiphertextFracture &operator*=(const PolynomialFracture &poly);
        const CiphertextFracture &operator+=(const CiphertextFracture &ctxf);
        const CiphertextFracture &operator*=(const CiphertextFracture &y);

        // TODO: verify.
        CiphertextFracture operator*(const PolynomialFracture &poly) const;
        CiphertextFracture operator*(const CiphertextFracture &y) const;
        CiphertextFracture operator+(const CiphertextFracture &ctxf);

    private:
        inline void add(const CiphertextFracture &ctxf, uint64_t i, uint64_t rns_num);
    };

    /**
     * This class is responsible for fracturing ciphertext into multiple parts.
     * Each fractured part should be able to be multiplied with a fractured plainetext.
     * Using fractured ciphertexts and plaintexts, we can perform multiplication in parallel on multiple machines.
     */
    class CiphertextShredder
    {
    public:
        CiphertextShredder(Essence e, std::uint64_t num_fractures) noexcept;

        explicit CiphertextShredder(const seal::Ciphertext &ctx, Essence e, std::uint64_t num_fractures) noexcept;

        CiphertextFracture &operator[](std::uint64_t index)
        {
            return ctx_parts[index];
        }

        // todo: func that turns this shredder back to full ciphertext.
        seal::Ciphertext into_ciphertext() const;

        uint64_t num_fractures();

        void set_fracture(uint64_t i, CiphertextFracture fracture);

    private:
        seal::fractures::Essence essence;
        std::vector<CiphertextFracture> ctx_parts;
    };
} // namespace seal::fractures