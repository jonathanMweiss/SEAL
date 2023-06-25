#pragma once

#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "fractured_rns.h"
#include "matrix.h"

namespace seal::fractures
{

    struct CiphertextFracture;
    /**
     * A class that represents a part of a polynomial that had been split.
     */
    struct PolynomialFracture
    {
        PolynomialFracture(std::uint64_t index, std::uint64_t _coeff_count, std::uint64_t modulus_size)
            : fracture_index(index), coeff_count(_coeff_count), rns_coefficients(coeff_count, modulus_size)
        {}

        inline seal::util::CoeffIter rns_poly_iter(std::uint64_t rns_num)
        {
            return { &rns_coefficients(0, rns_num) };
        }

        inline seal::util::ConstCoeffIter const_rns_poly_iter(std::uint64_t rns_num) const
        {
            return { &rns_coefficients(0, rns_num) };
        }
        CiphertextFracture operator*(const CiphertextFracture &ctxf) const;

        std::uint64_t fracture_index;
        std::uint64_t coeff_count;
        seal::util::matrix<std::uint64_t> rns_coefficients;
    };

    /**
     * Fractures a polynomial into multiple PolynomialFracture.
     * these polynomials are used to perform multiplication between a ciphertext and a plaintext in a distributed
     * manner.
     */
    class Polynomial
    {
    public:
        explicit Polynomial(const seal::util::ConstRNSIter &p, Essence e, std::uint64_t _num_fractures) noexcept;

        explicit Polynomial(const seal::Plaintext &p, const Essence &e, std::uint64_t _num_fractures) noexcept
            : Polynomial(seal::util::ConstRNSIter(p.data(), e.coeff_count), e, _num_fractures)
        {}

        const PolynomialFracture &get_fracture(std::uint64_t index);

    private:
        static PolynomialFracture compute_fracture(
            const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
            uint64_t num_fractures, uint64_t index);

        std::vector<PolynomialFracture> fractures;
        std::uint64_t num_fractures;
        const Essence essence;
    };
} // namespace seal::fractures