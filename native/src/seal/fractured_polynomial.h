#pragma once

#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "matrix.h"
#include "rns_fracture.h"

namespace seal::fractures
{
    /**
     * A class that represents a part of a polynomial that had been split.
     */
    struct PolynomialFracture
    {
        PolynomialFracture(std::uint64_t index, std::uint64_t coeff_count, std::uint64_t modulus_size)
            : fracture_index(index), coeff_count(coeff_count), rns_coefficients(coeff_count, modulus_size)
        {}

        seal::util::CoeffIter rns_poly_iter(std::uint64_t rns_num)
        {
            return { &rns_coefficients(0, rns_num) };
        }

        seal::util::ConstCoeffIter const_rns_poly_iter(std::uint64_t rns_num) const
        {
            return { &rns_coefficients(0, rns_num) };
        }

        std::uint64_t fracture_index;
        std::uint64_t coeff_count;
        matrix<std::uint64_t> rns_coefficients;
    };

    /**
     * Fractures a polynomial into multiple PolynomialFracture.
     * these polynomials are used to perform multiplication between a ciphertext and a plaintext in a distributed
     * manner.
     */
    class Polynomial
    {
    public:
        explicit Polynomial(const seal::util::ConstRNSIter &p, Essence e, std::uint64_t num_fractures) noexcept;

        explicit Polynomial(seal::Plaintext &p, const Essence &e, std::uint64_t num_fractures) noexcept
            : Polynomial(seal::util::ConstRNSIter(p.data(), e.coeff_count), e, num_fractures)
        {}

        const PolynomialFracture &get_fracture(std::uint64_t index);

    private:
        static PolynomialFracture compute_fracture(
            const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
            uint64_t num_fractures, uint64_t index);

        std::vector<PolynomialFracture> fractures;
        matrix<std::uint64_t> poly_data;
        std::uint64_t num_fractures;
        const Essence essence;
    };
} // namespace seal::fractures