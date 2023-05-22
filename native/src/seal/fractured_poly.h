#pragma once

#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "matrix.h"
#include "rns_fracture.h"

namespace seal::fractures
{

    class Poly
    {
        matrix<std::uint64_t> poly_data;
        std::uint64_t num_fractures;
        Essence essence;

        // gets a polynomial iterator, that starts at the correct position in the matrix.
        seal::util::CoeffIter rns_poly_iter(std::uint64_t rns_num);
        seal::util::ConstCoeffIter const_rns_poly_iter(std::uint64_t rns_num);

    public:
        // Creates a plaintext with the given number of coefficients and modulus size.
        static Poly from_plaintext(
            const seal::Evaluator &ev, const seal::fractures::Essence &e, std::uint16_t num_fractures,
            const seal::Plaintext &ptx);

        explicit Poly(seal::util::ConstRNSIter rns_iter, const Essence &e, std::uint64_t num_fractures) noexcept;


        static void memcpy(util::CoeffIter &cpy_to, util::ConstCoeffIter &cpy_from, std::uint64_t coeffcount);

    private:
        // Expects p to be in NTT form.
        explicit Poly(seal::Plaintext &p, const Essence &e, std::uint64_t num_fractures) noexcept
            : Poly(seal::util::RNSIter(p.data(), e.coeff_count), e, num_fractures){};
    };

} // namespace seal::fractures
