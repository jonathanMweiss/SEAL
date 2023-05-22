

#pragma once
#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "matrix.h"

namespace seal::fractures
{

    // Represents the key parameters to fracture a polynomial.
    struct Essence
    {
        //            const seal::SEALContext::ContextData context_data;
        seal::EncryptionParameters parms;
        std::vector<seal::Modulus> coeff_modulus;
        std::uint64_t coeff_count;
        std::uint64_t coeff_modulus_size;

        explicit Essence(const seal::SEALContext &context_data, seal::EncryptionParameters params)
            : parms(std::move(params)), coeff_modulus(parms.coeff_modulus()), coeff_count(parms.poly_modulus_degree()),
              coeff_modulus_size(coeff_modulus.size()){};
    };
    class RnsFracture
    {
    public:
        RnsFracture(
            const matrix<std::uint64_t> &rns_coefficients, std::uint64_t num_fractures, std::uint64_t frac_index);
    };
} // namespace seal::fractures