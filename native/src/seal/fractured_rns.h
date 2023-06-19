

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
    public:
        //            const seal::SEALContext::ContextData context_data;
        seal::SEALContext ctx;
        seal::EncryptionParameters parms;
        std::vector<seal::Modulus> coeff_modulus;
        std::uint64_t coeff_count;
        std::uint64_t coeff_modulus_size;

        explicit Essence(seal::SEALContext _ctx, seal::EncryptionParameters params)
            : ctx(std::move(_ctx)), parms(std::move(params)), coeff_modulus(parms.coeff_modulus()),
              coeff_count(parms.poly_modulus_degree()), coeff_modulus_size(coeff_modulus.size()){};
    };
    class RnsFracture
    {
    public:
        RnsFracture(
            const seal::util::matrix<std::uint64_t> &rns_coefficients, std::uint64_t num_fractures,
            std::uint64_t frac_index);
    };
} // namespace seal::fractures