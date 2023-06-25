

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
        std::uint64_t coeff_count;

        explicit Essence(seal::SEALContext _ctx, seal::EncryptionParameters params)
            : ctx(std::move(_ctx)), parms(std::move(params)), coeff_count(parms.poly_modulus_degree()){};
    };
} // namespace seal::fractures