// Test for fractured operations by JonathanWeiss.

#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractured_poly.h"
#include "seal/keygenerator.h"
#include "seal/modulus.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <string>
#include "gtest/gtest.h"

using namespace seal;
using namespace std;
namespace sealtest
{

    TEST(FracturedOps, Poly)
    {
        // enc params
        std::uint64_t N = 4096;
        int logt = 20;
        seal::EncryptionParameters enc_params(seal::scheme_type::bgv);
        enc_params.set_poly_modulus_degree(N);
        enc_params.set_coeff_modulus(seal::CoeffModulus::BFVDefault(N));
        enc_params.set_plain_modulus(seal::PlainModulus::Batching(N, logt + 1));

        // context

        seal::Plaintext p("1x^2 + 2x^1 + 3x^0");
        seal::SEALContext context_data(enc_params, true);
        seal::fractures::Essence e(context_data, enc_params);

        seal::fractures::Poly p_frac(p, 3, 3, 2);

        std::cout << "hello friend" << std::endl;
    }
} // namespace sealtest