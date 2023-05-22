// Test for fractured operations by JonathanWeiss.

#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractured_ciphertext.h"
#include "seal/fractured_poly.h"
#include "seal/fractured_polynomial.h"
#include "seal/keygenerator.h"
#include "seal/modulus.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <string>
#include <utility>
#include "gtest/gtest.h"

using namespace seal;
using namespace std;
namespace sealtest
{
    struct SetupObjs
    {
        seal::EncryptionParameters enc_params;
        seal::SEALContext context;
        seal::fractures::Essence essence;
        seal::Evaluator evaluator;
        seal::KeyGenerator keygen;
        seal::SecretKey secret_key;
        seal::Encryptor encryptor;
        seal::Decryptor decryptor;

        static SetupObjs New(std::uint64_t N = 4096, int logt = 20)
        {
            seal::EncryptionParameters enc(seal::scheme_type::bgv);

            enc.set_poly_modulus_degree(N);
            enc.set_coeff_modulus(seal::CoeffModulus::BFVDefault(N));
            enc.set_plain_modulus(seal::PlainModulus::Batching(N, logt + 1));

            return SetupObjs(enc);
        }

        explicit SetupObjs(seal::EncryptionParameters encryption_params)
            : enc_params(std::move(encryption_params)), context(enc_params, true), essence(context, enc_params),
              evaluator(context), keygen(context), secret_key(keygen.secret_key()), encryptor(context, secret_key),
              decryptor(context, secret_key){};

        seal::Plaintext random_plaintext() const
        {
            seal::Plaintext p(enc_params.poly_modulus_degree());

            // avoiding encoder usage - to prevent unwanted transformation to the ptx underlying elements
            std::vector<std::uint64_t> v(enc_params.poly_modulus_degree(), 0);

            seal::Blake2xbPRNGFactory factory;
            std::array<std::uint64_t, seal::prng_seed_uint64_count> seed{ 1, 2, 3, 4, 5, 6, 7, 8 };
            auto gen = factory.create(seed);

            auto mod = enc_params.plain_modulus().value();
            std::generate(v.begin(), v.end(), [gen = std::move(gen), &mod]() { return gen->generate() % mod; });

            for (std::uint64_t i = 0; i < enc_params.poly_modulus_degree(); ++i)
            {
                p[i] = v[i];
            }
            return p;
        }

        seal::Ciphertext random_ciphertext() const
        {
            auto ptx = random_plaintext();
            seal::Ciphertext ctx;
            encryptor.encrypt_symmetric(ptx, ctx);
            return ctx;
        }
    };

    TEST(FracturedOps, MultiplyAndReconstructSingleShard)
    {
        auto all = SetupObjs::New();
        std::uint64_t num_fractures = 256;

        seal::Ciphertext encrypted_ntt = all.random_ciphertext();
        auto ptx = all.random_plaintext();
        all.evaluator.transform_to_ntt_inplace(ptx, encrypted_ntt.parms_id());

        seal::Ciphertext encrypted_ntt_cpy(encrypted_ntt);
        seal::Plaintext ptx_cpy(ptx);

        seal::fractures::CiphertextShredder cshredder(
            encrypted_ntt_cpy, all.essence.coeff_modulus_size, all.essence.coeff_count, all.essence.coeff_modulus,
            num_fractures);
        seal::fractures::PolynomialShredder pshredder(
            ptx_cpy, all.essence.coeff_modulus_size, all.essence.coeff_count, num_fractures);

        // perform fractured multiplication:
        for (int i = 0; i < num_fractures; ++i)
        {
            cshredder[i] *= pshredder.get_fracture(i);
        } // regular multiplication and then we'll compare them:
        all.evaluator.multiply_plain_inplace(encrypted_ntt, ptx);

        // reconstruct into a valid ctx in ntt form.
        auto reconstructed = all.random_ciphertext();
        cshredder.into_ciphertext(reconstructed);

        assert(!encrypted_ntt.is_transparent());
        all.evaluator.sub_inplace(encrypted_ntt, reconstructed);
        assert(encrypted_ntt.is_transparent());

        seal::Plaintext ptx_res;
        all.decryptor.decrypt(reconstructed, ptx_res);
        std::cout << "reconstructed: " << ptx_res.to_string() << std::endl;
    }
} // namespace sealtest