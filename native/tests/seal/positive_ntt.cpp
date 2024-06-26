// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/keygenerator.h"
#include "seal/modulus.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <string>
#include "gtest/gtest.h"
#include "helpers.h"

using namespace seal;
using namespace std;

namespace sealtest
{
    /**
     * @param N the poly_degree of non padded polynomials.
     * @param expansion_ratio
     */
    void validate_rns_poly_adding(
        seal::util::PtrIter<std::uint64_t *> padded, seal::util::PtrIter<std::uint64_t *> regular, std::uint64_t N,
        std::uint64_t expansion_ratio = 4)
    {
        auto i1(padded);
        auto i2(regular);

        for (std::uint64_t i = 0; i < N; ++i)
        {
            ASSERT_EQ(*i1, *i2);
            i1++;
            i2++;
        }

        // Now assert the inbetween are zeros.
        for (std::uint64_t i = N; i < expansion_ratio * N; ++i)
        {
            if (*i1 != 0)
            {
                std::cout << "i1: " << *i1 << std::endl;
            }
            ASSERT_EQ(*i1, 0);
            i1++;
        }
    }
    TEST(EvaluatorTest, PostiveWrappedNTTPlaintextPadding)
    {
        std::uint64_t N = 128;
        // The common parameters: the plaintext and the polynomial moduli
        Modulus plain_modulus(65);

        // The parameters and the context of the higher level
        EncryptionParameters parms(scheme_type::bgv);
        parms.set_poly_modulus_degree(128);
        parms.set_plain_modulus(plain_modulus);
        parms.set_coeff_modulus(CoeffModulus::Create(N, { 30, 30, 30, 30 }));

        SEALContext context(parms, false, sec_level_type::none, 2);
        Evaluator evaluator(context);

        Plaintext ptx = random_plain(parms);
        evaluator.plain_to_coeff_space(ptx, context.first_parms_id());
        Plaintext cpy = ptx;

        evaluator.zero_pad(ptx, context.first_parms_id());

        std::uint64_t padding_expansation_ration = 4;
        // assert correct padding.
        ASSERT_EQ(
            ptx.dyn_array().size(),
            padding_expansation_ration * 3 * N); // 4x is from the padding, 3x is from the number of moduli.
        ASSERT_TRUE(ptx.coeff_count() != cpy.coeff_count());

        seal::util::RNSIter ptx_itr(ptx.data(), N * 4);
        seal::util::RNSIter cpy_itr(cpy.data(), N);
        SEAL_ITERATE(
            seal::util::iter(ptx_itr, cpy_itr), context.first_context_data()->parms().coeff_modulus().size(),
            [&](tuple<seal::util::PtrIter<std::uint64_t *>, seal::util::PtrIter<std::uint64_t *>> I) {
                validate_rns_poly_adding(std::get<0>(I), std::get<1>(I), N, padding_expansation_ration);
            });
    }

    TEST(EvaluatorTest, PostiveWrappedNTTCiphertextPadding)
    {
        std::uint64_t N = 128;
        // The common parameters: the plaintext and the polynomial moduli
        Modulus plain_modulus(65);

        // The parameters and the context of the higher level
        EncryptionParameters parms(scheme_type::bgv);
        parms.set_poly_modulus_degree(128);
        parms.set_plain_modulus(plain_modulus);
        parms.set_coeff_modulus(CoeffModulus::Create(N, { 30, 30, 30, 30 }));

        SEALContext context(parms, false, sec_level_type::none, 1);
        Evaluator evaluator(context);

        Plaintext ptx = random_plain(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());
        auto parms_id = context.first_parms_id();

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);

        seal::Ciphertext cpy;
        cpy = encrypted;

        evaluator.zero_pad(encrypted);

        // create two polynomial iterators:
        seal::util::PolyIter reg_iter(cpy);
        seal::util::PolyIter padded_iter(encrypted);

        std::vector<std::uint64_t> pd;

        SEAL_ITERATE(
            seal::util::iter(padded_iter, reg_iter), cpy.size(),
            [&](tuple<seal::util::RNSIter, seal::util::RNSIter> I) {
                SEAL_ITERATE(seal::util::iter(std::get<0>(I), std::get<1>(I)), cpy.coeff_modulus_size(), [&](auto J) {
                    // now check they are equal for the frst N elements., and the rest are zeros.
                    validate_rns_poly_adding(std::get<0>(J), std::get<1>(J), N, 2);
                });
            });
    }
} // namespace sealtest
