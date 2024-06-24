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

        // assert correct padding.
        ASSERT_EQ(ptx.dyn_array().size(), 4 * 3 * N); // 4x is from the padding, 3x is from the number of moduli.
        ASSERT_TRUE(ptx.coeff_count() != cpy.coeff_count());

        seal::util::RNSIter ptx_itr(ptx.data(), N * 4);
        seal::util::RNSIter cpy_itr(cpy.data(), N);
        SEAL_ITERATE(
            seal::util::iter(ptx_itr, cpy_itr), context.first_context_data()->parms().coeff_modulus().size(),
            [&](auto I) {
                auto i1(std::get<0>(I));
                auto i2(std::get<1>(I));

                for (std::uint64_t i = 0; i < N; ++i)
                {
                    ASSERT_TRUE(*i1 == *i2);
                    i1++;
                    i2++;
                }

                // Now assert the inbetween are zeros.
                for (std::uint64_t i = N; i < 4 * N; ++i)
                {
                    ASSERT_TRUE(*i1 == 0);
                    i1++;
                }
            });
    }
} // namespace sealtest
