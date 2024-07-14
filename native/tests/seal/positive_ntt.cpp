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
#include <unordered_set>
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
    void validate_rns_poly_padding(
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

    void assert_eq_ciphers(const Ciphertext &e1, const Ciphertext &e2)
    {
        auto e1_vc = dynarray_to_vector(e1.dyn_array());
        auto e2_vc = dynarray_to_vector(e2.dyn_array());
        ASSERT_EQ(e1_vc.size(), e2_vc.size());
        for (std::uint64_t i = 0; i < e1_vc.size(); ++i)
        {
            ASSERT_EQ(e1_vc[i], e2_vc[i]);
        }
    }

    SEALContext make_unsecure_context(int poly_deg, int nummults)
    {
        uint64_t N = poly_deg;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        auto tmp = vector<int>{ 54, 41, 42 };
        auto o = CoeffModulus::Create(N, tmp);
        parms.set_coeff_modulus(o);
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));
        return SEALContext(parms, false, sec_level_type::none, nummults);
    }

    SEALContext make_128deg_context(int nummults = 2)
    {
        std::uint64_t N = 128;
        // The common parameters: the plaintext and the polynomial moduli
        Modulus plain_modulus(65);

        // The parameters and the context of the higher level
        EncryptionParameters parms(scheme_type::bgv);
        parms.set_poly_modulus_degree(128);
        parms.set_plain_modulus(plain_modulus);
        parms.set_coeff_modulus(CoeffModulus::Create(N, { 30, 30, 30, 30 }));

        return { parms, false, sec_level_type::none, nummults };
    }

    TEST(EvaluatorTest, PostiveWrappedNTTPlaintextPadding)
    {
        std::uint64_t N = 128;
        auto context = make_128deg_context();
        Evaluator evaluator(context);

        Plaintext ptx = random_plain(context.first_context_data()->parms());
        evaluator.plain_to_coeff_space(ptx, context.first_parms_id());
        Plaintext cpy = ptx;

        evaluator.zero_pad(ptx, context.first_parms_id());

        std::uint64_t padding_expansation_ration = 4;

        // 4x is from the padding, 3x is from the number of moduli.
        ASSERT_EQ(ptx.dyn_array().size(), padding_expansation_ration * 3 * N);
        ASSERT_TRUE(ptx.coeff_count() != cpy.coeff_count());

        seal::util::RNSIter ptx_itr(ptx.data(), N * 4);
        seal::util::RNSIter cpy_itr(cpy.data(), N);
        SEAL_ITERATE(
            seal::util::iter(ptx_itr, cpy_itr), context.first_context_data()->parms().coeff_modulus().size(),
            [&](tuple<seal::util::PtrIter<std::uint64_t *>, seal::util::PtrIter<std::uint64_t *>> I) {
                validate_rns_poly_padding(std::get<0>(I), std::get<1>(I), N, padding_expansation_ration);
            });
    }

    TEST(EvaluatorTest, PostiveWrappedNTTCiphertextPadding)
    {
        std::uint64_t N = 128;
        auto context = make_128deg_context();
        auto parms = context.first_context_data()->parms();

        Evaluator evaluator(context);

        Plaintext ptx = random_plain(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);
        evaluator.transform_from_ntt_inplace(encrypted);

        seal::Ciphertext cpy;
        cpy = encrypted;

        evaluator.zero_pad(encrypted);

        // create two polynomial iterators:
        seal::util::PolyIter reg_iter(cpy);
        seal::util::PolyIter padded_iter(encrypted);

        SEAL_ITERATE(
            seal::util::iter(padded_iter, reg_iter), cpy.size(),
            [&](tuple<seal::util::RNSIter, seal::util::RNSIter> I) {
                SEAL_ITERATE(seal::util::iter(std::get<0>(I), std::get<1>(I)), cpy.coeff_modulus_size(), [&](auto J) {
                    // now check they are equal for the frst N elements., and the rest are zeros.
                    validate_rns_poly_padding(std::get<0>(J), std::get<1>(J), N, 2);
                });
            });
    }
    TEST(EvaluatorTest, PostiveWrappedNTTCiphertextPadding2)
    {
        std::uint64_t N = 1024;
        auto context = make_unsecure_context(N, 1);
        auto parms = context.first_context_data()->parms();

        Evaluator evaluator(context);

        Plaintext ptx = random_plain(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);
        evaluator.transform_from_ntt_inplace(encrypted);

        seal::Ciphertext cpy;
        cpy = encrypted;

        evaluator.zero_pad(encrypted);

        // create two polynomial iterators:
        seal::util::PolyIter reg_iter(cpy);
        seal::util::PolyIter padded_iter(encrypted);

        auto vc = dynarray_to_vector(encrypted.dyn_array());

        SEAL_ITERATE(
            seal::util::iter(padded_iter, reg_iter), cpy.size(),
            [&](tuple<seal::util::RNSIter, seal::util::RNSIter> I) {
                SEAL_ITERATE(seal::util::iter(std::get<0>(I), std::get<1>(I)), cpy.coeff_modulus_size(), [&](auto J) {
                    // now check they are equal for the frst N elements., and the rest are zeros.
                    validate_rns_poly_padding(std::get<0>(J), std::get<1>(J), N, 1);
                });
            });
    }

    TEST(EvaluatorTest, PositiveNttConstantMonomyialMultCTx)
    {
        auto context = make_unsecure_context(8192, 1);
        auto parms = context.first_context_data()->parms();
        Evaluator evaluator(context);

        Plaintext ptx("2"); //(parms);

        // creating encryption:
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        Ciphertext ctx(context);
        encryptor.encrypt(ptx, ctx);
        evaluator.transform_from_ntt_inplace(ctx);

        // multiplying them:
        ptx = evaluator.transform_to_positive_ntt(ptx);
        ctx = evaluator.transform_to_positive_ntt(ctx);

        seal::Ciphertext dst;
        evaluator.multiply_plain(ctx, ptx, dst);

        // Testsing for correct values after multiplying. (ptx after ntt = [2,2,2,2,2.....,2]
        auto positive_wrapped_params = context.get_context_data(context.positive_wrapped_parms_id())->parms();
        // Create poly_iter
        seal::util::PolyIter s_itr(ctx);
        seal::util::PolyIter d_itr(dst);
        SEAL_ITERATE(
            seal::util::iter(s_itr, d_itr), ctx.size(), [&](tuple<seal::util::RNSIter, seal::util::RNSIter> I) {
                std::uint64_t mod_index = 0;
                SEAL_ITERATE(
                    seal::util::iter(std::get<0>(I), std::get<1>(I)), positive_wrapped_params.coeff_modulus().size(),
                    [&](auto coefItersJ) {
                        auto mod = parms.coeff_modulus().at(mod_index++);
                        auto s(std::get<0>(coefItersJ));
                        auto double_s(std::get<1>(coefItersJ));

                        for (std::uint64_t i = 0; i < positive_wrapped_params.poly_modulus_degree(); ++i)
                        {
                            auto expected_value = seal::util::multiply_uint_mod(*s, 2, mod);
                            ASSERT_EQ(expected_value, *double_s);
                        }
                    });
            });
    }
    TEST(EvaluatorTest, PositiveINTTCTX)
    {
        GTEST_SKIP();
    }

    TEST(EvaluatorTest, postivie_ntt_constant_mult)
    {
        std::uint64_t N = 8192;
        seal::EncryptionParameters parms(seal::scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        //        enc.set_coeff_modulus(seal::CoeffModulus::BFVDefault(N, sec_level_type::tc192));

        auto tmp = std::vector<int>{ 54, 41, 42 };
        auto o = seal::CoeffModulus::Create(N, tmp);
        parms.set_coeff_modulus(o);
        parms.set_plain_modulus(seal::PlainModulus::Batching(N, 16 + 1));

        SEALContext context(parms, false, sec_level_type::none, 1);
        Evaluator evaluator(context);

        Plaintext ptx("2"); //(parms);

        // creating encryption:
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        Ciphertext ctx(context);
        encryptor.encrypt(ptx, ctx);
        evaluator.transform_from_ntt_inplace(ctx);

        // multiplying them:
        ptx = evaluator.transform_to_positive_ntt(ptx);
        ctx = evaluator.transform_to_positive_ntt(ctx);

        seal::Ciphertext dst;
        evaluator.multiply_plain(ctx, ptx, dst);

        // Testsing for correct values after multiplying. (ptx after ntt = [2,2,2,2,2.....,2]
        auto positive_wrapped_params = context.get_context_data(context.positive_wrapped_parms_id())->parms();
        // Create poly_iter
        seal::util::PolyIter s_itr(ctx);
        seal::util::PolyIter d_itr(dst);
        SEAL_ITERATE(
            seal::util::iter(s_itr, d_itr), ctx.size(), [&](tuple<seal::util::RNSIter, seal::util::RNSIter> I) {
                std::uint64_t mod_index = 0;
                SEAL_ITERATE(
                    seal::util::iter(std::get<0>(I), std::get<1>(I)), positive_wrapped_params.coeff_modulus().size(),
                    [&](auto coefItersJ) {
                        auto mod = parms.coeff_modulus().at(mod_index++);
                        auto s(std::get<0>(coefItersJ));
                        auto double_s(std::get<1>(coefItersJ));

                        for (std::uint64_t i = 0; i < positive_wrapped_params.poly_modulus_degree(); ++i)
                        {
                            auto expected_value = seal::util::multiply_uint_mod(*s, 2, mod);
                            ASSERT_EQ(expected_value, *double_s);
                        }
                    });
            });
    }

    TEST(EvaluatorTest, polynomialMod1)
    {
        std::uint64_t max_poly_size = 128 * 2;
        auto context = make_128deg_context();
        auto parms = context.first_context_data()->parms();
        Evaluator evaluator(context);

        Plaintext ptx("0");
        evaluator.plain_to_coeff_space(ptx, context.first_parms_id());
        evaluator.zero_pad(ptx, context.first_parms_id());

        for (std::uint64_t i = 0; i < ptx.dyn_array().size(); ++i)
        {
            ptx[i] = i % 512;
        }

        // done setup, now we will test the function.
        evaluator.polynomial_mod(ptx);

        //  result needs to be divided into three parts. each part should be of 128 elements, and have the same number
        // throughout. all 3 parts should have different values.
        auto result = plain_to_vector(ptx);

        for (std::uint64_t i = 0; i < 3; ++i)
        {
            for (auto j = i * 128; j < i * 128 + 128; ++j)
            {
                ASSERT_EQ(result[i * 128], result[j]);
            }
        }

        std::unordered_set<std::uint64_t> difference_test;
        for (std::uint64_t i = 0; i < 3; ++i)
        {
            auto j = i * 128;
            difference_test.insert(result[j]);
        }
        ASSERT_EQ(difference_test.size(), 3);
        ASSERT_EQ(result.size(), 128 * 3);
    }

    TEST(EvaluatorTest, polynomialMod2)
    {
        std::uint64_t max_poly_size = 128 * 2;
        auto context = make_128deg_context(1);
        auto parms = context.first_context_data()->parms();
        Evaluator evaluator(context);

        Plaintext ptx("0");
        evaluator.plain_to_coeff_space(ptx, context.first_parms_id());
        evaluator.zero_pad(ptx, context.first_parms_id());

        for (std::uint64_t i = 0; i < ptx.dyn_array().size(); ++i)
        {
            if (i % 128 != i % 256)
            {
                ptx[i] = 0;
                continue;
            }
            ptx[i] = i % 128;
        }

        auto vc = plain_to_vector(ptx);
        // done setup, now we will test the function.
        evaluator.polynomial_mod(ptx);

        //  result needs to be divided into three parts. each part should be of 128 elements, and have the same number
        // throughout. all 3 parts should have different values.
        auto result = plain_to_vector(ptx);
        uint64_t j = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int k = 0; k < 128; ++k)
            {
                ASSERT_EQ(result[j], j % 128);
                j++;
            }
        }
    }

    TEST(EvaluatorTest, ctxMod)
    {
        // shallow test, just ensuring the structure of the ctx remains correct.
        // polynomial mod is checked in polynomialMod
        std::uint64_t max_poly_size = 128 * 2;
        auto context = make_128deg_context(1);
        auto parms = context.first_context_data()->parms();

        Evaluator evaluator(context);
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);

        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);
        evaluator.transform_from_ntt_inplace(encrypted);

        auto non_padded = encrypted;
        evaluator.zero_pad(encrypted);
        evaluator.polynomial_mod(encrypted);

        assert_eq_ciphers(encrypted, non_padded);
    }

    TEST(EvaluatorTest, postiveInverseNtt)
    {
        uint64_t N = 8192;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        parms.set_coeff_modulus(CoeffModulus::Create(N, vector<int>{ 54, 41, 42 }));
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));
        seal::SEALContext context(parms, false, sec_level_type::tc192, 1);
        ASSERT_TRUE(context.parameters_set());

        Evaluator evaluator(context);

        Plaintext ptx = random_plain(context.first_context_data()->parms());

        // creating encryption:
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);
        Decryptor decryptor(context, keygen.secret_key());

        Ciphertext ctx(context);
        encryptor.encrypt(ptx, ctx);
        evaluator.transform_from_ntt_inplace(ctx);

        auto original_ctx = ctx;
        auto ctx_cpy = ctx;
        evaluator.zero_pad(ctx_cpy);

        ctx = evaluator.transform_to_positive_ntt(ctx);
        evaluator.transform_from_positive_ntt_inplace(ctx);

        assert_eq_ciphers(ctx, ctx_cpy);

        evaluator.polynomial_mod(ctx);
        assert_eq_ciphers(ctx, original_ctx);
    }

    TEST(EvaluatorTest, postiveNttWorkingContextExample)
    {
        uint64_t N = 8192;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        parms.set_coeff_modulus(CoeffModulus::Create(N, vector<int>{ 54, 41, 42 }));
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));

        ASSERT_TRUE(seal::SEALContext(parms, false, sec_level_type::tc128, 2).parameters_set());
        ASSERT_TRUE(seal::SEALContext(parms, false, sec_level_type::tc192, 2).parameters_set());
        ASSERT_FALSE(seal::SEALContext(parms, false, sec_level_type::tc256, 2).parameters_set());

        //        seal::SEALContext ctx(parms, false, sec_level_type::tc128, 1);
        //        auto ct = ctx.get_context_data(ctx.positive_wrapped_parms_id());
        //        std::cout << ct->parms().poly_modulus_degree();
    }

    TEST(EvaluatorTest, PaddedCtxAddition)
    {
        uint64_t N = 8192;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        parms.set_coeff_modulus(CoeffModulus::Create(N, vector<int>{ 54, 41, 42 }));
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));
        seal::SEALContext context(parms, false, sec_level_type::tc192, 1);
        ASSERT_TRUE(context.parameters_set());

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);

        Evaluator ev(context);

        seal::Ciphertext res1;
        auto cpy1 = encrypted;
        // TODO: ensure add is correct for our set of parameters.
        ev.add(cpy1, cpy1, res1);
        ev.transform_from_ntt_inplace(res1);

        ev.transform_from_ntt_inplace(encrypted);
        auto cpy2 = ev.transform_to_positive_ntt(encrypted);

        seal::Ciphertext res2;
        ev.add(cpy2, cpy2, res2);

        ev.transform_from_positive_ntt_inplace(res2);
        ev.polynomial_mod(res2);

        assert_eq_ciphers(res1, res2);
    }

    TEST(EvaluatorTest, paddedMultiplication)
    {
        //        GTEST_SKIP();
        uint64_t N = 8192;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        parms.set_coeff_modulus(CoeffModulus::Create(N, SetupObjs::get_8192_positive_ntt_moduli()));
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));
        seal::SEALContext context(parms, false, sec_level_type::tc192, 1);
        ASSERT_TRUE(context.parameters_set());

        auto ptx = random_plain(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);

        auto plain = random_plain(parms);

        Ciphertext encrypted(context);
        encryptor.encrypt(plain, encrypted);
        Evaluator ev(context);

        seal::Ciphertext res1;
        ev.multiply_plain(encrypted, ptx, res1);
        ev.transform_from_ntt_inplace(res1);

        ev.transform_from_ntt_inplace(encrypted);
        ev.transform_to_positive_ntt_inplace(encrypted);

        ev.transform_to_positive_ntt_inplace(ptx);

        seal::Ciphertext res2;
        ev.multiply_plain(encrypted, ptx, res2);

        ev.transform_from_positive_ntt_inplace(res2);
        ev.polynomial_mod(res2);

        assert_eq_ciphers(res1, res2);
    }

    // test whether multiplying ctx with ptx then with ctx is the same with padded_ntt.
    TEST(EvaluatorTest, paddedMultiplication2)
    {
        //        GTEST_SKIP();
        uint64_t N = 8192;
        EncryptionParameters parms(scheme_type::bgv);

        parms.set_poly_modulus_degree(N);
        parms.set_coeff_modulus(CoeffModulus::Create(N, SetupObjs::get_8192_positive_ntt_moduli()));
        parms.set_plain_modulus(PlainModulus::Batching(N, 16 + 1));
        seal::SEALContext context(parms, false, sec_level_type::tc128, 2);
        ASSERT_TRUE(context.parameters_set());

        auto ptx = random_plain(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey pk;
        keygen.create_public_key(pk);
        Encryptor encryptor(context, pk);

        Ciphertext encrypted(context);
        encryptor.encrypt(random_plain(parms), encrypted);
        Evaluator ev(context);

        Ciphertext tmp, res1;
        ev.multiply_plain(encrypted, ptx, tmp);
        ev.multiply(encrypted, tmp, res1);
        ev.transform_from_ntt_inplace(res1);

        ev.transform_from_ntt_inplace(encrypted);
        ev.transform_to_positive_ntt_inplace(encrypted);

        ev.transform_to_positive_ntt_inplace(ptx);

        seal::Ciphertext tmp2, res2;
        ev.multiply_plain(encrypted, ptx, tmp2);
        ev.multiply(encrypted, tmp2, res2);

        ev.transform_from_positive_ntt_inplace(res2);
        ev.polynomial_mod(res2);

        assert_eq_ciphers(res1, res2);
    }
} // namespace sealtest
