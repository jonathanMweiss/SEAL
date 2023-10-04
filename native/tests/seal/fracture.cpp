// Test for fractured operations by JonathanWeiss.

#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractured_ciphertext.h"
#include "seal/fractured_polynomial.h"
#include "seal/keygenerator.h"
#include "seal/modulus.h"
#include "seal/poly_eval.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <sstream>
#include <string>
#include <utility>
#include "frac_test_helpers.cpp"
#include "gtest/gtest.h"

using namespace seal;
using namespace std;
namespace sealtest::fracture
{

    TEST(PolyEvaluate, PtxIsScalar)
    {
        auto all = SetupObjs::New(1 << 12);
        seal::fractures::PolynomialEvaluator pe(all.essence);

        for (int i = 0; i < 10; ++i)
        {
            auto ctx = all.random_ciphertext();

            all.evaluator.transform_from_ntt_inplace(ctx);
            auto actual = pe.evaluate(ctx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            seal::Plaintext p1(all.essence.parms.poly_modulus_degree());
            p1[0] = std::uint64_t(all.prng()->generate()) % all.essence.parms.plain_modulus().value();

            all.evaluator.plain_to_coeff_space(p1, all.essence.ctx.first_parms_id());
            auto point = pe.evaluate(p1, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            all.evaluator.transform_plain_in_coeff_space_to_ntt_inplace(p1, all.essence.ctx.first_parms_id());

            all.evaluator.transform_to_ntt_inplace(ctx);
            all.evaluator.multiply_plain_inplace(ctx, p1);
            all.evaluator.transform_from_ntt_inplace(ctx);

            auto expected = pe.evaluate(ctx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            ASSERT_TRUE(expected == (actual * point));
        }
    }

    TEST(PolyEvaluate, ctxMultPtx)
    {
        // TODO: maybe im not getting the right value?
        //   it seems to work when RNS-size ==1.
        auto all = SetupObjs::New(1 << 11);
        seal::fractures::PolynomialEvaluator pe(all.essence);
        auto prng = all.prng();

        auto ctx1 = all.random_ciphertext();
        auto ptx = all.random_plaintext();

        seal::Ciphertext mult_res;

        auto ctx_data = all.context.get_context_data(all.context.first_parms_id());
        auto root = ctx_data->small_ntt_tables()->get_root();
        auto rns_tool = ctx_data->rns_tool();
        auto root_rns = std::vector<std::uint64_t>(all.enc_params.coeff_modulus().size(), root);
        rns_tool->base_q()->decompose(&(root_rns[0]), MemoryManager::GetPool());
        std::cout << "root_rns size: " << root_rns.size() << std::endl;

        all.evaluator.plain_to_coeff_space(ptx, all.essence.ctx.first_parms_id());
        auto p1 = pe.evaluate(ptx, root_rns);

        all.evaluator.transform_plain_in_coeff_space_to_ntt_inplace(ptx, all.essence.ctx.first_parms_id());
        all.evaluator.multiply_plain(ctx1, ptx, mult_res);

        all.evaluator.transform_from_ntt_inplace(ctx1);
        auto c1 = pe.evaluate(ctx1, root_rns);

        auto actual = c1 * p1;

        //
        all.evaluator.transform_from_ntt_inplace(mult_res);
        auto expected = pe.evaluate(mult_res, root_rns);
        ASSERT_TRUE(expected == actual);
    }

    TEST(PolyEvaluate, ctxAddCtx)
    {
        auto all = SetupObjs::New(1 << 10);
        seal::fractures::PolynomialEvaluator pe(all.essence);

        //        for (int i = 0; i < 500; ++i)
        //        {

        auto ctx_data = all.context.get_context_data(all.context.first_parms_id());
        auto root = ctx_data->small_ntt_tables()->get_root();
        auto rns_tool = ctx_data->rns_tool();
        auto root_rns = std::vector<std::uint64_t>(all.enc_params.coeff_modulus().size(), root);
        rns_tool->base_q()->decompose(&(root_rns[0]), MemoryManager::GetPool());
        std::cout << "root_rns size: " << root_rns.size() << std::endl;

        auto ctx1 = all.random_ciphertext();
        auto ctx2 = all.random_ciphertext();
        seal::Ciphertext addition_res;

        all.evaluator.multiply(ctx1, ctx2, addition_res);

        all.evaluator.transform_from_ntt_inplace(ctx1);
        all.evaluator.transform_from_ntt_inplace(ctx2);

        auto c1 = pe.evaluate(ctx1, root_rns);
        auto c2 = pe.evaluate(ctx2, root_rns);
        auto actual = c1 * c2;

        all.evaluator.transform_from_ntt_inplace(addition_res);
        auto expected = pe.evaluate(addition_res, root_rns);
        ASSERT_TRUE(expected == actual);
        //        }
    }

    namespace operations
    {
        TEST(FracturedOps, PolyEvalEquals)
        {
            auto all = SetupObjs::New(1 << 12);
            seal::fractures::PolynomialEvaluator pe(all.essence);

            auto ptx = all.random_plaintext();
            auto ctx = all.random_ciphertext();

            seal::Ciphertext res;
            all.evaluator.multiply_plain(ctx, ptx, res);

            all.evaluator.plain_to_coeff_space(ptx, all.essence.parms.parms_id());
            auto point1 = pe.evaluate(ptx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            all.evaluator.transform_from_ntt_inplace(ctx);
            auto cpoint = pe.evaluate(ctx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            all.evaluator.transform_from_ntt_inplace(res);
            auto rpoint = pe.evaluate(res, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

            cpoint *= point1;
            // compare:
            ASSERT_TRUE(cpoint == rpoint);
            std::cout << "lols" << std::endl;
        }

        TEST(FracturedOps, MultPlain)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext encrypted_ntt = all.random_ciphertext();
            auto ptx = all.random_ntt_plaintext();

            seal::fractures::CiphertextShredder cshredder(encrypted_ntt, all.essence, num_fractures);
            seal::fractures::Polynomial pshredder(ptx, all.essence, num_fractures);

            // perform fractured multiplication:
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                cshredder[i] *= pshredder.get_fracture(i);
            }

            // regular multiplication and then we'll compare them:
            all.evaluator.multiply_plain_inplace(encrypted_ntt, ptx);
            ASSERT_TRUE(!encrypted_ntt.is_transparent());

            all.evaluator.sub_inplace(encrypted_ntt, cshredder.into_ciphertext());
            ASSERT_TRUE(encrypted_ntt.is_transparent());

            seal::Plaintext ptx_res;
            all.decryptor.decrypt(cshredder.into_ciphertext(), ptx_res);
        }

        TEST(FracturedOps, addCtxSize2)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = all.random_ciphertext();
            seal::Ciphertext ctx2 = all.random_ciphertext();

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            add_inplace_ctx_fractures(ctxshred1, ctxshred2);

            all.evaluator.add_inplace(ctx1, ctx2);
            ASSERT_TRUE(!ctx1.is_transparent());

            all.evaluator.sub_inplace(ctx1, ctxshred1.into_ciphertext());
            ASSERT_TRUE(ctx1.is_transparent());
        }

        TEST(FracturedOps, addCtxSize3)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = ctxWithSize3(all);
            seal::Ciphertext ctx2 = ctxWithSize3(all);

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            add_inplace_ctx_fractures(ctxshred1, ctxshred2);

            all.evaluator.add_inplace(ctx1, ctx2);
            ASSERT_TRUE(!ctx1.is_transparent());

            all.evaluator.sub_inplace(ctx1, ctxshred1.into_ciphertext());
            ASSERT_TRUE(ctx1.is_transparent());
        }

        TEST(FracturedOps, AddTwoCtxFails)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = all.random_ciphertext();
            seal::Ciphertext ctx2 = all.random_ciphertext();

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            for (std::uint64_t i = 0; i < num_fractures - 1; ++i)
            {
                ctxshred1[i] += ctxshred2[i];
            }

            all.evaluator.add_inplace(ctx1, ctx2);
            ASSERT_TRUE(!ctx1.is_transparent());

            all.evaluator.sub_inplace(ctx1, ctxshred1.into_ciphertext());
            ASSERT_FALSE(ctx1.is_transparent());
        }

        TEST(FracturedOps, MultiplyFails)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext encrypted_ntt = all.random_ciphertext();
            auto ptx = all.random_ntt_plaintext();

            seal::fractures::CiphertextShredder cshredder(encrypted_ntt, all.essence, num_fractures);
            seal::fractures::Polynomial pshredder(ptx, all.essence, num_fractures);

            // perform fractured multiplication:
            for (std::uint64_t i = 0; i < num_fractures - 1; ++i)
            {
                cshredder[i] *= pshredder.get_fracture(i);
            }

            // regular multiplication and then we'll compare them:
            all.evaluator.multiply_plain_inplace(encrypted_ntt, ptx);
            ASSERT_TRUE(!encrypted_ntt.is_transparent());

            all.evaluator.sub_inplace(encrypted_ntt, cshredder.into_ciphertext());
            ASSERT_FALSE(encrypted_ntt.is_transparent());
        }

        TEST(FracturedOps, CtxVectorPeoductWithPlaintextsVector)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;
            auto num_ctxs = 50;

            std::vector<seal::Ciphertext> ctxs = random_ctx_vector(all, num_ctxs);
            std::vector<seal::Plaintext> ptxs = random_ptxs(all, num_ctxs);

            std::vector<seal::fractures::CiphertextShredder> ctxs_fracs;
            std::vector<seal::fractures::Polynomial> ptxs_fracs;

            for (std::uint64_t i = 0; i < std::uint64_t(num_ctxs); ++i)
            {
                ctxs_fracs.emplace_back(ctxs[i], all.essence, num_fractures);
                ptxs_fracs.emplace_back(ptxs[i], all.essence, num_fractures);
            }

            // multiply fractured:
            for (std::uint64_t i = 0; i < ctxs.size(); ++i)
            {
                all.evaluator.multiply_plain_inplace(ctxs[i], ptxs[i]);
                for (std::uint64_t j = 0; j < num_fractures; ++j)
                {
                    ctxs_fracs[i][j] *= ptxs_fracs[i].get_fracture(j);
                }
            }

            // summing into the first element of the vector.
            for (std::uint64_t i = 1; i < ctxs.size(); ++i)
            {
                all.evaluator.add_inplace(ctxs[0], ctxs[i]);
                add_inplace_ctx_fractures(ctxs_fracs[0], ctxs_fracs[i]);
            }

            auto actual_result = ctxs_fracs[0].into_ciphertext();
            ASSERT_TRUE(!actual_result.is_transparent());
            ASSERT_TRUE(!ctxs[0].is_transparent());

            all.evaluator.sub_inplace(actual_result, ctxs[0]);
            ASSERT_TRUE(actual_result.is_transparent());
        }

        TEST(FracturedOps, ctxCtxMult)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = all.random_ciphertext();
            seal::Ciphertext ctx2 = all.random_ciphertext();

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            seal::Ciphertext expected;
            all.evaluator.multiply(ctx1, ctx2, expected);
            ASSERT_TRUE(!expected.is_transparent());

            for (std::uint64_t i = 0; i < ctxshred1.num_fractures(); ++i)
            {
                ctxshred1[i] *= ctxshred2[i];
            }

            // into ctx is probably okay because i've seen it decompose correctly when it has size >2.
            all.evaluator.sub_inplace(expected, ctxshred1.into_ciphertext());
            ASSERT_TRUE(expected.is_transparent());
        }

        TEST(FracturedOps, multCtxThenAdd)
        {
            // ctx1*ctx2 and then ctx1+ctx1.

            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = all.random_ciphertext();
            seal::Ciphertext ctx2 = all.random_ciphertext();

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            all.evaluator.multiply_inplace(ctx1, ctx2);
            all.evaluator.add_inplace(ctx1, ctx1);

            // ctx1*ctx2
            for (std::uint64_t i = 0; i < ctxshred1.num_fractures(); ++i)
            {
                ctxshred1[i] *= ctxshred2[i];
            }

            // ctx1+ctx1
            for (std::uint64_t i = 0; i < ctxshred1.num_fractures(); ++i)
            {
                ctxshred1[i] += ctxshred1[i];
            }

            // compare:
            all.evaluator.sub_inplace(ctx1, ctxshred1.into_ciphertext());
            ASSERT_TRUE(ctx1.is_transparent());
        }

        TEST(FracturedOps, ctxVectorDotProd)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;
            std::uint64_t r = 1;
            std::uint64_t c = 50;
            auto n = int(r * c);

            seal::util::matrix<Ciphertext> left(r, c, random_ctx_vector(all, n));
            seal::util::matrix<Ciphertext> right(c, r, random_ctx_vector(all, n));

            auto shred_left = fracture_matrix(all, left, num_fractures);
            auto shred_right = fracture_matrix(all, right, num_fractures);

            // multiply and store the results in the entries of the left vector.
            auto expected = multiplyMatrices(all, left, right);
            std::vector<seal::util::matrix<seal::fractures::CiphertextFracture>> result_fracs;

            seal::fractures::CiphertextShredder result(all.essence, num_fractures);
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                auto tmp = multiplyMatrices(shred_left[i], shred_right[i]);
                ASSERT_TRUE(tmp.cols == 1 && tmp.rows == 1);
                result.set_fracture(i, tmp(0, 0));
            }

            //    compare:
            all.evaluator.sub_inplace(expected(0, 0), result.into_ciphertext());
            ASSERT_TRUE(expected(0, 0).is_transparent());
            std::cout << "YAY" << std::endl;
        }

        TEST(FracturedOps, manyAdditions)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            auto expected = random_ctx_vector(all, 200);

            std::vector<seal::fractures::CiphertextShredder> actual;
            std::vector<seal::fractures::CiphertextShredder> actual2;

            // fill shredders:
            for (auto &ctx : expected)
            {
                actual.emplace_back(ctx, all.essence, num_fractures);
            }

            // sum into the first element of the vector.
            for (std::uint64_t i = 1; i < expected.size(); ++i)
            {
                all.evaluator.add_inplace(expected[0], expected[i]);
                for (std::uint64_t j = 0; j < num_fractures; ++j)
                {
                    actual[0][j] += actual[i][j];
                }
            }

            // compare:
            all.evaluator.sub_inplace(expected[0], actual[0].into_ciphertext());
            ASSERT_TRUE(expected[0].is_transparent());
        }

        TEST(FracturedOps, ptxCtxMultThenCtxMult)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            auto ctx1 = all.random_ciphertext();
            auto ptx = all.random_ntt_plaintext();
            auto ctx2 = all.random_ciphertext();

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.essence, num_fractures);
            seal::fractures::Polynomial pshred(ptx, all.essence, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.essence, num_fractures);

            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                ctxshred1[i] *= pshred.get_fracture(i);
                ctxshred1[i] *= ctxshred2[i];
            }

            all.evaluator.multiply_plain_inplace(ctx1, ptx);
            all.evaluator.multiply_inplace(ctx1, ctx2);

            auto result = ctxshred1.into_ciphertext();
            all.evaluator.sub_inplace(result, ctx1);
            ASSERT_TRUE(result.is_transparent());
        }

        TEST(FracturedOps, ptxMatrixMultWithQueryiesFromLeftAndRight)
        {
            ASSERT_TRUE(false);
        }
    } // namespace operations

    namespace serialization
    {
        TEST(FracturedOps, serializePolyFrac)
        {
            auto all = SetupObjs::New();
            auto gen = all.prng();

            auto index = gen->generate();
            auto coeff_count = gen->generate() % 15;
            auto modulus_rns_size = gen->generate() % 6;
            fractures::PolynomialFracture pf(index, coeff_count, modulus_rns_size);

            for (auto &d : pf.rns_coefficients.data)
            {
                d = gen->generate();
            }

            stringstream stream;
            pf.save(stream);

            fractures::PolynomialFracture pf2(index, coeff_count, modulus_rns_size);
            auto st = stream.str();
            pf2.load((seal::seal_byte *)&st[0], stream.str().length());
            ASSERT_TRUE(pf2 == pf);

            fractures::PolynomialFracture pf3(index, coeff_count, modulus_rns_size);
            pf3.load(stream);
            ASSERT_TRUE(pf3 == pf);
        }

        TEST(FracturedOps, serializeCipherFrac)
        {
            auto all = SetupObjs::New();
            auto gen = all.prng();

            std::uint64_t index = gen->generate();
            std::uint64_t coeff_count = gen->generate() % 256;
            auto modulus = all.essence.parms.coeff_modulus();
            std::uint64_t num_polys = gen->generate() % 10;
            auto cf = fractures::CiphertextFracture::Empty(num_polys, index, coeff_count, modulus);
            for (std::uint64_t i = 0; i < num_polys; ++i)
            {
                fractures::PolynomialFracture pf(index, coeff_count, modulus.size());
                for (auto &d : pf.rns_coefficients.data)
                {
                    d = gen->generate();
                }
                cf.poly_fracs[i] = pf;
            }

            stringstream stream;
            cf.save(stream);

            auto cf2 = fractures::CiphertextFracture::Empty(num_polys, index, coeff_count, modulus);

            auto st = stream.str();
            cf2.load((seal::seal_byte *)&st[0], stream.str().length());
            ASSERT_TRUE(cf2 == cf);

            auto cf3 = fractures::CiphertextFracture::Empty(num_polys, index, coeff_count, modulus);
            cf3.load(stream);
            ASSERT_TRUE(cf3 == cf);
        }

    } // namespace serialization
} // namespace sealtest::fracture