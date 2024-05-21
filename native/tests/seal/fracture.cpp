
#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractured_ciphertext.h"
#include "seal/fractured_polynomial.h"
#include "seal/keygenerator.h"
#include "seal/matmul.hpp"
#include "seal/modulus.h"
#include "seal/poly_eval.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include "frac_test_helpers.cpp"
#include "gtest/gtest.h"

using namespace seal;
using namespace std;
namespace sealtest::fracture
{

    namespace polyval
    {

        TEST(PolyEvaluate, ctxXctx)
        {
            auto all = SetupObjs::New();
            auto r = root;

            ASSERT_TRUE(random_ctx_ctx_sz(all, r));
        }

        TEST(PolyEvaluate, ctxXptx)
        {
            auto all = SetupObjs::New();

            auto r = root;
            ASSERT_TRUE(random_ptx_ctx_sz(all, r));
        }
        TEST(PolyEvaluate, exhaustiveSZCtxCtx)
        {
            GTEST_SKIP();

            auto all = SetupObjs::New();
            auto roots = generate_roots(all, 1 << 14);

            for (auto &r : roots)
            {
                ASSERT_TRUE(random_ctx_ctx_sz(all, r));
            }
        }

        TEST(PolyEvaluate, exhaustiveSZptxCtx)
        {
            GTEST_SKIP();

            auto all = SetupObjs::New();
            auto roots = generate_roots(all, 1 << 14);

            for (auto &r : roots)
            {
                ASSERT_TRUE(random_ptx_ctx_sz(all, r));
            }
        }

        /**
         * Generate a root of unity for a given modulus and n.
         * @param m The modulus.
         * @param n The n-th root of unity.
         * @return The root of unity.
         */
        std::uint64_t gen_root_of_unity(const seal::Modulus &m, std::uint64_t n)
        {
            // $a$ is primitive element (root) if each element that isn't $0$ can be expressed as
            // $a^k$ for some $k$.
            std::uint64_t a;
            seal::util::try_minimal_primitive_root(n, m, a);

            // a^(q-1)/n is an n-th root of unity.
            // proof: we know that $a^(q-1)$ is 1 (fermat's little theorem),
            // as a result, $a^k\ne 1$ for all $k < q-1$. otherwise there is a cycle that doesn't go through every
            // element other than 0 in the field, in contradiction to $a$ being a generator/ primitive element.
            // Thus, $w=a^(q-1)/n$ is an n-th root of unity, because w^n = a^(q-1) = 1 and w^k != 1 for all k < n.
            return seal::util::exponentiate_uint_mod(a, (m.value() - 1 / n), m);
        }
        
        std::vector<std::uint64_t> generate_root_of_unity(const seal::SEALContext &cntx, seal::parms_id_type pid)
        {
            auto parms = cntx.get_context_data(pid)->parms();

            std::vector<std::uint64_t> r;
            for (std::uint64_t i = 0; i < parms.coeff_modulus().size(); ++i)
            {
                auto m = parms.coeff_modulus()[i];
                r.emplace_back(gen_root_of_unity(m, 1 << 14));
            }
            return r;
        }

        std::vector<std::uint64_t> generate_root_of_unity(const seal::SEALContext &cntx)
        {
            return generate_root_of_unity(cntx, cntx.first_parms_id());
        }

        TEST(PolyEvaluate, genpoint)
        {
            auto all = SetupObjs::New();
            auto parms = all.context.first_context_data()->parms();

            auto r = generate_root_of_unity(all.context);
            if (random_ctx_ctx_sz(all, r))
            {
                std::cout << "Success" << std::endl;
            }
            else
            {
                std::cout << "Fail" << std::endl;
            }
        }
        TEST(PolyEvaluate, PIR_SZ)
        {
            auto all = SetupObjs::New();

            std::uint64_t vec_size = 5;

            seal::util::matrix<Ciphertext> query_right(vec_size, 1, random_ctx_vector(all, int(vec_size)));
            seal::util::matrix<Ciphertext> query_left(1, vec_size, random_ctx_vector(all, int(vec_size)));
            apply_on_each_element<Ciphertext>(
                query_right, [&](Ciphertext &c) -> void { all.evaluator.transform_from_ntt_inplace(c); });
            //            //            seal::util::matrix<Plaintext> db(vec_size, vec_size, random_ptxs(all,
            //            int(vec_size *
            //            //            vec_size)));
            //
            //            evaluate_matrix<Ciphertext, seal::fractures::EvaluatedCipherPoint>(all, query_right);
            //            // Create evaluated matrices:
        }
    } // namespace polyval

    namespace operations
    {

        TEST(FracturedOps, plaintextFracture)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext encrypted_ntt = all.random_ciphertext();
            auto ptx = all.random_ntt_plaintext();

            seal::fractures::Polynomial pshredder(ptx, all.context, num_fractures);
            seal::fractures::CiphertextShredder cshredder(encrypted_ntt, all.context, num_fractures);
        }

        TEST(FracturedOps, MultPlain)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext encrypted_ntt = all.random_ciphertext();
            auto ptx = all.random_ntt_plaintext();

            seal::fractures::CiphertextShredder cshredder(encrypted_ntt, all.context, num_fractures);
            seal::fractures::Polynomial pshredder(ptx, all.context, num_fractures);

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

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

            add_inplace_ctx_fractures(ctxshred1, ctxshred2);

            all.evaluator.add_inplace(ctx1, ctx2);
            ASSERT_TRUE(!ctx1.is_transparent());

            all.evaluator.sub_inplace(ctx1, ctxshred1.into_ciphertext());
            ASSERT_TRUE(ctx1.is_transparent());
        }

        TEST(FracturedOps, ctxSize3Budget)
        {
            auto all = SetupObjs::New();

            auto ptx = all.random_ntt_plaintext();
            auto ctx1 = all.random_ciphertext();

            all.evaluator.multiply_plain_inplace(ctx1, ptx);
            all.evaluator.multiply_inplace(ctx1, ctx1);

            //            all.decryptor.decrypt(ctx1, p);
            std::cout << "=====" << std::endl;
            //            std::cout << p[0] << std::endl;
            std::cout << all.decryptor.invariant_noise_budget(ctx1) << std::endl;
            ASSERT_TRUE(all.decryptor.invariant_noise_budget(ctx1) > 0);
        }

        TEST(FracturedOps, addCtxSize3)
        {
            auto all = SetupObjs::New();
            std::uint64_t num_fractures = 256;

            seal::Ciphertext ctx1 = ctxWithSize3(all);
            seal::Ciphertext ctx2 = ctxWithSize3(all);

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

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

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

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

            seal::fractures::CiphertextShredder cshredder(encrypted_ntt, all.context, num_fractures);
            seal::fractures::Polynomial pshredder(ptx, all.context, num_fractures);

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
                ctxs_fracs.emplace_back(ctxs[i], all.context, num_fractures);
                ptxs_fracs.emplace_back(ptxs[i], all.context, num_fractures);
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

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

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

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

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

            seal::fractures::CiphertextShredder result(all.context, num_fractures);
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                auto tmp = multiply_matrices(shred_left[i], shred_right[i]);
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
                actual.emplace_back(ctx, all.context, num_fractures);
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

            seal::fractures::CiphertextShredder ctxshred1(ctx1, all.context, num_fractures);
            seal::fractures::Polynomial pshred(ptx, all.context, num_fractures);
            seal::fractures::CiphertextShredder ctxshred2(ctx2, all.context, num_fractures);

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

        /**
         * This is the final test for fracture update. Ensuring we can safely perform a fractured PIR query.
         * Assume the fractures are split across 256 machines, where the ctxs are changing once every x minutes (e.g.,
         *  10 < x < 30).
         */
        TEST(FracturedOps, PIRQuery)
        {
            auto all = SetupObjs::New();

            std::uint64_t vec_size = 5;

            seal::util::matrix<Ciphertext> query_right(vec_size, 1, random_ctx_vector(all, int(vec_size)));
            seal::util::matrix<Ciphertext> query_left(1, vec_size, random_ctx_vector(all, int(vec_size)));

            seal::util::matrix<Plaintext> db(vec_size, vec_size, random_ptxs(all, int(vec_size * vec_size)));

            auto ctx_vector = multiplyMatrices(all, db, query_right);
            auto response = multiplyMatrices(all, query_left, ctx_vector);

            // fracture:
            std::uint64_t num_fractures = 8;
            auto frac_query_right = fracture_matrix(all, query_right, num_fractures);
            auto frac_query_left = fracture_matrix(all, query_left, num_fractures);
            auto frac_db = fracture_matrix(all, db, num_fractures);

            // multiply:
            seal::fractures::CiphertextShredder composit(all.context, num_fractures);
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                auto frac_ctx_vector = multiply_matrices(frac_db[i], frac_query_right[i]);
                auto frac_response = multiply_matrices(frac_query_left[i], frac_ctx_vector);
                composit.set_fracture(i, frac_response(0, 0));
            }

            // compare:
            all.evaluator.sub_inplace(response(0, 0), composit.into_ciphertext());
            ASSERT_TRUE(response(0, 0).is_transparent());
            std::cout << "Noise budget: " << all.decryptor.invariant_noise_budget(response(0, 0)) << std::endl;
            ASSERT_TRUE(all.decryptor.invariant_noise_budget(response(0, 0)) > 0);
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

            fractures::PolynomialFracture pf3(0, 0, modulus_rns_size);
            pf3.load(stream);
            ASSERT_TRUE(pf3 == pf);
        }

        TEST(FracturedOps, serializeSparsePtx)
        {
            auto all = SetupObjs::New();

            auto ptx = all.random_ntt_plaintext();
            for (unsigned long i = 0; i < ptx.coeff_count(); ++i)
            {
                ptx[i] = i % 256;
            }

            for (uint64_t i = 4; i <= 1024; i = i * 2)
            {
                seal::fractures::Polynomial pshredder1(ptx, all.context, i);

                for (uint64_t frac = 0; frac < i; ++frac)
                {
                    auto pf = pshredder1[frac];
                    stringstream stream;
                    pf.save(stream);

                    fractures::PolynomialFracture pf2(frac, 0, 0);
                    auto st = stream.str();
                    pf2.load((seal::seal_byte *)&st[0], stream.str().length());
                    ASSERT_TRUE(pf2 == pf);
                }
            }
        }

        TEST(FracturedOps, serializeTestSizes)
        {
            auto all = SetupObjs::New();

            auto ptx = all.random_ntt_plaintext();
            std::uint64_t prev_size = 1;
            prev_size = prev_size << 63;
            for (int i = 4; i <= 1024; i = i * 2)
            {
                seal::fractures::Polynomial pshredder1(ptx, all.context, i);

                stringstream stream1;
                pshredder1[0].save(stream1);

                std::cout << "stream1 (frac 0 of " << i << " fractures): " << stream1.str().length() << std::endl;
                ASSERT_TRUE(stream1.str().length() < prev_size);
                prev_size = stream1.str().length();
            }
        }

        TEST(FracturedOps, serializeCipherFrac)
        {
            auto all = SetupObjs::New();
            auto gen = all.prng();

            std::uint64_t index = gen->generate();
            std::uint64_t coeff_count = gen->generate() % 256;
            auto modulus = all.context.first_context_data()->parms().coeff_modulus();
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