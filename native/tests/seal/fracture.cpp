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

        static SetupObjs New(std::uint64_t N = 4096 * 2, int logt = 20)
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

            shared_ptr<UniformRandomGenerator> gen = prng();

            auto mod = enc_params.plain_modulus().value();
            std::generate(v.begin(), v.end(), [gen = std::move(gen), &mod]() { return gen->generate() % mod; });

            for (std::uint64_t i = 0; i < enc_params.poly_modulus_degree(); ++i)
            {
                p[i] = v[i];
            }
            return p;
        }

        shared_ptr<UniformRandomGenerator> prng() const
        {
            Blake2xbPRNGFactory factory;
            array<uint64_t, prng_seed_uint64_count> seed{ 1, 2, 3, 4, 5, 6, 7, 8 };
            auto gen = factory.create(seed);
            return gen;
        }

        seal::Plaintext random_ntt_plaintext() const
        {
            auto p = random_plaintext();
            auto pid = essence.ctx.first_parms_id();

            evaluator.plain_to_coeff_space(p, pid);
            evaluator.transform_plain_in_coeff_space_to_ntt_inplace(p, pid);

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

    seal::Ciphertext ctxWithSize3(const SetupObjs &all)
    {
        seal::Ciphertext ctx1 = all.random_ciphertext();
        seal::Ciphertext ctx2 = all.random_ciphertext();
        all.evaluator.multiply_inplace(ctx1, ctx2);
        return ctx1;
    }

    void add_inplace_ctx_fractures(
        seal::fractures::CiphertextShredder &ctxshred1, seal::fractures::CiphertextShredder &ctxshred2)
    {
        for (std::uint64_t i = 0; i < ctxshred1.num_fractures(); ++i)
        {
            ctxshred1[i] += ctxshred2[i];
        }
    }

    template <typename T, typename U>
    void verify_dims(const seal::util::matrix<T> &a, const seal::util::matrix<U> &b)
    {
        if (a.cols != b.rows)
        {
            throw std::invalid_argument("MatrixOperations::multiply: left matrix cols != right matrix rows");
        }
    }

    // Doing matmul...
    template <typename T>
    seal::util::matrix<fractures::CiphertextFracture> multiplyMatrices(
        const seal::util::matrix<fractures::CiphertextFracture> &a, const seal::util::matrix<T> &b)
    {
        verify_dims(a, b);

        seal::util::matrix<fractures::CiphertextFracture> result(a.rows, b.cols);

        std::uint64_t poly_sz = 2;
        if constexpr ((std::is_same_v<T, seal::fractures::CiphertextFracture>))
        {
            poly_sz = 3;
        }

        auto index = a(0, 0).index;
        auto coeff_modulus = a(0, 0).coeff_modulus;
        auto coeff_count = a(0, 0).poly_fracs[0].coeff_count;

        for (std::uint64_t i = 0; i < result.rows; ++i)
        {
            for (std::uint64_t j = 0; j < result.cols; ++j)
            {
                // create an empty sum variable (using the parameters of a.data[0])
                auto sum = seal::fractures::CiphertextFracture::Empty(poly_sz, index, coeff_count, coeff_modulus);
                for (std::uint64_t k = 0; k < a.cols; ++k)
                {
                    sum += a(i, k) * b(k, j);
                }
                result(i, j) = sum;
            }
        }

        return result;
    }

    template <typename T>
    seal::util::matrix<Ciphertext> multiplyMatrices(
        SetupObjs &all, const seal::util::matrix<Ciphertext> &a, const seal::util::matrix<T> &b)
    {
        verify_dims(a, b);

        seal::util::matrix<Ciphertext> result(a.rows, b.cols);

        for (std::uint64_t i = 0; i < result.rows; ++i)
        {
            for (std::uint64_t j = 0; j < result.cols; ++j)
            {
                // create an empty sum variable (using the parameters of a.data[0])
                seal::Ciphertext sum(all.context);
                seal::Ciphertext tmp;
                all.evaluator.transform_to_ntt_inplace(sum);

                for (std::uint64_t k = 0; k < a.cols; ++k)
                {
                    if constexpr ((std::is_same_v<T, seal::Plaintext>))
                    {
                        all.evaluator.multiply_plain(a(i, k), b(k, j), tmp);
                    }
                    else
                    {
                        all.evaluator.multiply(a(i, k), b(k, j), tmp);
                    }
                    all.evaluator.add(tmp, sum, sum);
                }

                result(i, j) = sum;
            }
        }

        return result;
    }

    vector<util::matrix<fractures::CiphertextFracture>> fracture_matrix(
        const SetupObjs &all, util::matrix<Ciphertext> &ctx_mat, std::uint64_t num_fractures)
    {
        // take each ctx and shred it using ctxShredder. then create matrices with the same dims.
        std::vector<fractures::CiphertextShredder> ctx_shredders;
        for (auto &ctx_iter : ctx_mat.data)
        {
            ctx_shredders.emplace_back(ctx_iter, all.essence, num_fractures);
        }
        seal::util::matrix<fractures::CiphertextShredder> shredding(ctx_mat.rows, ctx_mat.cols, ctx_shredders);

        std::vector<seal::util::matrix<fractures::CiphertextFracture>> fractured_matrices;
        for (std::uint64_t frac_num = 0; frac_num < num_fractures; ++frac_num)
        {
            std::vector<fractures::CiphertextFracture> tmp;
            for (std::uint64_t row = 0; row < ctx_mat.rows; ++row)
            {
                for (std::uint64_t col = 0; col < ctx_mat.cols; ++col)
                {
                    tmp.push_back(std::move(shredding(row, col)[frac_num]));
                }
            }
            fractured_matrices.emplace_back(ctx_mat.rows, ctx_mat.cols, tmp);
        }
        return fractured_matrices;
    }

    TEST(FracturedOps, PolynomialEvaluation)
    {
        auto all = SetupObjs::New(1 << 12);
        seal::fractures::PolynomialEvaluator pe(all.essence);

        auto ptx = all.random_plaintext();
        auto ctx = all.random_ciphertext();

        seal::Ciphertext res;
        all.evaluator.multiply_plain(ctx, ptx, res);

        all.evaluator.plain_to_coeff_space(ptx, all.essence.ctx.first_parms_id());
        auto point1 = pe.evaluate(ptx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

        all.evaluator.transform_from_ntt_inplace(ctx);
        auto cpoint = pe.evaluate(ctx, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

        all.evaluator.transform_from_ntt_inplace(res);
        auto rpoint = pe.evaluate(res, std::vector<std::uint64_t>{ 1, 2, 3, 4, 5 });

        auto actual = point1 * cpoint;
        // compare:
        ASSERT_TRUE(actual == rpoint);
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

    std::vector<seal::Ciphertext> random_ctxs(const SetupObjs &all, int size)
    {
        std::vector<seal::Ciphertext> ctxs;
        for (int i = 0; i < size; ++i)
        {
            ctxs.push_back(all.random_ciphertext());
        }
        return ctxs;
    }
    std::vector<seal::Plaintext> random_ptxs(const SetupObjs &all, int size)
    {
        std::vector<seal::Plaintext> ptxs;
        for (int i = 0; i < size; ++i)
        {
            ptxs.push_back(all.random_ntt_plaintext());
        }
        return ptxs;
    }

    TEST(FracturedOps, CtxVectorPeoductWithPlaintextsVector)
    {
        auto all = SetupObjs::New();
        std::uint64_t num_fractures = 256;
        auto num_ctxs = 50;

        std::vector<seal::Ciphertext> ctxs = random_ctxs(all, num_ctxs);
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

        seal::util::matrix<Ciphertext> left(r, c, random_ctxs(all, n));
        seal::util::matrix<Ciphertext> right(c, r, random_ctxs(all, n));

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

        auto expected = random_ctxs(all, 200);

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
} // namespace sealtest