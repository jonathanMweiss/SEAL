#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractures.h"
#include "seal/keygenerator.h"
#include "seal/matmul.hpp"
#include "seal/modulus.h"
#include "seal/poly_eval.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
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
        seal::Evaluator evaluator;
        seal::KeyGenerator keygen;
        seal::SecretKey secret_key;
        seal::Encryptor encryptor;
        seal::Decryptor decryptor;

        /**
         * The following parameters are chosen to support ctx-ctx multiplication and still have budget left.
         * @param N
         * @param logt
         * @return
         */
        static SetupObjs New(std::uint64_t N = 4096 * 2, int logt = 16)
        {
            seal::EncryptionParameters enc(seal::scheme_type::bgv);

            enc.set_poly_modulus_degree(N);
            //            enc.set_coeff_modulus(seal::CoeffModulus::BFVDefault(N, sec_level_type::tc192));

            auto tmp = std::vector<int>{ 54, 41, 42 };
            auto o = seal::CoeffModulus::Create(N, tmp);
            enc.set_coeff_modulus(o);
            enc.set_plain_modulus(seal::PlainModulus::Batching(N, logt + 1));
            // 218 is 128-bit secure.
            // 152 is 192-bit secure.
            //            auto tmp = std::vector<int>{ 54, 50,40};
            //            auto tmp = std::vector<int>{ 54, 41, 42 };
            //            auto o = seal::CoeffModulus::Create(N, tmp);
            //            enc.set_coeff_modulus(o);
            //            enc.set_plain_modulus(seal::PlainModulus::Batching(N, logt + 1));

            return SetupObjs(enc);
        }

        explicit SetupObjs(seal::EncryptionParameters encryption_params)
            : enc_params(std::move(encryption_params)), context(enc_params, true), evaluator(context), keygen(context),
              secret_key(keygen.secret_key()), encryptor(context, secret_key), decryptor(context, secret_key){};

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
            auto pid = context.first_parms_id();

            evaluator.plain_to_coeff_space(p, pid);
            evaluator.transform_plain_in_coeff_space_to_ntt_inplace(p, pid);

            return p;
        }

        seal::Ciphertext random_ciphertext() const
        {
            auto ptx = random_plaintext();
            seal::Ciphertext ctx;
            encryptor.encrypt_symmetric(ptx, ctx);
            if (ctx.is_ntt_form())
            {
                return ctx;
            }

            evaluator.transform_to_ntt_inplace(ctx);
            return ctx;
        }

        seal::Ciphertext random_nontt_ciphertext()
        {
            auto ptx = random_plaintext();
            seal::Ciphertext ctx;
            encryptor.encrypt_symmetric(ptx, ctx);
            evaluator.transform_from_ntt_inplace(ctx);
            return ctx;
        }
    };

    std::vector<seal::Ciphertext> random_ctx_vector(const SetupObjs &all, int size)
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
    template <typename U, typename V>
    seal::Ciphertext mult_row(
        const SetupObjs &all, uint64_t i, uint64_t j, const seal::util::matrix<U> &left,
        const seal::util::matrix<V> &right)
    {
        seal::Ciphertext tmp;
        seal::Ciphertext tmp_result(all.context);
        tmp_result.is_ntt_form() = true;

        // assume that in seal::Plaintext case we don't want to turn into splitPlaintexts
        for (uint64_t k = 0; k < left.cols; ++k)
        {
            if constexpr ((std::is_same_v<V, seal::Plaintext>))
            {
                all.evaluator.multiply_plain(left(i, k), right(k, j), tmp);
            }
            else if constexpr (std::is_same_v<U, seal::Plaintext>)
            {
                all.evaluator.multiply_plain(right(k, j), left(i, k), tmp);
            }
            else
            {
                all.evaluator.multiply(left(i, k), right(k, j), tmp);
            }
            all.evaluator.add(tmp, tmp_result, tmp_result);
        }
        return tmp_result;
    }

    template <typename U, typename V>
    void mat_mult(
        const SetupObjs &all, const seal::util::matrix<U> &left, const seal::util::matrix<V> &right,
        seal::util::matrix<seal::Ciphertext> &result)
    {
        for (uint64_t i = 0; i < left.rows; ++i)
        {
            for (uint64_t j = 0; j < right.cols; ++j)
            {
                result(i, j) = mult_row(all, i, j, left, right);
            }
        }
    }

    template <typename T, typename U>
    seal::util::matrix<Ciphertext> multiplyMatrices(
        SetupObjs &all, const seal::util::matrix<T> &a, const seal::util::matrix<U> &b)
    {
        if constexpr (std::is_same_v<T, seal::Plaintext> && std::is_same_v<U, seal::Plaintext>)
        {
            throw std::invalid_argument("MatrixOperations::multiply: cannot multiply plaintext by plaintext");
        }
        verify_correct_dimension(a, b);
        verify_not_empty_matrices(a, b);

        seal::util::matrix<Ciphertext> result(a.rows, b.cols);
        mat_mult(all, a, b, result);
        return result;
    }

    vector<util::matrix<fractures::CiphertextFracture>> fracture_matrix(
        const SetupObjs &all, util::matrix<Ciphertext> &ctx_mat, std::uint64_t num_fractures)
    {
        // take each ctx and shred it using ctxShredder. then create matrices with the same dims.
        std::vector<fractures::CiphertextShredder> ctx_shredders;
        for (auto &ctx_iter : ctx_mat.data)
        {
            ctx_shredders.emplace_back(ctx_iter, all.context, num_fractures);
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

    vector<util::matrix<fractures::PolynomialFracture>> fracture_matrix(
        const SetupObjs &all, util::matrix<Plaintext> &ctx_mat, std::uint64_t num_fractures)
    {
        // take each ctx and shred it using ctxShredder. then create matrices with the same dims.
        std::vector<fractures::Polynomial> ptx_shredder;
        for (auto &ctx_iter : ctx_mat.data)
        {
            ptx_shredder.emplace_back(ctx_iter, all.context, num_fractures);
        }
        seal::util::matrix<fractures::Polynomial> shredding(ctx_mat.rows, ctx_mat.cols, ptx_shredder);

        std::vector<seal::util::matrix<fractures::PolynomialFracture>> fractured_matrices;
        for (std::uint64_t frac_num = 0; frac_num < num_fractures; ++frac_num)
        {
            std::vector<fractures::PolynomialFracture> tmp;
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

    void multiply_scalar(
        const std::vector<std::uint64_t> &a, std::vector<std::uint64_t> &dest,
        const std::vector<seal::Modulus> &modulus)
    {
        if (a.size() != dest.size())
        {
            throw std::invalid_argument("multiply_scalar: a and dest must have the same size");
        }
        if (a.size() > modulus.size())
        {
            throw std::invalid_argument("multiply_scalar: a and modulus must have the same size");
        }

        for (std::uint64_t i = 0; i < a.size(); ++i)
        {
            dest[i] = util::multiply_uint_mod(a[i], dest[i], modulus[i]);
        }
    }
    std::string rns_number_to_string(const std::vector<std::uint64_t> &rns_number)
    {
        std::string s;
        s += "{ ";
        for (auto &i : rns_number)
        {
            s += std::to_string(i);
            if (&i != &rns_number.back())
            {
                s += ", ";
            }
        }
        s += " }";
        return s;
    }

    void ctx_to_json(std::stringstream &ss, const SetupObjs &all, const seal::Ciphertext &ctx)
    {
        //        auto &ss = std::cout;
        ss << "{" << std::endl;

        auto modulus = all.context.first_context_data()->parms().coeff_modulus();
        std::uint64_t poly_number = 0;
        SEAL_ITERATE(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter rns_iter_per_poly) {
            std::uint64_t rns_number = 0;
            ss << "\""
               << "poly" << poly_number++ << "\":{";
            SEAL_ITERATE(seal::util::iter(rns_iter_per_poly), modulus.size(), [&](auto coef_iter) {
                ss << "\"" << modulus[rns_number++].value() << "\":[";
                if (poly_number == 3)
                {
                    std::cout << "";
                }
                for (uint64_t j = 0; j < all.context.first_context_data()->parms().poly_modulus_degree(); ++j)
                {
                    ss << *coef_iter << ",";
                    coef_iter++;
                }
                ss << "]," << std::endl;
            });
            ss << "}," << std::endl;
        });
        ss << "}" << std::endl;
    }

    void ctx_json_into_file(std::stringstream &ss, const std::string &filename)
    {
        //        code that dumps all contents of string stream into a file:
        std::ofstream out(filename);
        out << ss.str();
        out.close();
    }

    void ctx_into_file(const SetupObjs &all, const seal::Ciphertext &ctx, const std::string &filename)
    {
        std::stringstream ss;
        ctx_to_json(ss, all, ctx);
        ctx_json_into_file(ss, filename);
    }

} // namespace sealtest