#pragma once

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

    shared_ptr<UniformRandomGenerator> prng();

    seal::Plaintext random_plain(const seal::EncryptionParameters &enc_params);
    void find_possible_parameters();
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
            enc.set_coeff_modulus(seal::CoeffModulus::BFVDefault(N, sec_level_type::tc192));

//            auto tmp = std::vector<int>{ 54, 41, 42 };
            auto tmp = std::vector<int>{ 51, 57, 57 };
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
            : enc_params(std::move(encryption_params)), context(enc_params, true, sec_level_type::tc192, 2),
              evaluator(context), keygen(context), secret_key(keygen.secret_key()), encryptor(context, secret_key),
              decryptor(context, secret_key){};

        seal::Plaintext random_plaintext() const
        {
            return random_plain(enc_params);
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

    std::vector<seal::Ciphertext> random_ctx_vector(const SetupObjs &all, int size);
    std::vector<seal::Plaintext> random_ptxs(const SetupObjs &all, int size);

    seal::Ciphertext ctxWithSize3(const SetupObjs &all);

    void add_inplace_ctx_fractures(
        seal::fractures::CiphertextShredder &ctxshred1, seal::fractures::CiphertextShredder &ctxshred2);
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
        const SetupObjs &all, util::matrix<Ciphertext> &ctx_mat, std::uint64_t num_fractures);

    vector<util::matrix<fractures::PolynomialFracture>> fracture_matrix(
        const SetupObjs &all, util::matrix<Plaintext> &ctx_mat, std::uint64_t num_fractures);

    std::string rns_number_to_string(const std::vector<std::uint64_t> &rns_number);

    void ctx_to_json(std::stringstream &ss, const SetupObjs &all, const seal::Ciphertext &ctx);

    void ctx_json_into_file(std::stringstream &ss, const std::string &filename);

    void ctx_into_file(const SetupObjs &all, const seal::Ciphertext &ctx, const std::string &filename);

    // root of X^n+1 in RNS form. for specific ring of BGV used in all of these tests.
    const std::vector<std::uint64_t> root{ 9354911369072846, 1245024710537 };

    template <typename MatrixType, typename EvaluatedType>
    util::matrix<EvaluatedType> evaluate_matrix(const SetupObjs &all, seal::util::matrix<MatrixType> &m)
    {
        std::vector<EvaluatedType> inner_vec;
        seal::fractures::PolynomialEvaluator pe(all.context);

        for (std::size_t i = 0; i < m.rows; ++i)
        {
            for (std::size_t j = 0; j < m.cols; ++j)
            {
                inner_vec.emplace_back(pe.evaluate(m(i, j), root));
            }
        }
        return seal::util::matrix<EvaluatedType>(m.rows, m.cols, std::move(inner_vec));
    }

    bool random_ctx_ctx_sz(const SetupObjs &all, const vector<std::uint64_t> &r);
    bool random_ptx_ctx_sz(const SetupObjs &all, const vector<std::uint64_t> &r);

    template <typename T>
    void apply_on_each_element(seal::util::matrix<T> &m, std::function<void(T &)> f)
    {
        for (auto &elem : m.data)
        {
            f(elem);
        }
    }

    template <typename T>
    void apply_on_each_element(seal::util::matrix<T> &m, void (*f)(T &))
    {
        for (auto &elem : m.data)
        {
            f(elem);
        }
    }
    std::vector<std::uint64_t> dynarray_to_vector(const seal::DynArray<std::uint64_t> &dyn);
    std::vector<std::uint64_t> plain_to_vector(const Plaintext &ptx);

} // namespace sealtest