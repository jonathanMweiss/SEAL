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
#include "helpers.h"
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
            ss << "\"" << "poly" << poly_number++ << "\":{";
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

    bool random_ctx_ctx_sz(const SetupObjs &all, const vector<std::uint64_t> &r)
    {
        auto ctx1 = all.random_ciphertext();
        auto ctx2 = all.random_ciphertext();

        seal::Ciphertext res;
        all.evaluator.multiply(ctx1, ctx2, res);

        all.evaluator.transform_from_ntt_inplace(res);
        all.evaluator.transform_from_ntt_inplace(ctx1);
        all.evaluator.transform_from_ntt_inplace(ctx2);

        seal::fractures::PolynomialEvaluator pe(all.context);

        auto ev1 = pe.evaluate(ctx1, r);
        auto ev2 = pe.evaluate(ctx2, r);
        auto ev_res = pe.evaluate(res, r);

        // evaluate the multiplication:
        auto ev_mul = ev1 * ev2;

        return (ev_res == ev_mul);
    }
    bool random_ptx_ctx_sz(const SetupObjs &all, const vector<std::uint64_t> &r)
    {
        auto ptx = all.random_plaintext();
        auto ctx = all.random_ciphertext();

        seal::Plaintext cpy;
        all.evaluator.transform_to_ntt(ptx, ctx.parms_id(), cpy);

        seal::Ciphertext res;
        all.evaluator.multiply_plain(ctx, cpy, res);

        all.evaluator.transform_from_ntt_inplace(res);
        all.evaluator.transform_from_ntt_inplace(ctx);

        // ptxroot:
        seal::fractures::PolynomialEvaluator pe(all.context);
        auto ev1 = pe.evaluate(ctx, r);
        auto ev2 = pe.evaluate(ptx, r);
        auto ev_res = pe.evaluate(res, r);

        // evaluate the multiplication:
        auto ev_mul = ev1 * ev2;
        return (ev_res == ev_mul);
    }

    std::vector<std::uint64_t> dynarray_to_vector(const seal::DynArray<std::uint64_t> &dyn)
    {
        std::vector<std::uint64_t> v;
        v.reserve(dyn.size()); // x4 is for ntt pad. x2 is for modulus size
        for (std::uint64_t i = 0; i < dyn.size(); ++i)
        {
            v.emplace_back(dyn.at(i));
        }
        return v;
    }

    std::vector<std::uint64_t> plain_to_vector(const Plaintext &ptx)
    {
        return dynarray_to_vector(ptx.dyn_array());
    }

    shared_ptr<UniformRandomGenerator> prng()
    {
        Blake2xbPRNGFactory factory;
        array<uint64_t, prng_seed_uint64_count> seed{ 1, 2, 3, 4, 5, 6, 7, 8 };
        auto gen = factory.create(seed);
        return gen;
    }

    seal::Plaintext random_plain(const EncryptionParameters &enc_params)
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

    std::uint64_t vector_gcd(std::vector<std::uint64_t> v)
    {
        std::uint64_t _gcd = std::gcd(v[0], v[1]);
        for (std::uint64_t i = 2; i < v.size(); ++i)
        {
            _gcd = gcd(_gcd, v[i]);
        }

        return _gcd;
    }

    bool has_proper_roots(uint64_t wanted_ntt_deg, const vector<Modulus> &o)
    {
        vector<uint64_t> tmp;
        // NTT L(n).
        for (uint64_t l = 0; l < o.size(); ++l)
        {
            tmp.push_back((o[l].value() - 1));
        }

        // to get NTT to work (need to have roots of unity) we need the following to be true:
        return 0 == (vector_gcd(tmp) % wanted_ntt_deg);
    }

    bool is_valid_ntt_params_for_deg(int i, int j, int k, std::uint64_t N, std::uint64_t wanted_ntt_deg)
    {
        if (i + j + k < 152)
        {
            return false;
        }
        auto o = seal::CoeffModulus::Create(N, std::vector<int>{ i, j, k });
        return has_proper_roots(wanted_ntt_deg, o);
    }

    //        the following works i: 41 j: 57 k: 57
    //        the following works i: 51 j: 57 k: 57
    //        the following works i: 57 j: 41 k: 57
    //        the following works i: 57 j: 51 k: 57
    //        the following works i: 57 j: 57 k: 41
    //        the following works i: 57 j: 57 k: 51
    void find_possible_parameters()
    {
        std::uint64_t N = 8192;
        seal::EncryptionParameters parms(seal::scheme_type::bgv);
        parms.set_poly_modulus_degree(N);

        parms.set_plain_modulus(seal::PlainModulus::Batching(N, 16 + 1));

        for (int i = 41; i < 60; ++i)
        {
            for (int j = 41; j < 60; j++)
            {
                for (int k = 41; k < 60; k++)
                {
                    if (!is_valid_ntt_params_for_deg(i, j, k, N, 1 << 16))
                    {
                        continue;
                    };
                    auto tmp = std::vector<int>{ i, j, k };
                    auto o = seal::CoeffModulus::Create(N, tmp);

                    parms.set_coeff_modulus(o);
                    auto cntx = SEALContext(parms, false, sec_level_type::tc128, 2);
                    if (!cntx.parameters_set())
                    {
                        continue;
                    }

                    std::cout << "the following works i: " << i << " j: " << j << " k: " << k << std::endl;
                }
            }
        }
    }
} // namespace sealtest