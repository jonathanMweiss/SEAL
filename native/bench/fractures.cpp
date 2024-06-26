#include "seal/seal.h"
#include "seal/util/rlwe.h"
#include "bench.h"

using namespace benchmark;
using namespace sealbench;
using namespace seal;
using namespace std;

namespace sealbench
{

    void polynomial_fracture(State &state, shared_ptr<BMEnv> bm_env)
    {
        std::uint64_t num_fractures = 256;
        seal::fractures::Essence essence(bm_env->context(), bm_env->parms());
        seal::Plaintext ptx;
        bm_env->randomize_pt_bgv(ptx);

        for (auto _ : state)
        {
            seal::fractures::Polynomial pshredder(ptx, essence, num_fractures);
        }
    }

    void fracture_ctx(State &state, shared_ptr<BMEnv> bm_env)
    {
        std::uint64_t num_fractures = 256;
        seal::fractures::Essence essence(bm_env->context(), bm_env->parms());

        seal::Ciphertext ctx;
        bm_env->randomize_ct_bgv(ctx);
        for (auto _ : state)
        {
            seal::fractures::CiphertextShredder cshred(ctx, essence, num_fractures);
        }
    }

    void fractured_plain_mult(State &state, shared_ptr<BMEnv> bm_env)
    {
        std::uint64_t num_fractures = 256;
        seal::fractures::Essence essence(bm_env->context(), bm_env->parms());

        seal::Ciphertext ctx;
        bm_env->randomize_ct_bgv(ctx);

        seal::Plaintext ptx;
        bm_env->randomize_pt_bgv(ptx);

        seal::fractures::Polynomial pshredder(ptx, essence, num_fractures);
        seal::fractures::CiphertextShredder cshredder(ctx, essence, num_fractures);

        for (auto _ : state)
        {
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                cshredder[i] *= pshredder.get_fracture(i);
            }
        }
    }

    void fractured_ctx_mult(State &state, shared_ptr<BMEnv> bm_env)
    {
        std::uint64_t num_fractures = 256;
        seal::fractures::Essence essence(bm_env->context(), bm_env->parms());

        seal::Ciphertext ctx1;
        bm_env->randomize_ct_bgv(ctx1);

        seal::Ciphertext ctx2;
        bm_env->randomize_ct_bgv(ctx2);

        for (auto _ : state)
        {
            state.PauseTiming();
            seal::fractures::CiphertextShredder cshredder1(ctx1, essence, num_fractures);
            seal::fractures::CiphertextShredder cshredder2(ctx2, essence, num_fractures);
            state.ResumeTiming();

            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                cshredder1[i] *= cshredder2[i];
            }
        }
    }
} // namespace sealbench