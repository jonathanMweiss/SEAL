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
        //        auto all = SetupObjs::New();
        std::uint64_t num_fractures = 256;

        //        seal::Ciphertext encrypted_ntt = all.random_ciphertext();
        //        auto ptx = all.random_plaintext();
        //        all.evaluator.transform_to_ntt_inplace(ptx, encrypted_ntt.parms_id());

        //        seal::Ciphertext encrypted_ntt_cpy(encrypted_ntt);
        //        seal::Plaintext ptx_cpy(ptx);
        seal::fractures::Essence essence(bm_env->context(), bm_env->parms());
        seal::Plaintext ptx;
        bm_env->randomize_pt_bgv(ptx);

        for (auto _ : state)
        {
            seal::fractures::Polynomial pshredder(ptx, essence, num_fractures);
        }
    }
} // namespace sealbench