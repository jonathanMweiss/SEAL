// Test for fractured operations by JonathanWeiss.

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

using namespace seal;
using namespace std;
namespace sealtest
{
    TEST(FracturedOps, Poly)
    {
        std::cout << "hello friend" << std::endl;
    }
} // namespace sealtest