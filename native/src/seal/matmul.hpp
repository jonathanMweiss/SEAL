#pragma once
#include "seal/batchencoder.h"
#include "seal/ckks.h"
#include "seal/context.h"
#include "seal/decryptor.h"
#include "seal/encryptor.h"
#include "seal/evaluator.h"
#include "seal/fractures.h"
#include "seal/keygenerator.h"
#include "seal/matrix.h"
#include "seal/modulus.h"
#include "seal/poly_eval.h"
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <string>
#include <utility>

template <typename T, typename U>
void verify_not_empty_matrices(const seal::util::matrix<T> &left, const seal::util::matrix<U> &right)
{
    if (left.data.empty())
    {
        throw std::invalid_argument("MatrixOperations::multiply: received empty left matrix");
    }
    if (right.data.empty())
    {
        throw std::invalid_argument("MatrixOperations::multiply: received empty right matrix");
    }
}

template <typename T, typename U>
void verify_correct_dimension(const seal::util::matrix<T> &left, const seal::util::matrix<U> &right)
{
    if (left.cols != right.rows)
    {
        throw std::invalid_argument("MatrixOperations::multiply: left matrix cols != right matrix rows");
    }
}

// Doing matmul...
template <typename U, typename T>
seal::util::matrix<seal::fractures::CiphertextFracture> multiply_matrices(
    const seal::util::matrix<U> &a, const seal::util::matrix<T> &b)
{
    if constexpr ((std::is_same_v<T, seal::Plaintext>)&&std::is_same_v<U, seal::Plaintext>)
    {
        throw std::invalid_argument("MatrixOperations::multiply: cannot multiply ptx by ptx");
    }
    verify_correct_dimension(a, b);

    seal::util::matrix<seal::fractures::CiphertextFracture> result(a.rows, b.cols);

    std::uint64_t poly_sz = 2;

    if constexpr ((std::is_same_v<T, seal::fractures::CiphertextFracture>)&&std::is_same_v<
                      U, seal::fractures::CiphertextFracture>)
    {
        poly_sz = 3;
    }

    const seal::fractures::CiphertextFracture *ctx = nullptr;
    if constexpr (std::is_same_v<T, seal::fractures::CiphertextFracture>)
    {
        ctx = &b(0, 0);
    }
    else
    {
        ctx = &a(0, 0);
    }

    if (ctx == nullptr)
    {
        throw std::invalid_argument("MatrixOperations::multiply: cannot multiply ptx by ptx");
    }

    auto index = ctx->index;
    auto coeff_modulus = ctx->coeff_modulus;
    auto coeff_count = ctx->poly_fracs[0].coeff_count;

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
