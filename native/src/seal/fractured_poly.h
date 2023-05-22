#pragma once

#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "matrix.h"

namespace seal
{
    namespace fractures
    {

        // Represents the key parameters to fracture a polynomial.
        struct Essence
        {
            //            const seal::SEALContext::ContextData context_data;
            seal::EncryptionParameters parms;
            std::vector<seal::Modulus> coeff_modulus;
            std::uint64_t coeff_count;
            std::uint64_t coeff_modulus_size;

            explicit Essence(const seal::SEALContext &context_data, seal::EncryptionParameters params)
                : parms(std::move(params)), coeff_modulus(parms.coeff_modulus()),
                  coeff_count(parms.poly_modulus_degree()), coeff_modulus_size(coeff_modulus.size()){};
        };

        class Poly
        {
            matrix<std::uint64_t> poly_data;
            std::uint64_t num_fractures;

            // gets a polynomial iterator, that starts at the correct position in the matrix.
            seal::util::CoeffIter rns_poly_iter(std::uint64_t rns_num)
            {
                return { &poly_data(0, rns_num) };
            }

        public:
            // Creates a plaintext with the given number of coefficients and modulus size.
            static Poly crate_plaintext(
                const seal::Evaluator &ev, const seal::fractures::Essence &e, std::uint16_t num_fractures,
                const seal::Plaintext &ptx)
            {
                if (ptx.is_ntt_form())
                {
                    throw std::invalid_argument("Poly is already in NTT form");
                }

                // copying ptx, so we can transform it to NTT form.
                seal::Plaintext p(ptx);
                ev.transform_to_ntt_inplace(p, e.parms.parms_id());

                return seal::fractures::Poly(p, e, num_fractures);
            }

            explicit Poly(seal::util::ConstRNSIter rns_iter, const Essence &e, std::uint64_t num_fractures) noexcept
                : poly_data(e.coeff_count, e.coeff_modulus_size), num_fractures(num_fractures)
            {
                // 1. create a matrix of size (colums) e.coeff_count x (rows) e.coeff_modulus_size
                // 2. fill the matrix with the coefficients of the plaintext
                // Create fractures that are capable to run through the matrix ( they go through each row for N values)
                // When we need to marshal such values, we need them to be responsbile for multiple things -
                //   they need to multiply with others that have the same positions in the point-value representation.
                //   they need to be able to unmarshal themselves, along with correct "iterators".

                // consists of e.coeff_modulus_size steps, where each step is of size e.coeff_count.
                std::uint64_t rns_num = 0;
                std::for_each_n(seal::util::iter(rns_iter), e.coeff_modulus_size, [&](auto coef_iter) {
                    //                                        coef_iter += index * e.coeff_count / num_fractures;
                    auto frac_coef_iter = rns_poly_iter(rns_num);
                    rns_num++;

                    for (uint64_t j = 0; j < e.coeff_count; ++j)
                    {
                        *frac_coef_iter = *coef_iter;
                        frac_coef_iter++;
                        coef_iter++;
                    }
                });
            }

        private:
            // Expects p to be in NTT form.
            explicit Poly(seal::Plaintext &p, const Essence &e, std::uint64_t num_fractures) noexcept
                : Poly(seal::util::ConstRNSIter(p.data(), e.coeff_count), e, num_fractures){};
        };

    } // namespace fractures
} // namespace seal
