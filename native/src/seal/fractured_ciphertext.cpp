#include "fractured_ciphertext.h"

namespace seal::fractures
{

    CiphertextShredder::CiphertextShredder(seal::Ciphertext &ctx, Essence e, std::uint64_t num_fractures) noexcept
        : num_coefficients(num_coefficients)
    {
        std::vector<Polynomial> poly_shredders;
        poly_shredders.reserve(ctx.size());
        // TODO: create a copy once, and give it the positions it can use.
        std::for_each_n(seal::util::PolyIter(ctx), ctx.size(), [&](seal::util::RNSIter rns_iter_per_poly) {
            poly_shredders.emplace_back(rns_iter_per_poly, e, num_fractures);
        });

        ctx_parts.reserve(num_fractures);

        for (std::uint64_t i = 0; i < num_fractures; ++i)
        {
            CiphertextFracture fracture{ {}, i, e.coeff_modulus };

            fracture.poly_fracs.reserve(ctx.size());

            for (std::uint64_t j = 0; j < ctx.size(); ++j)
            {
                fracture.poly_fracs.push_back(poly_shredders[j].get_fracture(i));
            }

            ctx_parts.push_back(fracture);
        }
    }

    CiphertextFracture CiphertextShredder::get_fracture(std::uint64_t index)
    {
        if (index >= ctx_parts.size())
        {
            throw std::invalid_argument("index out of range");
        }
        return ctx_parts[index];
    }

    void CiphertextShredder::into_ciphertext(seal::Ciphertext &ctx)
    {
        auto num_fractures = ctx_parts.size();

        auto nelements = 0;
        auto poly_iter_num = -1;
        std::for_each_n(seal::util::PolyIter(ctx), ctx.size(), [&](seal::util::RNSIter p_i) {
            poly_iter_num++;

            auto poly_modulus_degree = p_i.poly_modulus_degree();
            auto fracture_size = poly_modulus_degree / num_fractures;
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                auto rns_iter_num = -1;
                // because of RNS we perform this for-each. (coeff_modulus_size elements per "coef")
                std::for_each_n(iter(p_i), ctx_parts[0].coeff_modulus.size(), [&](seal::util::CoeffIter write_to) {
                    rns_iter_num++;

                    write_to += i * fracture_size;
                    auto read_from = ctx_parts[i].poly_fracs[poly_iter_num].const_rns_poly_iter(rns_iter_num);
                    for (std::uint64_t j = 0; j < fracture_size; ++j)
                    {
                        *write_to = *read_from;
                        write_to++;
                        read_from++;
                        nelements++;
                    }
                });
            }
        });
    }

    const CiphertextFracture &CiphertextFracture::operator*=(const PolynomialFracture &poly)
    {
        if (index != poly.fracture_index)
        {
            throw std::invalid_argument("fracture indices do not match");
        }

        for (auto &frac : poly_fracs)
        {
            for (std::uint64_t rns_num = 0; rns_num < coeff_modulus.size(); ++rns_num)
            {
                seal::util::dyadic_product_coeffmod(
                    frac.const_rns_poly_iter(rns_num), poly.const_rns_poly_iter(rns_num), frac.coeff_count,
                    coeff_modulus[rns_num], frac.rns_poly_iter(rns_num));
            }
        }

        return *this;
    }
} // namespace seal::fractures