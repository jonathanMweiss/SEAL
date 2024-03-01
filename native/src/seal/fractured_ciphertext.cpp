#include "fractured_ciphertext.h"
#include "defines.h"
#include <utility>

namespace seal::fractures
{
    CiphertextShredder::CiphertextShredder(Essence e, std::uint64_t num_fractures) noexcept
        : essence(std::move(e)), ctx_parts(num_fractures)
    {}

    CiphertextShredder::CiphertextShredder(const seal::Ciphertext &ctx, Essence e, std::uint64_t num_fractures) noexcept
        : essence(e)
    {
        std::vector<Polynomial> poly_shredders;
        poly_shredders.reserve(ctx.size());

        SEAL_ITERATE(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter rns_iter_per_poly) {
            poly_shredders.emplace_back(rns_iter_per_poly, e, num_fractures);
        });

        ctx_parts.reserve(num_fractures);

        for (std::uint64_t i = 0; i < num_fractures; ++i)
        {
            CiphertextFracture fracture{ {}, i, essence.parms.coeff_modulus() };

            fracture.poly_fracs.reserve(ctx.size());

            for (std::uint64_t j = 0; j < ctx.size(); ++j)
            {
                fracture.poly_fracs.push_back(poly_shredders[j].get_fracture(i));
            }

            ctx_parts.push_back(fracture);
        }
    }

    seal::Ciphertext CiphertextShredder::into_ciphertext() const
    {
        auto num_fractures = ctx_parts.size();

        seal::Ciphertext ctx(essence.parms);
        ctx.resize(ctx_parts[0].poly_fracs.size());
        ctx.is_ntt_form() = true;

        auto poly_iter_num = -1;
        SEAL_ITERATE(seal::util::PolyIter(ctx), ctx.size(), [&](seal::util::RNSIter p_i) {
            poly_iter_num++;

            auto poly_modulus_degree = p_i.poly_modulus_degree();
            auto fracture_size = poly_modulus_degree / num_fractures;
            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                auto rns_iter_num = -1;
                // because of RNS we perform this for-each. (coeff_modulus_size elements per "coef")
                // We go up to -1 because the inner loop will go through the last element
                SEAL_ITERATE(iter(p_i), ctx_parts[0].coeff_modulus.size() - 1, [&](seal::util::CoeffIter write_to) {
                    rns_iter_num++;

                    write_to += i * fracture_size;
                    auto read_from = ctx_parts[i].poly_fracs[poly_iter_num].const_rns_poly_iter(rns_iter_num);
                    for (std::uint64_t j = 0; j < fracture_size; ++j)
                    {
                        *write_to = *read_from;
                        write_to++;
                        read_from++;
                    }
                });
            }
        });

        return ctx;
    }

    uint64_t CiphertextShredder::num_fractures()
    {
        return ctx_parts.size();
    }

    void CiphertextShredder::set_fracture(uint64_t i, CiphertextFracture fracture)
    {
        ctx_parts[i] = std::move(fracture);
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

    const CiphertextFracture &CiphertextFracture::operator+=(const CiphertextFracture &ctxf)
    {
        if (index != ctxf.index)
        {
            throw std::invalid_argument("fracture indices do not match");
        }
        if (ctxf.poly_fracs.size() != poly_fracs.size())
        {
            throw std::invalid_argument("fracture sizes do not match");
        }

        for (std::uint64_t i = 0; i < poly_fracs.size(); ++i)
        {
            for (std::uint64_t rns_num = 0; rns_num < coeff_modulus.size(); ++rns_num)
            {
#ifdef SEAL_DEBUG
                if (poly_fracs[i].fracture_index != ctxf.poly_fracs[i].fracture_index)
                {
                    throw std::invalid_argument("fracture indices do not match");
                }
#endif

                seal::util::add_poly_coeffmod(
                    poly_fracs[i].const_rns_poly_iter(rns_num), ctxf.poly_fracs[i].const_rns_poly_iter(rns_num),
                    poly_fracs[i].coeff_count, coeff_modulus[rns_num], poly_fracs[i].rns_poly_iter(rns_num));
            }
        }

        return *this;
    }

    const CiphertextFracture &CiphertextFracture::operator*=(const CiphertextFracture &ctxf)
    {
        if (index != ctxf.index)
        {
            throw std::invalid_argument("fracture indices do not match");
        }
        if (ctxf.poly_fracs.size() != 2)
        {
            throw std::invalid_argument("not supported");
        }
        if (ctxf.poly_fracs.size() != poly_fracs.size())
        {
            throw std::invalid_argument("fracture sizes do not match");
        }

        this->poly_fracs.emplace_back(index, ctxf.poly_fracs[0].coeff_count, coeff_modulus.size());
        PolynomialFracture tmp(index, ctxf.poly_fracs[0].coeff_count, coeff_modulus.size());
        for (std::uint64_t rns_num = 0; rns_num < coeff_modulus.size(); ++rns_num)
        {
            // need to do some ops between the polynomials...

            // x[2] = x[1] * y[1]
            seal::util::dyadic_product_coeffmod(
                this->poly_fracs[1].const_rns_poly_iter(rns_num), ctxf.poly_fracs[1].const_rns_poly_iter(rns_num),
                this->poly_fracs[2].coeff_count, coeff_modulus[rns_num], this->poly_fracs[2].rns_poly_iter(rns_num));

            // Compute second output polynomial, overwriting input
            // temp = x[1] * y[0]
            seal::util::dyadic_product_coeffmod(
                this->poly_fracs[1].const_rns_poly_iter(rns_num), ctxf.poly_fracs[0].const_rns_poly_iter(rns_num),
                tmp.coeff_count, coeff_modulus[rns_num], tmp.rns_poly_iter(rns_num));
            // x[1] = x[0] * y[1]
            seal::util::dyadic_product_coeffmod(
                this->poly_fracs[0].const_rns_poly_iter(rns_num), ctxf.poly_fracs[1].const_rns_poly_iter(rns_num),
                this->poly_fracs[1].coeff_count, coeff_modulus[rns_num], this->poly_fracs[1].rns_poly_iter(rns_num));
            // x[1] += temp
            seal::util::add_poly_coeffmod(
                this->poly_fracs[1].const_rns_poly_iter(rns_num), tmp.const_rns_poly_iter(rns_num), tmp.coeff_count,
                coeff_modulus[rns_num], this->poly_fracs[1].rns_poly_iter(rns_num));

            // Compute first output polynomial, overwriting input
            // x[0] = x[0] * y[0]
            seal::util::dyadic_product_coeffmod(
                this->poly_fracs[0].const_rns_poly_iter(rns_num), ctxf.poly_fracs[0].const_rns_poly_iter(rns_num),
                this->poly_fracs[0].coeff_count, coeff_modulus[rns_num], this->poly_fracs[0].rns_poly_iter(rns_num));
        }

        return *this;
    }

    CiphertextFracture CiphertextFracture::operator*(const CiphertextFracture &y) const
    {
        CiphertextFracture cpy(*this);
        cpy *= y;
        return cpy;
    }

    CiphertextFracture CiphertextFracture::operator+(const CiphertextFracture &y)
    {
        CiphertextFracture cpy(*this);
        cpy += y;
        return cpy;
    }

    CiphertextFracture CiphertextFracture::operator*(const PolynomialFracture &poly) const
    {
        CiphertextFracture cpy(*this);
        cpy *= poly;
        return cpy;
    }

    CiphertextFracture CiphertextFracture::Empty(
        std::uint64_t num_polys, std::uint64_t index, std::uint64_t num_coefs, std::vector<seal::Modulus> coef_mod)
    {
        CiphertextFracture tmp{ {}, index, coef_mod };
        for (std::uint64_t i = 0; i < num_polys; ++i)
        {
            tmp.poly_fracs.emplace_back(index, num_coefs, coef_mod.size());
        }
        return tmp;
    }

    CiphertextFracture CiphertextFracture::Empty(const CiphertextFracture &ctxf)
    {
        return Empty(ctxf.poly_fracs.size(), ctxf.index, ctxf.poly_fracs[0].coeff_count, ctxf.coeff_modulus);
    }

    bool CiphertextFracture::operator==(const CiphertextFracture &ctxf)
    {
        if (index != ctxf.index)
        {
            return false;
        }

        if (poly_fracs.size() != ctxf.poly_fracs.size())
        {
            return false;
        }

        for (std::uint64_t i = 0; i < poly_fracs.size(); ++i)
        {
            if (poly_fracs[i] == ctxf.poly_fracs[i])
            {
                continue;
            }

            return false;
        }

        return true;
    }

    // expects the module to be sent some other way.
    std::streamoff CiphertextFracture::save_size(compr_mode_type compr_mode) const
    {
        auto members_size = util::add_safe(
            sizeof(std::uint64_t), // index (the index value of this fracture).
            sizeof(std::uint64_t) // num polynomial fractions (can be 2-3).
        );

        // then each fracture and its size.
        for (auto &pf : poly_fracs)
        {
            auto pf_size = pf.save_size(compr_mode);
            members_size = util::add_safe(
                static_cast<std::size_t>(members_size), // adding to oneself.
                static_cast<std::size_t>(sizeof(std::uint64_t)), // setting size for frac size
                static_cast<std::size_t>(pf_size) // size the frac needs.
            );
        }

        return util::safe_cast<std::streamoff>(Serialization::ComprSizeEstimate(
            util::add_safe(
                members_size, // adding to oneself.
                sizeof(Serialization::SEALHeader)),
            compr_mode));
    }

    void CiphertextFracture::save_members(std::ostream &stream) const
    {
        auto old_except_mask = stream.exceptions();
        auto num_fracs = poly_fracs.size();
        try
        {
            stream.write(reinterpret_cast<const char *>(&index), sizeof(std::uint64_t));
            stream.write(reinterpret_cast<const char *>(&num_fracs), sizeof(std::uint64_t));

            for (auto &pf : poly_fracs)
            {
                //                auto save_size = pf.save_size(compr_mode_type::none);
                //                stream.write(reinterpret_cast<const char *>(&save_size), sizeof(std::uint64_t));
                pf.save(stream);
            }
        }
        catch (const std::ios_base::failure &)
        {
            stream.exceptions(old_except_mask);
            throw std::runtime_error("I/O error");
        }
        catch (...)
        {
            stream.exceptions(old_except_mask);
            throw;
        }

        stream.exceptions(old_except_mask);
    }
    void CiphertextFracture::load_members(std::istream &stream, SEAL_MAYBE_UNUSED SEALVersion version)
    {
        auto old_except_mask = stream.exceptions();
        try
        {
            stream.read(reinterpret_cast<char *>(&index), sizeof(std::uint64_t));
            std::uint64_t num_fracs;
            stream.read(reinterpret_cast<char *>(&num_fracs), sizeof(std::uint64_t));

            poly_fracs.reserve(num_fracs);
            for (std::uint64_t i = 0; i < num_fracs; ++i)
            {
                //                std::uint64_t pf_size;
                //                stream.read(reinterpret_cast<char *>(&pf_size), sizeof(std::uint64_t));

                PolynomialFracture pf(index, 0, coeff_modulus.size());
                pf.load(stream);
                //                std::vector<seal::seal_byte> bts(pf_size);
                //                stream.read(reinterpret_cast<char *>(&bts[0]), static_cast<long>(pf_size));
                //                pf.load(reinterpret_cast<seal_byte *>(&bts[0]), pf_size);
                poly_fracs[i] = std::move(pf);
            }
        }
        catch (const std::ios_base::failure &)
        {
            stream.exceptions(old_except_mask);
            throw std::runtime_error("I/O error");
        }
        catch (...)
        {
            stream.exceptions(old_except_mask);
            throw;
        }
        stream.exceptions(old_except_mask);
    }
} // namespace seal::fractures