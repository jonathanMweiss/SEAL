#pragma once

#include "seal/fractures.h"
#include <seal/util/polyarithsmallmod.h>

namespace seal::fractures
{
    struct CiphertextFracture
    {
        std::vector<PolynomialFracture> poly_fracs;
        std::uint64_t index;
        std::vector<seal::Modulus> coeff_modulus;

        static CiphertextFracture Empty(
            std::uint64_t num_polys, std::uint64_t index, std::uint64_t num_coefs, std::vector<seal::Modulus> coef_mod);

        static CiphertextFracture Empty(const CiphertextFracture &ctxf);

        const CiphertextFracture &operator*=(const PolynomialFracture &poly);
        const CiphertextFracture &operator+=(const CiphertextFracture &ctxf);
        const CiphertextFracture &operator*=(const CiphertextFracture &y);

        bool operator==(const CiphertextFracture &ctxf);

        // TODO: verify.
        CiphertextFracture operator*(const PolynomialFracture &poly) const;
        CiphertextFracture operator*(const CiphertextFracture &y) const;
        CiphertextFracture operator+(const CiphertextFracture &ctxf);

        /** Returns the size of the object in bytes, as if written to a stream. **/
        std::streamoff save_size(compr_mode_type compr_mode) const;

        inline std::streamoff save(
            std::ostream &stream, compr_mode_type compr_mode = Serialization::compr_mode_default) const
        {
            using namespace std::placeholders;
            return Serialization::Save(
                std::bind(&CiphertextFracture::save_members, this, _1), save_size(compr_mode_type::none), stream,
                compr_mode, false);
        }

        inline std::streamoff load(const seal_byte *in, std::size_t size)
        {
            using namespace std::placeholders;
            return Serialization::Load(std::bind(&CiphertextFracture::load_members, this, _1, _2), in, size, false);
        }

        inline std::streamoff load(std::istream &in)
        {
            using namespace std::placeholders;
            return Serialization::Load(std::bind(&CiphertextFracture::load_members, this, _1, _2), in, false);
        }

    private:
        // Saves the (uncompressed) object to an output stream.
        void save_members(std::ostream &stream) const;
        void load_members(std::istream &stream, SEALVersion version);
    };

    class CiphertextFracturingTool
    {
    public:
        static CiphertextFracture compute_fracture(
            const seal::Ciphertext &ctx, const seal::SEALContext &context, uint64_t num_fractures, uint64_t index)
        {
            auto context_data = context.get_context_data(ctx.parms_id());
            auto &parms = context_data->parms();
            auto &coeff_modulus = parms.coeff_modulus();

            CiphertextFracture fracture{ {}, index, coeff_modulus };
            fracture.poly_fracs.reserve(ctx.size());
            
            SEAL_ITERATE(seal::util::ConstPolyIter(ctx), ctx.size(), [&](seal::util::ConstRNSIter rns_iter_per_poly) {
                // iterating through the ctx.size() number of polynomials.
                fracture.poly_fracs.emplace_back(PolynomialFracturingTool::compute_fracture(
                    num_fractures, index, parms.poly_modulus_degree(), coeff_modulus.size(), rns_iter_per_poly));
            });
            return fracture;
        }
    };

    class CiphertextBuilderTool
    {
    public:
        static void validate_parms_id(
            const std::vector<CiphertextFracture> &ctx_parts, const SEALContext &context, const parms_id_type &parms_id)
        {
            auto context_data = context.get_context_data(parms_id);
            auto &parms = context_data->parms();
            auto &coeff_modulus = parms.coeff_modulus();

            if (ctx_parts[0].coeff_modulus.size() != coeff_modulus.size())
            {
                throw std::invalid_argument("validating parms_id failed in CiphertextBuilderTool, coeff_modulus size "
                                            "of ctx fractures mismatch parms coeff_modulus size");
            }
        }

        static seal::Ciphertext into_ciphertext(
            const std::vector<CiphertextFracture> &ctx_parts, const seal::SEALContext &context)
        {
            return into_ciphertext(ctx_parts, context, context.first_parms_id());
        }

        /**
         * This function will build a full ciphertext from the fractured parts.
         * @return
         **/
        static seal::Ciphertext into_ciphertext(
            const std::vector<CiphertextFracture> &ctx_parts, const seal::SEALContext &context,
            const parms_id_type parms_id)
        {
            validate_parms_id(ctx_parts, context, parms_id);

            auto context_data = context.get_context_data(parms_id);
            auto &parms = context_data->parms();
            auto &coeff_modulus = parms.coeff_modulus();

            seal::Ciphertext ctx(context, parms_id);
            ctx.resize(ctx_parts[0].poly_fracs.size());
            ctx.is_ntt_form() = true;

            auto num_fractures = ctx_parts.size();
            auto poly_iter_num = -1;
            SEAL_ITERATE(seal::util::PolyIter(ctx), ctx.size(), [&](seal::util::RNSIter p_i) {
                poly_iter_num++;

                auto poly_modulus_degree = p_i.poly_modulus_degree();
                auto fracture_size = poly_modulus_degree / num_fractures;
                for (std::uint64_t i = 0; i < num_fractures; ++i)
                {
                    auto rns_iter_num = -1;
                    SEAL_ITERATE(iter(p_i), coeff_modulus.size(), [&](seal::util::CoeffIter write_to) {
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
    };

    /**
     * This class is responsible for fracturing ciphertext into multiple parts.
     * Each fractured part should be able to be multiplied with a fractured plainetext.
     * Using fractured ciphertexts and plaintexts, we can perform multiplication in parallel on multiple machines.
     */
    class CiphertextShredder
    {
    public:
        CiphertextShredder(const seal::SEALContext &, std::uint64_t num_fractures) noexcept;

        explicit CiphertextShredder(
            const seal::Ciphertext &ctx, const SEALContext &context, std::uint64_t num_fractures) noexcept;

        CiphertextFracture &operator[](std::uint64_t index)
        {
            return ctx_parts[index];
        }

        // todo: func that turns this shredder back to full ciphertext.
        seal::Ciphertext into_ciphertext() const;

        uint64_t num_fractures();

        void set_fracture(uint64_t i, CiphertextFracture fracture);

    private:
        seal::SEALContext context;
        std::vector<CiphertextFracture> ctx_parts;
    };
} // namespace seal::fractures