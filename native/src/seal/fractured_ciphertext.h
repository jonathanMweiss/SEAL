#pragma once

#include "seal/fractured_ciphertext.h"
#include "seal/fractured_polynomial.h"
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

        bool operator==(const CiphertextFracture &ctxf) const;
        bool operator!=(const CiphertextFracture &ctxf) const;

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
            const seal::Ciphertext &ctx, const seal::SEALContext &context, uint64_t num_fractures, uint64_t index);
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
            const std::vector<CiphertextFracture> &ctx_parts, const SEALContext &context,
            const parms_id_type &parms_id);
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

//        inline seal::Ciphertext into_ciphertext() const
//        {
//            return into_ciphertext(context.first_parms_id());
//        }

        seal::Ciphertext into_ciphertext(const parms_id_type &parms_id) const;

        uint64_t num_fractures();

        void set_fracture(uint64_t i, CiphertextFracture fracture);

    private:
        seal::SEALContext context;
        std::vector<CiphertextFracture> ctx_parts;
    };
} // namespace seal::fractures