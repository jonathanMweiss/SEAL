#pragma once

#include <seal/util/polyarithsmallmod.h>
#include "seal/fractures.h"

namespace seal::fractures
{
    struct CiphertextFracture
    {
        static CiphertextFracture Empty(
            std::uint64_t num_polys, std::uint64_t index, std::uint64_t num_coefs, std::vector<seal::Modulus> coef_mod);

        static CiphertextFracture Empty(const CiphertextFracture &ctxf);

        std::vector<PolynomialFracture> poly_fracs;
        std::uint64_t index;
        std::vector<seal::Modulus> coeff_modulus;

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

    /**
     * This class is responsible for fracturing ciphertext into multiple parts.
     * Each fractured part should be able to be multiplied with a fractured plainetext.
     * Using fractured ciphertexts and plaintexts, we can perform multiplication in parallel on multiple machines.
     */
    class CiphertextShredder
    {
    public:
        CiphertextShredder(Essence e, std::uint64_t num_fractures) noexcept;

        explicit CiphertextShredder(const seal::Ciphertext &ctx, Essence e, std::uint64_t num_fractures) noexcept;

        CiphertextFracture &operator[](std::uint64_t index)
        {
            return ctx_parts[index];
        }

        // todo: func that turns this shredder back to full ciphertext.
        seal::Ciphertext into_ciphertext() const;

        uint64_t num_fractures();

        void set_fracture(uint64_t i, CiphertextFracture fracture);

    private:
        seal::fractures::Essence essence;
        std::vector<CiphertextFracture> ctx_parts;
    };
} // namespace seal::fractures