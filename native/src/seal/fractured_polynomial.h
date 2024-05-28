#pragma once

#include "seal/evaluator.h"
#include "seal/plaintext.h"
#include <utility>
#include "matrix.h"

namespace seal::fractures
{

    struct CiphertextFracture;
    /**
     * A class that represents a part of a polynomial that had been split.
     */
    struct PolynomialFracture
    {
        std::uint64_t fracture_index;
        std::uint64_t coeff_count;
        seal::util::matrix<std::uint64_t> rns_coefficients;

        // An empty polynomial. used to prepare in advance inside a vector for instance.
        PolynomialFracture() : fracture_index(0), coeff_count(0), rns_coefficients({ 0, 0, {} }){}

        PolynomialFracture(std::uint64_t index, std::uint64_t _coeff_count, std::uint64_t modulus_size)
            : fracture_index(index), coeff_count(_coeff_count), rns_coefficients(coeff_count, modulus_size)
        {}

        inline seal::util::CoeffIter rns_poly_iter(std::uint64_t rns_num)
        {
            return { &rns_coefficients(0, rns_num) };
        }

        inline seal::util::ConstCoeffIter const_rns_poly_iter(std::uint64_t rns_num) const
        {
            return { &rns_coefficients(0, rns_num) };
        }
        CiphertextFracture operator*(const CiphertextFracture &ctxf) const;
        bool operator==(const PolynomialFracture &other) const;

        /** Returns the size of the object in bytes, as if written to a stream. **/
        std::streamoff save_size(compr_mode_type compr_mode) const;

        inline std::streamoff save(
            std::ostream &stream, compr_mode_type compr_mode = Serialization::compr_mode_default) const
        {
            using namespace std::placeholders;
            return Serialization::Save(
                std::bind(&PolynomialFracture::save_members, this, _1), save_size(compr_mode_type::none), stream,
                compr_mode, false);
        }

        inline std::streamoff load(const seal_byte *in, std::size_t size)
        {
            using namespace std::placeholders;
            return Serialization::Load(std::bind(&PolynomialFracture::load_members, this, _1, _2), in, size, false);
        }

        inline std::streamoff load(std::istream &in)
        {
            using namespace std::placeholders;
            return Serialization::Load(std::bind(&PolynomialFracture::load_members, this, _1, _2), in, false);
        }

    private:
        // Saves the object to an output stream.
        void save_members(std::ostream &stream) const;
        void load_members(std::istream &stream, SEALVersion version);

        std::size_t compute_data_size() const;
    };

    /**
     * Can compute fractures on demand.
     * assumes that the polynomial is already in NTT form, and is in the correct moduli.
     * (Assumes first encryption parameters for the polynomial)
     */
    class PolynomialFracturingTool
    {
    public:
        /**
         * Extracts a single fracture out of numerous possible fractures over a polynomial.
         * @param p the plaintext to extract the fracture from.
         * @param ctx the SEAL context.
         * @param num_fractures is the number of expected fractures to be made over the polynomial.
         * @param index the index of the fracture (0 <= index < num_fractures).
         * @return a PolynomialFracture object.
         */
        static PolynomialFracture compute_fracture(
            const seal::Plaintext &p, const seal::SEALContext &ctx, uint64_t num_fractures, uint64_t index);

        /**
         * Computes a single fracture out of numerous possible fractures over a polynomial.
         * @param num_fractures is the number of expected fractures to be made over the polynomial.
         * @param index the index of the fracture (0 <= index < num_fractures).
         * @param coeff_count the number of coefficients in the polynomial.
         * @param coeff_modulus_size the number of moduli in the polynomial. (used to determine the number of RNS slots)
         * @param constiter an iterator over the polynomial coefficients.
         * @return a PolynomialFracture object.
         */
        static PolynomialFracture compute_fracture(
            uint64_t num_fractures, uint64_t index, size_t coeff_count, uint64_t coeff_modulus_size,
            const util::ConstRNSIter &constiter);
    };

    /**
     * Fractures a polynomial into multiple PolynomialFracture.
     * these polynomials are used to perform multiplication between a ciphertext and a plaintext in a distributed
     * manner.
     */
    class Polynomial
    {
    public:
        explicit Polynomial(
            const seal::Plaintext &p, const seal::SEALContext &context, std::uint64_t _num_fractures) noexcept
            : fractures(), num_fractures(_num_fractures)
        {
            seal::util::ConstRNSIter iter(p.data(), context.first_context_data()->parms().poly_modulus_degree());

            auto cntx_data = context.get_context_data(p.parms_id());
            const EncryptionParameters &params = cntx_data->parms();
            fractures.reserve(num_fractures);

            for (std::uint64_t i = 0; i < num_fractures; ++i)
            {
                fractures.push_back(compute_fracture(
                    iter, params.coeff_modulus().size(), params.poly_modulus_degree(), num_fractures, i));
            }
        }

        const PolynomialFracture &get_fracture(std::uint64_t index);
        const PolynomialFracture &operator[](std::uint64_t index);

        /**
         * states the size of 'save' for the polynomial to an output stream. doesn't saves the fracture::Essence object.
         * @param compr_mode
         * @return the save size.
         */
        std::streamoff save_size(compr_mode_type compr_mode) const;

    private:
        void save_members(std::ostream &stream) const;

        static PolynomialFracture compute_fracture(
            const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
            uint64_t num_fractures, uint64_t index);

        std::vector<PolynomialFracture> fractures;
        std::uint64_t num_fractures;
    };
} // namespace seal::fractures