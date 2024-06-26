#include "fractured_polynomial.h"
#include "fractured_ciphertext.h"

namespace seal::fractures
{

    PolynomialFracture Polynomial::compute_fracture(
        const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
        uint64_t num_fractures, uint64_t index)
    {
        return PolynomialFracturingTool::compute_fracture(num_fractures, index, num_coeffs, modulus_size, rns_iter);
    }

    const PolynomialFracture &Polynomial::get_fracture(std::uint64_t index)
    {
        if (index >= fractures.size())
        {
            throw std::invalid_argument("index out of bounds");
        }

        return fractures[index];
    }

    const PolynomialFracture &Polynomial::operator[](std::uint64_t index)
    {
        return fractures[index];
    }

    std::streamoff Polynomial::save_size(compr_mode_type compr_mode) const
    {
        auto members_size = sizeof(std::uint64_t); // num_fractures as a start.

        // then each fracture and its size.
        for (auto &pf : fractures)
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

    void Polynomial::save_members(std::ostream &stream) const
    {
        auto old_except_mask = stream.exceptions();
        try
        {
            stream.write(reinterpret_cast<const char *>(&num_fractures), sizeof(std::uint64_t));
            // each element should have some kind of size.
            for (auto &pf : fractures)
            {
                auto save_size = pf.save_size(compr_mode_type::none);
                stream.write(reinterpret_cast<const char *>(&save_size), sizeof(std::uint64_t));
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

    CiphertextFracture PolynomialFracture::operator*(const CiphertextFracture &ctxf) const
    {
        PolynomialFracture tmp(*this);
        return ctxf * tmp;
    }

    bool PolynomialFracture::operator==(const PolynomialFracture &other) const
    {
        if (rns_coefficients.data.size() != other.rns_coefficients.data.size())
        {
            return false;
        }
        if (coeff_count != other.coeff_count)
        {
            return false;
        }
        if (fracture_index != other.fracture_index)
        {
            return false;
        }

        for (std::uint64_t j = 0; j < rns_coefficients.data.size(); ++j)
        {
            if (rns_coefficients.data[j] != other.rns_coefficients.data[j])
            {
                return false;
            }
        }

        return true;
    }

    std::streamoff PolynomialFracture::save_size(compr_mode_type compr_mode) const
    {
        size_t members_size = Serialization::ComprSizeEstimate(
            util::add_safe(
                sizeof(std::uint64_t), // data_per_rns (rows)
                sizeof(std::uint64_t), // rns_coefficients (cols)
                sizeof(std::uint64_t), // fracture_index
                compute_data_size()),
            compr_mode);

        return util::safe_cast<std::streamoff>(util::add_safe(sizeof(Serialization::SEALHeader), members_size));
    }
    size_t PolynomialFracture::compute_data_size() const
    {
        std::size_t data_size = rns_coefficients.cols * rns_coefficients.rows * sizeof(uint64_t);
        return data_size;
    }

    void PolynomialFracture::save_members(std::ostream &stream) const
    {
        auto old_except_mask = stream.exceptions();
        try
        {
            // Throw exceptions on std::ios_base::badbit and std::ios_base::failbit
            stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);

            auto num_values_byte = this->rns_coefficients.rows;
            stream.write(reinterpret_cast<const char *>(&num_values_byte), sizeof(uint64_t));

            auto modulus_size_byte = this->rns_coefficients.cols;
            stream.write(reinterpret_cast<const char *>(&modulus_size_byte), sizeof(uint64_t));

            auto index_num_byte = this->fracture_index;
            stream.write(reinterpret_cast<const char *>(&index_num_byte), sizeof(uint64_t));

            auto data_size = compute_data_size();
            if (data_size > std::numeric_limits<std::streamsize>::max())
            {
                throw std::runtime_error("data_size is too large");
            }

            stream.write(reinterpret_cast<const char *>(&this->rns_coefficients.data[0]), static_cast<long>(data_size));
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

    void PolynomialFracture::load_members(std::istream &stream, SEAL_MAYBE_UNUSED SEALVersion version)
    {
        auto old_except_mask = stream.exceptions();
        try
        {
            std::uint64_t rows;
            std::uint64_t cols;
            std::uint64_t index;

            stream.read(reinterpret_cast<char *>(&rows), sizeof(uint64_t));
            stream.read(reinterpret_cast<char *>(&cols), sizeof(uint64_t));
            stream.read(reinterpret_cast<char *>(&index), sizeof(uint64_t));

            this->coeff_count = rows;
            this->fracture_index = index;

            this->rns_coefficients = seal::util::matrix<std::uint64_t>(rows, cols);

            stream.read(
                reinterpret_cast<char *>(&this->rns_coefficients.data[0]),
                static_cast<long>(rows * cols * sizeof(std::uint64_t)));
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

    PolynomialFracture PolynomialFracturingTool::compute_fracture(
        const Plaintext &p, const SEALContext &ctx, uint64_t num_fractures, uint64_t index)
    {
        auto coeff_count = p.coeff_count();

        auto contex_data = ctx.get_context_data(p.parms_id());
        auto enc_params = contex_data->parms();
        auto coeff_modulus = enc_params.coeff_modulus();

        seal::util::ConstRNSIter constiter(p.data(), coeff_count);

        return compute_fracture(num_fractures, index, coeff_count, coeff_modulus.size(), constiter);
    }

    PolynomialFracture PolynomialFracturingTool::compute_fracture(
        uint64_t num_fractures, uint64_t index, size_t coeff_count, uint64_t coeff_modulus_size,
        const util::ConstRNSIter &constiter)
    {
        PolynomialFracture fracture(index, coeff_count / num_fractures, coeff_modulus_size);

        uint64_t rns_index = 0;
        SEAL_ITERATE(seal::util::iter(constiter), coeff_modulus_size, [&](auto coef_iter) {
            coef_iter += index * coeff_count / num_fractures;
            auto write_into_iter = fracture.rns_poly_iter(rns_index++);

            for (uint64_t j = 0; j < coeff_count / num_fractures; ++j)
            {
                *write_into_iter = *coef_iter;
                //                fracture.rns_coefficients(j, rns_num) = *coef_iter;
                coef_iter++;
                write_into_iter++;
            }
        });
        return fracture;
    }
} // namespace seal::fractures