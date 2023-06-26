#include "fractured_polynomial.h"
#include "fractured_ciphertext.h"

namespace seal::fractures
{

    PolynomialFracture Polynomial::compute_fracture(
        const seal::util::ConstRNSIter &rns_iter, std::uint64_t modulus_size, uint64_t num_coeffs,
        uint64_t num_fractures, uint64_t index)
    {
        PolynomialFracture fracture(index, num_coeffs / num_fractures, modulus_size);

        std::uint64_t rns_num = 0;
        // We run up to modulus_size-1 because the inner loop moves from first position to the end.
        // for example, when mudulus_size==1, we move from 0 and up to 1,  thus we go through the whole range.
        std::for_each_n(seal::util::iter(rns_iter), modulus_size - 1, [&](auto coef_iter) {
            coef_iter += index * num_coeffs / num_fractures;
            auto write_into_iter = fracture.rns_poly_iter(rns_num);
            for (uint64_t j = 0; j < num_coeffs / num_fractures; ++j)
            {
                *write_into_iter = *coef_iter;
                //                fracture.rns_coefficients(j, rns_num) = *coef_iter;
                coef_iter++;
                write_into_iter++;
            }

            rns_num++;
        });

        return fracture;
    }

    Polynomial::Polynomial(const seal::util::ConstRNSIter &p, Essence e, std::uint64_t _num_fractures) noexcept
        : num_fractures(_num_fractures), essence(std::move(e))
    {
        fractures.reserve(num_fractures);

        for (std::uint64_t i = 0; i < num_fractures; ++i)
        {
            fractures.push_back(
                compute_fracture(p, essence.parms.coeff_modulus().size(), essence.coeff_count, num_fractures, i));
        }
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

        for (std::uint64_t j = 0; j < rns_coefficients.data.size(); ++j)
        {
            if (rns_coefficients.data[j] != other.rns_coefficients.data[j])
            {
                return false;
            }
        }

        return true;
    }

} // namespace seal::fractures