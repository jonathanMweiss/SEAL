
#include "fractured_poly.h"
seal::fractures::Poly seal::fractures::Poly::from_plaintext(
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

seal::fractures::Poly::Poly(
    seal::util::ConstRNSIter rns_iter, const seal::fractures::Essence &e, std::uint64_t num_fractures) noexcept
    : poly_data(e.coeff_count, e.coeff_modulus_size), num_fractures(num_fractures), essence(e)
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
        // coef_iter += index * e.coeff_count / num_fractures;
        auto frac_coef_iter = rns_poly_iter(rns_num);
        Poly::memcpy(frac_coef_iter, coef_iter, essence.coeff_count);
        rns_num++;
    });
}

void seal::fractures::Poly::memcpy(
    seal::util::CoeffIter &cpy_to, seal::util::ConstCoeffIter &cpy_from, std::uint64_t coeffcount)
{
    for (std::uint64_t j = 0; j < coeffcount; ++j)
    {
        *cpy_to = *cpy_from;
        cpy_from++;
        cpy_to++;
    }
}

seal::util::ConstCoeffIter seal::fractures::Poly::const_rns_poly_iter(std::uint64_t rns_num)
{
    return { &poly_data(0, rns_num) };
}

seal::util::CoeffIter seal::fractures::Poly::rns_poly_iter(std::uint64_t rns_num)
{
    return { &poly_data(0, rns_num) };
}
