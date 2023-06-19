

#pragma once
#include <cassert>
#include <cstdint>
#include <vector>

namespace seal::util
{
    template <typename T>
    class matrix
    {
    public:
        std::uint64_t rows;
        std::uint64_t cols;
        std::vector<T> data;

        matrix() : rows(0), cols(0), data()
        {}

        matrix(std::uint64_t _rows, std::uint64_t _cols) : rows(_rows), cols(_cols), data(rows * cols){};

        matrix(std::uint64_t _rows, std::uint64_t _cols, std::vector<T> _data) : rows(_rows), cols(_cols), data(_data)
        {}

        matrix(std::uint64_t _rows, std::uint64_t _cols, T cpy) : rows(_rows), cols(_cols), data(rows * cols, cpy)
        {}

        inline void resize(std::uint64_t _rows, std::uint64_t _cols)
        {
            rows = _rows;
            cols = _cols;
            data.resize(rows * cols);
        };


        [[nodiscard]] inline std::uint64_t pos(std::uint64_t row, std::uint64_t col) const
        {
            return row + col * rows;
        }

        inline const T &operator()(std::uint64_t row, std::uint64_t col) const
        {
#ifdef SEAL_DEBUG
            assert_pos(row, col);
#endif
            return data[pos(row, col)];
        }

        inline T &operator()(std::uint64_t row, std::uint64_t col)
        {
#ifdef SEAL_DEBUG
            assert_pos(row, col);
#endif
            return data[pos(row, col)];
        }

        inline void assert_pos(std::uint64_t row, std::uint64_t col) const
        {
            if (row + col * rows >= data.size())
            {
                //                std::cout << std::endl;
            }
            assert(row + col * rows < data.size());
        };

        std::uint64_t entries()
        {
            return rows * cols;
        }
    };
} // namespace
