#include "epidemiology/common.h"
#include <gtest/gtest.h>

TEST(TestCommon, tensor_dimension_prods)
{
    ASSERT_EQ(std::vector<size_t>({4, 12, 24, 24}), epi::tensor_dimension_prods({1, 2, 3, 4}));
    ASSERT_EQ(std::vector<size_t>({2, 4, 8, 16}), epi::tensor_dimension_prods({2, 2, 2, 2}));
}

TEST(TestCommon, flatten_index_order_two)
{
    size_t n = 3;
    size_t m = 4;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            auto flat_index = epi::flatten_index({i, j}, {n, m});
            ASSERT_LT(flat_index, n * m);
            ASSERT_EQ(m * i + j, flat_index);
        }
    }
}

TEST(TestCommon, flatten_index_order_three)
{
    size_t n = 3;
    size_t m = 4;
    size_t l = 5;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            for (size_t k = 0; k < l; ++k) {
                auto flat_index = epi::flatten_index({i, j, k}, {n, m, l});
                ASSERT_LT(flat_index, n * m * l);
                ASSERT_EQ(l * m * i + l * j + k, flat_index);
            }
        }
    }
}

TEST(TestCommon, unravel_index_order_two)
{
    size_t n = 3;
    size_t m = 4;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            auto indices = epi::unravel_index(m * i + j, {n, m});
            ASSERT_EQ(i, indices[0]);
            ASSERT_EQ(j, indices[1]);
        }
    }
}

TEST(TestCommon, unravel_index_order_three)
{
    size_t n = 3;
    size_t m = 4;
    size_t l = 5;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            for (size_t k = 0; k < l; ++k) {
                auto indices = epi::unravel_index(l * m * i + l * j + k, {n, m, l});
                ASSERT_EQ(i, indices[0]);
                ASSERT_EQ(j, indices[1]);
                ASSERT_EQ(k, indices[2]);
            }
        }
    }
}
