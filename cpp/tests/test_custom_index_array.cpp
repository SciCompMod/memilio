#include "epidemiology/utils/custom_index_array.h"
#include <gtest/gtest.h>

enum class Dim1
{
    Idx1,
    Idx2,
    Idx3,
    Count = 3
};

enum class Dim2
{
    Idx1,
    Idx2,
    Count = 2
};

enum class Dim3
{
    Idx1,
    Idx2,
    Idx3,
    Idx4,
    Count = 4
};

TEST(CustomIndexArray, sizesAndDimensions)
{
    using ArrayType1 = epi::CustomIndexArray<float, Dim1>;
    ASSERT_EQ(ArrayType1::size(), 3);
    ASSERT_EQ(ArrayType1::dimensions.size(), 1);
    ASSERT_EQ(ArrayType1::dimensions[0], 3);

    using ArrayType2 = epi::CustomIndexArray<double, Dim2>;
    ASSERT_EQ(ArrayType2::size(), 2);
    ASSERT_EQ(ArrayType2::dimensions.size(), 1);
    ASSERT_EQ(ArrayType2::dimensions[0], 2);

    using ArrayType3 = epi::CustomIndexArray<double, Dim3>;
    ASSERT_EQ(ArrayType3::size(), 4);
    ASSERT_EQ(ArrayType3::dimensions.size(), 1);
    ASSERT_EQ(ArrayType3::dimensions[0], 4);

    using ArrayType4 = epi::CustomIndexArray<float, Dim1, Dim2>;
    ASSERT_EQ(ArrayType4::size(), 6);
    ASSERT_EQ(ArrayType4::dimensions.size(), 2);
    ASSERT_EQ(ArrayType4::dimensions[0], 3);
    ASSERT_EQ(ArrayType4::dimensions[1], 2);

    using ArrayType5 = epi::CustomIndexArray<char, Dim2, Dim1>;
    ASSERT_EQ(ArrayType5::size(), 6);
    ASSERT_EQ(ArrayType5::dimensions.size(), 2);
    ASSERT_EQ(ArrayType5::dimensions[0], 2);
    ASSERT_EQ(ArrayType5::dimensions[1], 3);

    using ArrayType6 = epi::CustomIndexArray<char, Dim1, Dim2, Dim3>;
    ASSERT_EQ(ArrayType6::size(), 24);
    ASSERT_EQ(ArrayType6::dimensions.size(), 3);
    ASSERT_EQ(ArrayType6::dimensions[0], 3);
    ASSERT_EQ(ArrayType6::dimensions[1], 2);
    ASSERT_EQ(ArrayType6::dimensions[2], 4);
}

TEST(CustomIndexArray, ConstantInitialization)
{
    using ArrayType = epi::CustomIndexArray<double, Dim1, Dim2>;

    ArrayType array(42.);
    for (int i=0; i<(int)Dim1::Count; i++) {
        for (int j=0; j<(int)Dim2::Count; j++) {
             ASSERT_DOUBLE_EQ((array[{(Dim1)i, (Dim2)j}]), 42.);
        }
    }
}

TEST(CustomIndexArray, GetFlatIndex1D)
{
    using ArrayType1 = epi::CustomIndexArray<std::string, Dim1>;
    ASSERT_EQ(ArrayType1::get_flat_index(Dim1::Idx1), 0);
    ASSERT_EQ(ArrayType1::get_flat_index(Dim1::Idx2), 1);
    ASSERT_EQ(ArrayType1::get_flat_index(Dim1::Idx3), 2);

    using ArrayType2 = epi::CustomIndexArray<bool, Dim2>;
    ASSERT_EQ(ArrayType2::get_flat_index(Dim2::Idx1), 0);
    ASSERT_EQ(ArrayType2::get_flat_index(Dim2::Idx2), 1);
}

TEST(CustomIndexArray, GetFlatIndex2D)
{
    using ArrayType = epi::CustomIndexArray<double, Dim1, Dim2, Dim3>;
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx1, Dim3::Idx1),  0);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx1, Dim3::Idx2),  1);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx1, Dim3::Idx3),  2);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx1, Dim3::Idx4),  3);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx2, Dim3::Idx1),  4);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx2, Dim3::Idx2),  5);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx2, Dim3::Idx3),  6);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx1, Dim2::Idx2, Dim3::Idx4),  7);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx1, Dim3::Idx1),  8);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx1, Dim3::Idx2),  9);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx1, Dim3::Idx3), 10);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx1, Dim3::Idx4), 11);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx2, Dim3::Idx1), 12);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx2, Dim3::Idx2), 13);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx2, Dim3::Idx3), 14);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx2, Dim2::Idx2, Dim3::Idx4), 15);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx1, Dim3::Idx1), 16);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx1, Dim3::Idx2), 17);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx1, Dim3::Idx3), 18);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx1, Dim3::Idx4), 19);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx2, Dim3::Idx1), 20);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx2, Dim3::Idx2), 21);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx2, Dim3::Idx3), 22);
    ASSERT_EQ(ArrayType::get_flat_index(Dim1::Idx3, Dim2::Idx2, Dim3::Idx4), 23);
}

TEST(CustomIndexArray, get)
{
    using ArrayType = epi::CustomIndexArray<int, Dim1, Dim2, Dim3>;
    ArrayType array(42);

    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 42);
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 42);
    ASSERT_EQ(array.array()[6], 42);

    array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}] = 18;
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 18);
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 18);
    ASSERT_EQ(array.array()[6], 18);
}

TEST(CustomIndexArray, set)
{
    using ArrayType = epi::CustomIndexArray<int, Dim1, Dim2, Dim3>;
    ArrayType array(42);

    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 42);
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 42);
    ASSERT_EQ(array.array()[6], 42);

    array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}] = 18;
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 18);
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 18);
    ASSERT_EQ(array.array()[6], 18);

    array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}] = 99;
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 99);
    ASSERT_EQ((array[{Dim1::Idx1, Dim2::Idx2, Dim3::Idx3}]), 99);
    ASSERT_EQ(array.array()[6], 99);


}
