#include "epidemiology/utils/custom_index_array.h"
#include <gtest/gtest.h>

// A category can be a tag
struct Dim1;

// A category can be an enum class
enum class Dim2
{
    Male,
    Female,
    Count = 2
};

// We can derive a tag from epi::Index for a shorter notation
struct Dim3 : epi::Index<Dim3> {
    Dim3(size_t val) : epi::Index<Dim3>(val){}
};

TEST(CustomIndexArray, sizesAndDimensions)
{
    using ArrayType1 = epi::CustomIndexArray<float, Dim1>;
    ArrayType1 array1({epi::Index<Dim1>(3)});
    ASSERT_EQ(array1.numel(), 3);
    ASSERT_EQ(array1.size<Dim1>(), epi::Index<Dim1>(3));

    using ArrayType2 = epi::CustomIndexArray<double, Dim2>;
    ArrayType2 array2({Dim2::Count});
    ASSERT_EQ(array2.numel(), 2);
    ASSERT_EQ(array2.size<Dim2>(), Dim2(2));

    using ArrayType3 = epi::CustomIndexArray<double, Dim3>;
    ArrayType3 array3({Dim3(4)});
    ASSERT_EQ(array3.numel(), 4);
    ASSERT_EQ(array3.size<Dim3>(), Dim3(4));

    using ArrayType4 = epi::CustomIndexArray<float, Dim1, Dim2>;
    ArrayType4 array4({epi::Index<Dim1>(3), Dim2::Count});
    ASSERT_EQ(array4.numel(), 6);
    ASSERT_EQ(array4.size<Dim1>(), epi::Index<Dim1>(3));
    ASSERT_EQ(array4.size<Dim2>(), Dim2(2));

    using ArrayType5 = epi::CustomIndexArray<char, Dim2, Dim1>;
    ArrayType5 array5({Dim2::Count, epi::Index<Dim1>(3)});
    ASSERT_EQ(array5.numel(), 6);
    ASSERT_EQ(array5.size<Dim2>(), Dim2(2));
    ASSERT_EQ(array5.size<Dim1>(), epi::Index<Dim1>(3));

    using ArrayType6 = epi::CustomIndexArray<char, Dim1, Dim2, Dim3>;
    ArrayType6 array6({epi::Index<Dim1>(3), Dim2::Count, Dim3(4)});
    ASSERT_EQ(array6.numel(), 24);
    ASSERT_EQ(array6.size<Dim1>(), epi::Index<Dim1>(3));
    ASSERT_EQ(array6.size<Dim2>(), epi::Index<Dim2>(2));
    ASSERT_EQ(array6.size<Dim3>(), epi::Index<Dim3>(4));
}

TEST(CustomIndexArray, ConstantInitialization)
{
    using ArrayType = epi::CustomIndexArray<double, Dim1, Dim2>;

    ArrayType array({epi::Index<Dim1>(3), epi::Index<Dim2>(2)}, 42.);
    for (auto i = epi::Index<Dim1>(0); i < array.size<Dim1>(); ++i) {
        for (auto j = epi::Index<Dim2>(0); j < array.size<Dim2>(); j++) {
             ASSERT_DOUBLE_EQ((array[{i, j}]), 42.);
        }
    }
}

TEST(CustomIndexArray, forEach)
{
    using ArrayType = epi::CustomIndexArray<int, Dim1, Dim2>;

    int counter=0;
    ArrayType array({epi::Index<Dim1>(3), epi::Index<Dim2>(2)});
    for (auto i = epi::Index<Dim1>(0); i < array.size<Dim1>(); ++i) {
        for (auto j = epi::Index<Dim2>(0); j < array.size<Dim2>(); j++) {
            array[{i,j}] = counter++;
        }
    }
    counter = 0;
    for (auto& v : array){
        ASSERT_DOUBLE_EQ(v, counter++);
    }
}

TEST(CustomIndexArray, GetFlatIndex1D)
{
    using ArrayType1 = epi::CustomIndexArray<std::string, Dim1>;
    ArrayType1 array1({epi::Index<Dim1>(3)});
    ASSERT_EQ(array1.get_flat_index({epi::Index<Dim1>(0)}), 0);
    ASSERT_EQ(array1.get_flat_index({epi::Index<Dim1>(1)}), 1);
    ASSERT_EQ(array1.get_flat_index({epi::Index<Dim1>(2)}), 2);

    using ArrayType2 = epi::CustomIndexArray<bool, Dim2>;
    ArrayType2 array2({epi::Index<Dim2>(2)});
    ASSERT_EQ(array2.get_flat_index({epi::Index<Dim2>(0)}), 0);
    ASSERT_EQ(array2.get_flat_index({epi::Index<Dim2>(1)}), 1);
}

TEST(CustomIndexArray, GetFlatIndex2D)
{
    using ArrayType = epi::CustomIndexArray<double, Dim1, Dim2, Dim3>;
    ArrayType array({epi::Index<Dim1>(3), epi::Index<Dim2>(2), epi::Index<Dim3>(4)});
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(0), epi::Index<Dim3>(0)}),  0);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(0), epi::Index<Dim3>(1)}),  1);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(0), epi::Index<Dim3>(2)}),  2);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(0), epi::Index<Dim3>(3)}),  3);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(1), epi::Index<Dim3>(0)}),  4);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(1), epi::Index<Dim3>(1)}),  5);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(1), epi::Index<Dim3>(2)}),  6);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(0), epi::Index<Dim2>(1), epi::Index<Dim3>(3)}),  7);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(0), epi::Index<Dim3>(0)}),  8);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(0), epi::Index<Dim3>(1)}),  9);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(0), epi::Index<Dim3>(2)}), 10);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(0), epi::Index<Dim3>(3)}), 11);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(1), epi::Index<Dim3>(0)}), 12);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(1), epi::Index<Dim3>(1)}), 13);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(1), epi::Index<Dim3>(2)}), 14);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(1), epi::Index<Dim2>(1), epi::Index<Dim3>(3)}), 15);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(0), epi::Index<Dim3>(0)}), 16);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(0), epi::Index<Dim3>(1)}), 17);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(0), epi::Index<Dim3>(2)}), 18);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(0), epi::Index<Dim3>(3)}), 19);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(1), epi::Index<Dim3>(0)}), 20);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(1), epi::Index<Dim3>(1)}), 21);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(1), epi::Index<Dim3>(2)}), 22);
    ASSERT_EQ(array.get_flat_index({epi::Index<Dim1>(2), epi::Index<Dim2>(1), epi::Index<Dim3>(3)}), 23);
}

TEST(CustomIndexArray, get_and_set)
{
    using ArrayType = epi::CustomIndexArray<int, Dim1, Dim2, Dim3>;
    ArrayType array({epi::Index<Dim1>(3), epi::Index<Dim2>(2), epi::Index<Dim3>(4)}, 42);

    ASSERT_EQ((array[{(epi::Index<Dim1>)0, (epi::Index<Dim2>)1, (epi::Index<Dim3>)2}]), 42);
    ASSERT_EQ(array.array()[6], 42);

    array[{epi::Index<Dim1>(0), epi::Index<Dim2>(1), epi::Index<Dim3>(2)}] = 18;
    ASSERT_EQ((array[{(epi::Index<Dim1>)0, (epi::Index<Dim2>)1, (epi::Index<Dim3>)2}]), 18);
    ASSERT_EQ(array.array()[6], 18);

}

struct Tag0 : public epi::Index<Tag0> {
    Tag0(size_t val) : epi::Index<Tag0>(val){}
};
struct Tag1 : public epi::Index<Tag1> {
    Tag1(size_t val) : epi::Index<Tag1>(val){}
};
struct Tag2 : public epi::Index<Tag2> {
    Tag2(size_t val) : epi::Index<Tag2>(val){}
};
struct Tag3 : public epi::Index<Tag3> {
    Tag3(size_t val) : epi::Index<Tag3>(val){}
};

TEST(CustomIndexArray, slice)
{
    epi::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3>
            array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(1)}]), 3.1459);

    // change value of all entries with indices 1:2 along the third dimension Tag2.
    // This should change the flat incides 2,3,4,5,8,9,10,11
    int idx = 0;
    for (auto& v : array.slice<Tag2>({1,2}))
    {
        v = idx++;
    }
    ASSERT_EQ(array.slice<Tag2>({1,2}).numel(), 8);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(0)}]), 0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(1)}]), 1);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(0)}]), 2);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(1)}]), 3);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(0)}]), 4);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(1)}]), 5);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(0)}]), 6);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(1)}]), 7);
}

TEST(CustomIndexArray, slice_with_stride)
{
    epi::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3>
            array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change value of all entries with indices 0 or 2 along the third dimension Tag2.
    // This should change the flat incides 0,1,4,5,6,7,10,11
    int idx = 0;
    for (auto& v : array.slice<Tag2>({0,2,2}))
    {
        v = idx++;
    }
    ASSERT_EQ(array.slice<Tag2>({0,2,2}).numel(), 8);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(0)}]), 0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(1)}]), 1);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(0)}]), 2);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(1)}]), 3);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(0)}]), 4);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(1)}]), 5);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(0)}]), 6);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(1)}]), 7);
}

TEST(CustomIndexArray, slice_as_array)
{
    epi::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3>
            array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change value of all entries with indices 1:2 along the third dimension Tag2.
    // This should change the flat incides 2,3,4,5,8,9,10,11
    int idx = 0;
    for (auto& v : array)
    {
         v = idx++;
    }

    auto slice_array = array.slice<Tag2>({1,2}).as_array();
    ASSERT_EQ(slice_array.numel(), 8);
    ASSERT_EQ(slice_array.size<Tag0>(), epi::Index<Tag0>(1));
    ASSERT_EQ(slice_array.size<Tag1>(), epi::Index<Tag1>(2));
    ASSERT_EQ(slice_array.size<Tag2>(), epi::Index<Tag2>(2));
    ASSERT_EQ(slice_array.size<Tag3>(), epi::Index<Tag3>(2));

    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(0), Tag2(0), Tag3(0)}]), 2);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(0), Tag2(0), Tag3(1)}]), 3);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(0), Tag2(1), Tag3(0)}]), 4);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(0), Tag2(1), Tag3(1)}]), 5);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(1), Tag2(0), Tag3(0)}]), 8);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(1), Tag2(0), Tag3(1)}]), 9);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(1), Tag2(1), Tag3(0)}]), 10);
    ASSERT_DOUBLE_EQ((slice_array[{Tag0(0), Tag1(1), Tag2(1), Tag3(1)}]), 11);
}

TEST(CustomIndexArray, chained_slice)
{
    epi::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3>
            array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change all values that have an index 1 along Tag1 and an index 1 or 2 along Tag3
    // This corresponds to the last four flat indices 8,9,10,11
    for (auto& v : array.slice<Tag1>({1,1}).slice<Tag2>({1,2}))
    {
        v = 42;
    }

    ASSERT_EQ(array.slice<Tag1>({1,1}).slice<Tag2>({1,2}).numel(), 4);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(0)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0), Tag3(1)}]), 3.1459);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(0)}]), 42);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1), Tag3(1)}]), 42);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(0)}]), 42);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2), Tag3(1)}]), 42);
}
