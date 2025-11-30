/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Khoa Nguyen
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "memilio/utils/custom_index_array.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// A category can be a tag
struct Dim1;

// A category can be an enum class
enum class Dim2
{
    Male,
    Female,
    Count = 2
};

// We can derive a tag from mio::Index for a shorter notation
struct Dim3 : mio::Index<Dim3> {
    Dim3(size_t val)
        : mio::Index<Dim3>(val)
    {
    }
};

TEST(CustomIndexArray, sizesAndDimensions)
{
    using ArrayType1 = mio::CustomIndexArray<float, Dim1>;
    ArrayType1 array1({mio::Index<Dim1>(3)});
    ASSERT_EQ(array1.numel(), 3);
    ASSERT_EQ(array1.size<Dim1>(), mio::Index<Dim1>(3));

    using ArrayType2 = mio::CustomIndexArray<double, Dim2>;
    ArrayType2 array2({Dim2::Count});
    ASSERT_EQ(array2.numel(), 2);
    ASSERT_EQ(array2.size<Dim2>(), Dim2(2));

    using ArrayType3 = mio::CustomIndexArray<double, Dim3>;
    ArrayType3 array3({Dim3(4)});
    ASSERT_EQ(array3.numel(), 4);
    ASSERT_EQ(array3.size<Dim3>(), Dim3(4));

    using ArrayType4 = mio::CustomIndexArray<float, Dim1, Dim2>;
    ArrayType4 array4({mio::Index<Dim1>(3), Dim2::Count});
    ASSERT_EQ(array4.numel(), 6);
    ASSERT_EQ(array4.size<Dim1>(), mio::Index<Dim1>(3));
    ASSERT_EQ(array4.size<Dim2>(), Dim2(2));

    using ArrayType5 = mio::CustomIndexArray<char, Dim2, Dim1>;
    ArrayType5 array5({Dim2::Count, mio::Index<Dim1>(3)});
    ASSERT_EQ(array5.numel(), 6);
    ASSERT_EQ(array5.size<Dim2>(), Dim2(2));
    ASSERT_EQ(array5.size<Dim1>(), mio::Index<Dim1>(3));

    using ArrayType6 = mio::CustomIndexArray<char, Dim1, Dim2, Dim3>;
    ArrayType6 array6({mio::Index<Dim1>(3), Dim2::Count, Dim3(4)});
    ASSERT_EQ(array6.numel(), 24);
    ASSERT_EQ(array6.size<Dim1>(), mio::Index<Dim1>(3));
    ASSERT_EQ(array6.size<Dim2>(), mio::Index<Dim2>(2));
    ASSERT_EQ(array6.size<Dim3>(), mio::Index<Dim3>(4));
}

TEST(CustomIndexArray, ConstantInitialization)
{
    using ArrayType = mio::CustomIndexArray<double, Dim1, Dim2>;

    ArrayType array({mio::Index<Dim1>(3), mio::Index<Dim2>(2)}, 42.);
    for (auto i = mio::Index<Dim1>(0); i < array.size<Dim1>(); ++i) {
        for (auto j = mio::Index<Dim2>(0); j < array.size<Dim2>(); j++) {
            ASSERT_DOUBLE_EQ((array[{i, j}]), 42.);
        }
    }
}

TEST(CustomIndexArray, RangeInitialization)
{
    using ArrayType = mio::CustomIndexArray<double, Dim1, Dim2>;

    std::vector<double> values = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    ArrayType array({mio::Index<Dim1>(2), mio::Index<Dim2>(3)}, values.begin(), values.end());
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(0), mio::Index<Dim2>(0)}]), 0.1);
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(0), mio::Index<Dim2>(1)}]), 0.2);
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(0), mio::Index<Dim2>(2)}]), 0.3);
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(1), mio::Index<Dim2>(0)}]), 0.4);
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(1), mio::Index<Dim2>(1)}]), 0.5);
    ASSERT_DOUBLE_EQ((array[{mio::Index<Dim1>(1), mio::Index<Dim2>(2)}]), 0.6);
}

TEST(CustomIndexArray, equality)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2>;

    auto array = ArrayType({mio::Index<Dim1>(3), mio::Index<Dim2>(4)}, 5);
    ASSERT_EQ(array, ArrayType({mio::Index<Dim1>(3), mio::Index<Dim2>(4)}, 5));
    ASSERT_NE(array, ArrayType({mio::Index<Dim1>(3), mio::Index<Dim2>(4)}, 4));
    ASSERT_NE(array, ArrayType({mio::Index<Dim1>(4), mio::Index<Dim2>(3)}, 5));
}

TEST(CustomIndexArray, constantAssignment)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2>;

    auto array = ArrayType({mio::Index<Dim1>(3), mio::Index<Dim2>(4)}, 3);
    array      = 4;
    ASSERT_EQ(array, ArrayType({mio::Index<Dim1>(3), mio::Index<Dim2>(4)}, 4));
}

TEST(CustomIndexArray, forEach)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2>;

    int counter = 0;
    ArrayType array({mio::Index<Dim1>(3), mio::Index<Dim2>(2)});
    for (auto i = mio::Index<Dim1>(0); i < array.size<Dim1>(); ++i) {
        for (auto j = mio::Index<Dim2>(0); j < array.size<Dim2>(); j++) {
            array[{i, j}] = counter++;
        }
    }
    counter = 0;
    for (auto& v : array) {
        ASSERT_DOUBLE_EQ(v, counter++);
    }
}

TEST(CustomIndexArray, GetFlatIndex1D)
{
    using ArrayType1 = mio::CustomIndexArray<std::string, Dim1>;
    ArrayType1 array1({mio::Index<Dim1>(3)});
    ASSERT_EQ(array1.get_flat_index({mio::Index<Dim1>(0)}), 0);
    ASSERT_EQ(array1.get_flat_index({mio::Index<Dim1>(1)}), 1);
    ASSERT_EQ(array1.get_flat_index({mio::Index<Dim1>(2)}), 2);

    using ArrayType2 = mio::CustomIndexArray<bool, Dim2>;
    ArrayType2 array2({mio::Index<Dim2>(2)});
    ASSERT_EQ(array2.get_flat_index({mio::Index<Dim2>(0)}), 0);
    ASSERT_EQ(array2.get_flat_index({mio::Index<Dim2>(1)}), 1);
}

TEST(CustomIndexArray, GetFlatIndex2D)
{
    using ArrayType = mio::CustomIndexArray<double, Dim1, Dim2, Dim3>;
    ArrayType array({mio::Index<Dim1>(3), mio::Index<Dim2>(2), mio::Index<Dim3>(4)});
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(0), mio::Index<Dim3>(0)}), 0);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(0), mio::Index<Dim3>(1)}), 1);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(0), mio::Index<Dim3>(2)}), 2);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(0), mio::Index<Dim3>(3)}), 3);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(1), mio::Index<Dim3>(0)}), 4);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(1), mio::Index<Dim3>(1)}), 5);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(1), mio::Index<Dim3>(2)}), 6);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(0), mio::Index<Dim2>(1), mio::Index<Dim3>(3)}), 7);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(0), mio::Index<Dim3>(0)}), 8);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(0), mio::Index<Dim3>(1)}), 9);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(0), mio::Index<Dim3>(2)}), 10);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(0), mio::Index<Dim3>(3)}), 11);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(1), mio::Index<Dim3>(0)}), 12);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(1), mio::Index<Dim3>(1)}), 13);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(1), mio::Index<Dim3>(2)}), 14);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(1), mio::Index<Dim2>(1), mio::Index<Dim3>(3)}), 15);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(0), mio::Index<Dim3>(0)}), 16);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(0), mio::Index<Dim3>(1)}), 17);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(0), mio::Index<Dim3>(2)}), 18);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(0), mio::Index<Dim3>(3)}), 19);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(1), mio::Index<Dim3>(0)}), 20);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(1), mio::Index<Dim3>(1)}), 21);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(1), mio::Index<Dim3>(2)}), 22);
    ASSERT_EQ(array.get_flat_index({mio::Index<Dim1>(2), mio::Index<Dim2>(1), mio::Index<Dim3>(3)}), 23);
}

TEST(CustomIndexArray, get_and_set)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2, Dim3>;
    ArrayType array({mio::Index<Dim1>(3), mio::Index<Dim2>(2), mio::Index<Dim3>(4)}, 42);

    ASSERT_EQ((array[{(mio::Index<Dim1>)0, (mio::Index<Dim2>)1, (mio::Index<Dim3>)2}]), 42);
    ASSERT_EQ(array.array()[6], 42);

    array[{mio::Index<Dim1>(0), mio::Index<Dim2>(1), mio::Index<Dim3>(2)}] = 18;
    ASSERT_EQ((array[{(mio::Index<Dim1>)0, (mio::Index<Dim2>)1, (mio::Index<Dim3>)2}]), 18);
    ASSERT_EQ(array.array()[6], 18);
}

struct Tag0 : public mio::Index<Tag0> {
    Tag0(size_t val)
        : mio::Index<Tag0>(val)
    {
    }
};
struct Tag1 : public mio::Index<Tag1> {
    Tag1(size_t val)
        : mio::Index<Tag1>(val)
    {
    }
};
struct Tag2 : public mio::Index<Tag2> {
    Tag2(size_t val)
        : mio::Index<Tag2>(val)
    {
    }
};
struct Tag3 : public mio::Index<Tag3> {
    Tag3(size_t val)
        : mio::Index<Tag3>(val)
    {
    }
};

TEST(CustomIndexArray, slice)
{
    mio::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3> array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

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
    for (auto& v : array.slice<Tag2>({1, 2})) {
        v = idx++;
    }
    ASSERT_EQ(array.slice<Tag2>({1, 2}).numel(), 8);
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

TEST(CustomIndexArray, slice_single)
{
    mio::CustomIndexArray<int, Tag0, Tag1> array{{Tag0(2), Tag1(3)}};
    array[{mio::Index<Tag0>(0), mio::Index<Tag1>(0)}] = 1;
    array[{mio::Index<Tag0>(0), mio::Index<Tag1>(1)}] = 2;
    array[{mio::Index<Tag0>(0), mio::Index<Tag1>(2)}] = 3;
    array[{mio::Index<Tag0>(1), mio::Index<Tag1>(0)}] = 4;
    array[{mio::Index<Tag0>(1), mio::Index<Tag1>(1)}] = 5;
    array[{mio::Index<Tag0>(1), mio::Index<Tag1>(2)}] = 6;

    ASSERT_THAT(array.slice(mio::Index<Tag1>(1)), testing::ElementsAre(2, 5));
}

TEST(CustomIndexArray, slice_with_stride)
{
    mio::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3> array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change value of all entries with indices 0 or 2 along the third dimension Tag2.
    // This should change the flat incides 0,1,4,5,6,7,10,11
    int idx = 0;
    for (auto& v : array.slice<Tag2>({0, 2, 2})) {
        v = idx++;
    }
    ASSERT_EQ(array.slice<Tag2>({0, 2, 2}).numel(), 8);
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
    mio::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3> array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change value of all entries with indices 1:2 along the third dimension Tag2.
    // This should change the flat incides 2,3,4,5,8,9,10,11
    int idx = 0;
    for (auto& v : array) {
        v = idx++;
    }

    auto slice_array = array.slice<Tag2>({1, 2}).as_array();
    ASSERT_EQ(slice_array.numel(), 8);
    ASSERT_EQ(slice_array.size<Tag0>(), mio::Index<Tag0>(1));
    ASSERT_EQ(slice_array.size<Tag1>(), mio::Index<Tag1>(2));
    ASSERT_EQ(slice_array.size<Tag2>(), mio::Index<Tag2>(2));
    ASSERT_EQ(slice_array.size<Tag3>(), mio::Index<Tag3>(2));

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
    mio::CustomIndexArray<double, Tag0, Tag1, Tag2, Tag3> array({Tag0(1), Tag1(2), Tag2(3), Tag3(2)}, 3.1459);

    // change all values that have an index 1 along Tag1 and an index 1 or 2 along Tag3
    // This corresponds to the last four flat indices 8,9,10,11
    for (auto& v : array.slice<Tag1>({1, 1}).slice<Tag2>({1, 2})) {
        v = 42;
    }

    ASSERT_EQ(array.slice<Tag1>({1, 1}).slice<Tag2>({1, 2}).numel(), 4);
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

TEST(CustomIndexArray, sliceConstantAssignment)
{
    mio::CustomIndexArray<double, Tag0, Tag1, Tag2> array({Tag0(1), Tag1(2), Tag2(3)}, 3.1459);

    array.slice<Tag2>({1, 2}) = 42.0; //sequence of 2 indices along one dimension
    array.slice(Tag2(0))      = 17.0; //single index along one dimension

    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(0)}]), 17.0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(1)}]), 42.0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(0), Tag2(2)}]), 42.0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(0)}]), 17.0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(1)}]), 42.0);
    ASSERT_DOUBLE_EQ((array[{Tag0(0), Tag1(1), Tag2(2)}]), 42.0);
}

TEST(CustomIndexArray, resize_all)
{
    mio::CustomIndexArray<double, Tag0, Tag1> array{{Tag0(1), Tag1(2)}, 1.23};
    array.resize({Tag0(2), Tag1(4)});
    ASSERT_EQ(array.size().indices, std::make_tuple(Tag0{2}, Tag1{4}));
    ASSERT_EQ(array.numel(), 8);
}

TEST(CustomIndexArray, resize_one_dimension)
{
    mio::CustomIndexArray<double, Tag0, Tag1> array{{Tag0(1), Tag1(2)}, 1.23};
    array.resize(Tag0(3));
    ASSERT_EQ(array.size().indices, std::make_tuple(Tag0{3}, Tag1{2}));
    ASSERT_EQ(array.numel(), 6);
}

// Additional test for checking the set_multiple functionality
TEST(CustomIndexArray, setMultiple_validIndices)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2>;
    // Initialize with zeros
    ArrayType array({mio::Index<Dim1>(3), Dim2::Count}, 0);

    // Define indices to set
    std::vector<ArrayType::Index> indices = {{mio::Index<Dim1>(0), Dim2::Male}, {mio::Index<Dim1>(2), Dim2::Female}};

    // Set these indices to 42
    array.set_multiple(indices, 42);

    // Verify that the correct indices have been updated
    EXPECT_EQ((array[{mio::Index<Dim1>(0), Dim2::Male}]), 42);
    EXPECT_EQ((array[{mio::Index<Dim1>(2), Dim2::Female}]), 42);

    // Verify that other indices are unchanged
    EXPECT_EQ((array[{mio::Index<Dim1>(0), Dim2::Female}]), 0);
    EXPECT_EQ((array[{mio::Index<Dim1>(1), Dim2::Male}]), 0);
    EXPECT_EQ((array[{mio::Index<Dim1>(1), Dim2::Female}]), 0);
    EXPECT_EQ((array[{mio::Index<Dim1>(2), Dim2::Male}]), 0);
}

TEST(CustomIndexArray, setMultiple_emptyIndices)
{
    using ArrayType = mio::CustomIndexArray<int, Dim1, Dim2>;
    // Initialize with fives
    ArrayType array({mio::Index<Dim1>(2), Dim2::Count}, 5);

    // Empty vector of indices
    std::vector<ArrayType::Index> indices;

    // Attempt to set multiple indices to 42
    array.set_multiple(indices, 42);

    // Verify that all entries remain unchanged
    for (int age = 0; age < 2; ++age) {
        for (int gender = 0; gender < static_cast<int>(Dim2::Count); ++gender) {
            EXPECT_EQ((array[{mio::Index<Dim1>(age), static_cast<Dim2>(gender)}]), 5);
        }
    }
}

TEST(CustomIndexArray, convert_floating_point)
{
    // case: float arrays with mixed precision; expect same result after conversion through double precision
    float value_float = 0.10005f;
    auto size_dim1    = mio::Index<Dim1>(3);
    mio::CustomIndexArray<float, Dim1> array_float({size_dim1}, value_float);
    mio::CustomIndexArray<float, Dim1> array_float2 = array_float.convert<double>().convert<float>();
    for (auto i = mio::Index<Dim1>(0); i < size_dim1; ++i) {
        EXPECT_FLOAT_EQ(array_float2[i], value_float);
    }

    // case: double arrays with mixed precision; expect the value to be truncated during conversion to float
    double value_double = 1.0 + 5e-16;
    mio::CustomIndexArray<double, Dim1> array_double({size_dim1}, value_double);
    array_float = array_double.convert<float>();
    for (auto i = mio::Index<Dim1>(0); i < size_dim1; ++i) {
        // Check if array_double is unchanged
        EXPECT_DOUBLE_EQ(array_double[i], value_double);
        EXPECT_FLOAT_EQ(array_float[i], 1.f);
    }
}
