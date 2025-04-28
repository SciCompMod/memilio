/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
/**
 * This example demonstrates using the serialization framework
 * and extending it for your own types.
 * See memilio/io/README.md for more information.
 */

#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"

namespace ioex
{
struct Foo {
    std::string s;

    //serialize Foo
    //the IOContext knows the format and can handle errors, all it needs is the data
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        //create an object to receive the data for this class
        auto obj = io.create_object("Foo");
        //add s to the object with a key so it can found when deserializing
        //the framework knows how to handle strings by itself
        obj.add_element("s", s);
    }

    //deserialize Foo
    //the IOContext provides the data
    template <class IOContext>
    static mio::IOResult<Foo> deserialize(IOContext& io)
    {
        // Retrieve an object for this class.
        //not all formats use the type string, but some (e.g. XML) require it to tag
        //elements that don't have a specific name, e.g. in a list.
        auto obj = io.expect_object("Foo");
        // Retrieve the data element by key. The tag defines the type of the element.
        auto s_rslt = obj.expect_element("s", mio::Tag<std::string>{});
        // The retrieval of elements may fail if the key cannot be found or the element cannot be
        //converted to the right type. use apply to inspect one or more results
        //and create the object if all lookups were succesful.
        return mio::apply(
            io,
            [](auto&& s_) {
                return Foo{s_};
            },
            s_rslt);
    }
};

struct Bar {
    int i;
    std::vector<Foo> foos;

    //serialize Bar
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Bar");
        obj.add_element("i", i);
        //not all data elements have a name of their own,
        //some are part of a container. Use add_list to add
        //multiple elements of the same type under the same key.
        //the framework uses Foo::serialize internally to
        //serialize the Foo objects in the list.
        obj.add_list("foos", foos.begin(), foos.end());
    }

    //deserialize Bar
    template <class IOContext>
    static mio::IOResult<Bar> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Bar");
        //lookup data elements in the same order as they were added by serialize
        //some formats (e.g. binary) don't support random access lookup.
        auto i_rslt    = obj.expect_element("i", mio::Tag<int>{});
        auto foos_rslt = obj.expect_list("foos", mio::Tag<Foo>{});
        //the function passed to apply is allowed to do more than just create the object.
        //e.g. it can validate values and return an error
        return mio::apply(
            io,
            [](auto&& i_, auto&& foos_) -> mio::IOResult<Bar> {
                //use mio::success or mio::failure to return an IOResult
                if (i_ >= 0) {
                    return mio::success(Bar{i_, std::vector<Foo>{foos_.begin(), foos_.end()}});
                }
                return mio::failure(mio::StatusCode::OutOfRange, "i must be non-negative.");
            },
            i_rslt, foos_rslt);
    }
};
} // namespace ioex

mio::IOResult<void> print_json()
{
    ioex::Bar b{42, {{"Hello"}, {"World"}}};

    //Try to turn the Bar object into a json value.
    auto rslt = mio::serialize_json(b);

    //IOResult can be inspected manually e.g.
    //if (rslt) { do_something(rslt.value()); return success(); }
    //else { return rslt.as_failure(); }
    //For convenience, the BOOST_OUTCOME_TRY macro can be used.
    //If the operation failed, the error is returned immediately.
    //If the operation was succesful, the result is unpacked and assigned to a new variable.
    //e.g.
    BOOST_OUTCOME_TRY(auto&& js, rslt);
    //could also be BOOST_OUTCOME_TRY(auto&& js, mio::serialize_json(b)) in one line

    //print json (Json::Value) to console
    //could also write to file or do anything else.
    std::cout << js << std::endl;

    //operation succesful, return void
    return mio::success();
}

mio::IOResult<ioex::Bar> read_json()
{
    //create json to deserialize
    //could also be read from file or stream
    Json::Value js;
    js["i"]            = 42;
    js["foos"][0]["s"] = "Hello";
    js["foos"][1]["s"] = "World";

    return mio::deserialize_json(js, mio::Tag<ioex::Bar>{});
}

int main()
{
    std::cout << "Printing ioex::Bar object...\n";
    auto r = print_json();
    if (r) {
        std::cout << "Success.\n";
    }
    else {
        std::cout << "Error: " << r.error().formatted_message() << "\n";
    }

    std::cout << "Deserializing ioex::Bar object... \n";
    auto bar_rslt = read_json();
    if (bar_rslt) {
        std::cout << "Success.\n";
    }
    else {
        std::cout << "Error: " << r.error().formatted_message() << "\n";
    }
}
