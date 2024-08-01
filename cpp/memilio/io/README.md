# MEmilio C++ IO

This directory contains utilities for reading and writing data from and to files in different formats. The main part is a serialization framework that can be used to define the structure of data without using a specific file format. There are implementations of the framework for different formats. The framework is described in detail below, also see the [serialization example](../../examples/serialize.cpp).

## The Serialization framework

## Main functions and types

- functions serialize and deserialize:
     Main entry points to the framework to write and read values, respectively. The functions expect an IOContext
     (see Concepts below) that stores the serialized data. (De-)serialization can be customized by providing a
     (de-)serialize_internal overload or a (de-)serialize member function for the type. See the section "Implementing
     serialization for a new type" or the documentation for `serialize` and `deserialize`.
- IOStatus and IOResult:
     Used for error handling, see section "Error Handling" below.

## Concepts

1. IOContext
Stores data that describes serialized objects of any type in some unspecified format and provides structured
access to the data for deserialization. Implementations of this concept may store the data in any format
they want including binary. The data may also be written directly to disk. The context also keeps track
of errors. An IOContext object `io` allows the following operations:

    - `io.create_object("Type")`:
        Returns an IOObject for the type called `"Type"`. The IOObject (see below) allows adding data that describes
        the object to be serialized. The function must return something that can be assigned to a local
        variable, e.g., a temporary or copyable function. IOObject may store references to the context internally,
        so the lifetime of the local IOObject may not exceed the lifetime of the IOContext that created it.
    - `io.expect_object("Type")`:
        Returns an IOObject for the type called `"Type"`.
        The IOObject (see below) provides access to the data needed for deserialization.
    - `io.flags()`:
        Returns the flags that determine the behavior of serialization; see IOFlags.
    - `io.error()`:
        Returns an IOStatus object to check if there were any errors during serialization.
        Usually it is not necessary to check this manually but can be used to report the error faster and
        avoid expensive operations that would be wasted anyway.
    - `io.set_error(s)` with some IOStatus object:
        Stores an error that was generated outside of the IOContext, e.g., if a value that was deserialized
        is outside an allowed range.

1. IOObject
Gives structured access to serialized data. During serialization, data can be added with `add_...` operations.
During deserialization, data can be retrieved with `expect_...` operations. Data must be retrieved in the same order
as it was added since, e.g., binary format does not allow lookup by key. The following operations are supported
for an IOObject `obj`:

     - `obj.add_element("Name", t)`:
          Stores an object `t` in the IOObject under the key "Name". If `t` is of basic type (i.e., int, string), IOObject
          is expected to handle it directly. Otherwise, the object uses `mio::serialize` to get the data for `t`.
     - `obj.add_list("Name", b, e)`:
          Stores the elements in the range represented by iterators `b` and `e` under the key "Name". The individual elements are not named.
          The elements are either handled directly by the IOObject or using `mio::serialize` just like `add_element`.
     - `obj.add_optional("Name", p)`:
          Stores the element pointed to by pointer `p` under the key "Name". The pointer may be null. Otherwise identical to add_element.
     - `obj.expect_element("Name", Tag<T>{})`:
          If an object of type T can be found under the key "Name" and can be deserialized, returns the object. Otherwise returns
          an error. Analogously to serialization, the IOObject is expected to handle basic types directly and use `mio::deserialize`
          otherwise.
     - `obj.expect_list("Name", Tag<T>{})`:
          If a list of objects of type T can be found under the key "Name" and can be deserialized, returns a range that can be
          iterated over. Otherwise returns an error.
     - `obj.expect_optional("Name", Tag<T>{})`:
          Returns ```boost::optional<T>``` if an optional value of type T can be found under the key "Name". The optional may contain a
          value or it may be empty. Otherwise returns an error. Note that for some formats a wrong key is indistinguishable from
          an empty optional, so make sure to provide the correct key.

## Error handling

Errors are handled by returning error codes. The type IOStatus contains an error code and an optional string with additional
information. The type IOResult contains either a value or an IOStatus that describes an error. Operations that can fail return
an `IOResult<T>` where T is the type of the value that is produced by the operation if it is succesful. Except where necessary
because of dependencies, the framework does not throw nor catch any exceptions. IOContext and IOObject implementations are
expected to store errors. During serialization, `add_...` operations fail without returning errors, but the error is stored
in the IOObject and subsequent calls are usually no-ops. During deserialization, the values produced must usually be used or
inspected, so `expect_...` operations return an IOResult. The `apply` utility function provides a simple way to inspect the result
of multiple `expect_...` operations and use the values if all are succesful. See the documentation of `IOStatus`, `IOResult`
and `apply` below for more details.

## Adding a new data type to be serialized

Serialization of a new type T can be customized by providing _either_ member functions `serialize` and `deserialize` _or_ free functions
`serialize_internal` and `deserialize_internal`.

The `void serialize(IOContext& io)` member function takes an IO context and uses `create_object` and `add_...` operations
 to add data. The static `IOResult<T> deserialize(IOContext& io)` member function takes an IO context and uses `expect_...`
operations to retrieve the data. The `apply` utility function can be used to inspect the result of the `expect_...`
operations and construct the object of type T.
E.g.:

```cpp
struct Foo {
  int i;
  template<class IOContext>
  void serialize(IOContext& io) {
    auto obj = io.create_object("Foo");
    obj.add_element("i", i);
  }
  template<class IOContext>
  static IOResult<Foo> deserialize(IOContext& io) {
    auto obj = io.expect_object("Foo");
    auto i_result = obj.expect_element("i", mio::Tag<int>{});
    return mio::apply(io, [](auto&& i) { return Foo{i}; }, i_result);
  }
};
```

The free functions `serialize_internal` and `deserialize_internal` must be found with argument dependent lookup (ADL).
They can be used if no member function should or can be added to the type. See the code below for examples where
this was done for, e.g., Eigen3 matrices and STL containers.

## Adding a new format

Implement concepts IOContext and IOObject that provide the operations listed above. Your implemenation should handle
all built in types as well as std::string. It may handle other types (e.g., STL containers) as well if it can do so
more efficiently than the provided general free functions.

## Other IO modules

- HDF5 support classes for C++
- Reading of mobility matrix files

## Auto-serialization

This feature provides an easy and convenient method to serialize and deserialize classes, but with additional requirements and a reduced feature set. To give an example:

```cpp
struct Foo {
  int i;
  auto auto_serialize() {
    return make_auto_serialization("Foo", NVP("i", i));
  }
};
```

The auto-serialization effectively only supports the `add_element` and `expect_element` operations defined in the Concepts section, where the function arguments are provided by the name-value pairs (NVPs). Note that the value part of an NVP is also used to assign a value during deserialization, hence the class members must be used directly in the NVP constructor (i.e. as a non-const lvalue reference).

The requirements for auto-serialization are:
- The class must be trivially constructible.
  - Alternatively, you may provide a spezialisation of the struct `AutoSerializableFactory`.
- There is exactly one NVP for every class member (but the names and their order is arbitrary).
  - Values must be passed directly.
- Every class member itself is both (auto-)serializable and assignable.

## The command line interface

We provide a function `mio::command_line_interface` in the header `memilio/io/cli.h`, that can be used to write to or read from a parameter set. It can take parameters from command line arguments (i.e. the content of `argv` in the main function), and assign them to or get them from a `mio::ParameterSet`. A small example can be seen in `cpp/examples/cli.cpp`.

The command line interface (CLI) provides some non-parameter options listed below.

Name  (Alias) <div style="width:150px"></div>| Description
----|----
`--help` (`-h`) | Shows the basic usage of the CLI, and lists each parameter by name, as well as any alias and description. Takes priority before all other options and exits the programm.
`--print_option` | Can be used with a (space separated) list of parameter names or aliases (without dashes) to print the current values of each parameter to the terminal. This shows the correct json format used by the parameters. Exits after use.
`--read_from_json` | Allows reading parameters from a file instead of the command line. Both parameter names and aliases can be used, for example: <br>```{"<ParameterName>" : <value>, "<ParameterAlias>" : <value> }```
 `--write_to_json` | Writes _all_ parameters with their current values to a specified file.

In general, an option is defined as a string, which consists either of two dashes followed by a name (e.g. --help), or a single dash followed by an alias (e.g. -h). Apart from the built-in options, the names each refer to a parameter that can be set.

To set the value of a parameter from the command line, first type the corresponding parameter option (see --help), followed by the value that should be assigned (reference --print_option). Values are given as a Json value corresponding to the Type of the parameter. Note that some characters may need to be escaped or quoted. For example, the Json string `"some string"` must be entered as `\"some string\"` or `'"some string"'`.

## Working with the History object

The History object provides a way to save data throughout the simulation process. It offers an interface where users can define the data to be saved from a given object using Loggers and the method of saving it using Writers. Afterward, the user can access this data from the History object and manipulate it. For a basic Logger use case, refer to [this example](../../examples/history.cpp). For an example demonstrating using a Logger in the ABM, refer to [this example](../../examples/abm_history_example.cpp).

### Loggers

The `Logger` struct is a tool for logging data from a given object. Each user-implemented Logger must have a `Type` and implement two functions: `Type log(const T&)` and `bool should_log(const T&)`. The input `T` for these functions is the same as the one given to the `History` member-function `History::log`, e.g. `Model&` in the ABM.

- `Type`: Return Type of `log`.
- `log`: This function determines which data from the input `T` is saved. It must have the same return Type `Type` as the Loggers Type `Type`.
- `should_log`: This function must return a boolean to determine if data should be logged and can use the input `T` for this, e.g. if `T` fullfills some criteria.

Users can derive their Loggers from `LogOnce` or `LogAlways` to use a predefined `should_log` function. `LogOnce` logs only at the first call of `Logger::log()`, while `LogAlways` logs every time `log` is called. All implemented Loggers must be default constructible/destructible. For user-defined examples in the ABM, refer to [this file](../../models/abm/common_abm_loggers.h).

```cpp
struct LoggerExample{ /* :LogOnce/LogAlways if one wants to derive the should log from these. */
    using Type = /* type of the record */;
    /* Below, T must be replaced by the type T from History::log(t). */
    Type log(const T& t) 
    {
        return /* something of type Type */;
    }
    bool should_log(const T& t) 
    {
          /* Determine wether log and add_record should be called by History::log(t). */
          return /* true or false */
    }
};
```

### Writers

The `Writer` struct defines how to store the logged data from one or more implemented `Loggers`. Each user-implemented `Writer` must have a `Data` Type and implement the `template <class Logger> static void add_record(const typename Logger::Type& t, Data& data)` function.

- `Data`: This is some kind of container that stores the data returned by the Loggers. For example, this can be a TimeSeries or depend on the Loggers (like `std::tuple<std::vector<Logger::Type>...>`).
- `add_record`. This manipulates the passed Data member of the `History` class to store the value `t` returned by the Loggers. It is used whenever `History::log` is called and `Logger::should_log` is true.

A predefined universal `Writer` `DataWriterToMemory` is already implemented in [history.h](history.h). This stores the data from the loggers in a tuple of vectors every time the Logger is called. Another `Writer`  `TimeSeriesWriter` can be found in [this file](../../models/abm/common_abm_loggers.h), which saves data in a Timeseries. The according Logger has to have a suitable return type.

```cpp
template <class... Loggers>
struct DataWriterExample {
    using Data = /* Container for the stored data of the Loggers */;
    template <class Logger>
    static void add_record(const typename Logger::Type& t, Data& data)
    {
          /* Manipulation of data to store the value t returned by the Loggers */;
    }
};
```

### History

The `History` class manages the Writers and Loggers and provides an interface to log data. It is templated on one `Writer` and several suitable and unique `Loggers`. To use the Writer to log something, the `History` provides the function `void log(const T& t)` to call the `add_record` function of the `Writer` if the Logger function `should_log` returns true.

To access the data from the `History` class after logging, we provide the function `get_log` to access all records. For this, the lifetime of the `History` has to be as long as one wants to have access to the data, e.g. a history should not be constructed in the function it is called in when data is needed later.

To access data from a specific Logger, one can use `std::get<x>` where x is the position of the Logger in the template argument list of the `History` object. Refer to [this example](../../examples/history.cpp) for a simple implementation of a history object and [this full ABM example](../../simulation/abm.cpp) for a more advanced use case of the History object with several History objects in use.

As mentioned, if multiple Writers have to be used simultaneously, a separate History object is needed for each Writer. For a use case of this, refer to [the ABM Simulation advance function.](../../models/abm/simulation.cpp)
