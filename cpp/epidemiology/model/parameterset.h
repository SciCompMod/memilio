#ifndef PARAMETERSET_H
#define PARAMETERSET_H

namespace epi
{

//TODO: Implement this as a compile-time map.
// I don't know how to implement compile time maps... :(
template <class... Params>
class ParameterSet
{
public:
    ParameterSet()
    {
    }

    template <class Parameter>
    void set(ScalarType value)
    {
        //TODO!!!
    }

    template <class Parameter>
    typename Parameter::Type const get() const
    {
        //TODO!!!
    }
};

} // namespace epi

#endif // COMPARTMENTALMODEL_H
