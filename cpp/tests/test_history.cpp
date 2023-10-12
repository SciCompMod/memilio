#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memilio/io/history.h"

class example
{
public:
    int a = 1;
    int b = 2;
};

struct LogPair : mio::LogAlways {
    using Type = std::pair<int, int>;
    static Type log(const example& ex)
    {
        return {ex.a, ex.b};
    }
};
struct LogAOnce : mio::LogOnce {
    using Type = int;
    static Type log(const example& ex)
    {
        return ex.a;
    }
};

template <class... Loggers>
struct Dummy {
    template <class T>
    size_t type_index()
    {
        return mio::details::index_templ_pack<T, Loggers...>();
    }
};

TEST(HistoryObject, log)
{
    example ex;
    mio::HistoryWithMemoryWriter<LogPair, LogAOnce> history;
    int n_runs = 2;
    for (int i = 0; i < n_runs; i++) {
        history.log(ex);
    }
    auto data = history.get_log();
    ASSERT_EQ(std::get<0>(data).size(), n_runs);
    ASSERT_EQ(std::get<1>(data).size(), 1); //cause we log this only once
    ASSERT_EQ(std::get<0>(data)[0], std::get<0>(data)[1]);
    ASSERT_EQ(std::get<0>(data)[0], std::make_pair(ex.a, ex.b));
    ASSERT_EQ(std::get<1>(data)[0], ex.a);
}
