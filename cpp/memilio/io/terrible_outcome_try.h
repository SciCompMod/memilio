#ifndef TERRIBLE_OUTCOME_TRY_HACK
#define TERRIBLE_OUTCOME_TRY_HACK

// makro selector, aka pseudo overloads by number of args
// works by passing vaargs before name, shifting the selection of NAME to the left with additional args
#define GET_MACRO(_0, _1, _2, NAME, ...) NAME
// num args specializations
#define BOOST_OUTCOME_TRY_1(expr) if (!expr) return mio::failure(expr.error().code(), expr.error().message());
#define BOOST_OUTCOME_TRY_2(val, expr) BOOST_OUTCOME_TRY_1(expr) val = std::move(expr.value());
// main makro calling selector
#define BOOST_OUTCOME_TRY(...) GET_MACRO(_0, ##__VA_ARGS__, BOOST_OUTCOME_TRY_2, BOOST_OUTCOME_TRY_1, INVALID_0)(__VA_ARGS__)

#endif
