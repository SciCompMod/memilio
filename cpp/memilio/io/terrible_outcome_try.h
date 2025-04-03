#ifndef TERRIBLE_OUTCOME_TRY_HACK
#define TERRIBLE_OUTCOME_TRY_HACK

template <class T>
auto ___BOT_eval(T&& expr) {
    return (std::decay_t<decltype(expr.value())>) expr.value(); 
}

// makro selector, aka pseudo overloads by number of args
// works by passing vaargs before name, shifting the selection of NAME to the left with additional args
#define GET_MACRO(_0, _1, _2, NAME, ...) NAME
// num args specializations
#define BOOST_OUTCOME_TRY_1(expr) if (!expr) return mio::failure(expr.error().code(), expr.error().message());
#define BOOST_OUTCOME_TRY_2(val, expr) BOOST_OUTCOME_TRY_1(expr) val = ___BOT_eval(expr);
// main makro calling selector
#define BOOST_OUTCOME_TRY(...) GET_MACRO(_0, ##__VA_ARGS__, BOOST_OUTCOME_TRY_2, BOOST_OUTCOME_TRY_1, INVALID_0)(__VA_ARGS__)

#endif
