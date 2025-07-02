/*---------------------------------------------------------------------------*\
 * Algorithmic Differentiation by Overloading in C++
    Copyright (C) 2021 Software and Tools for Computational Engineering (STCE)
                  RWTH Aachen University               www.stce.rwth-aachen.de
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
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

#pragma once

#ifndef AD_HPP
#define AD_HPP

#define AD_VERSION_SHORT "AD v1.0"

#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <complex>
#include <string>
#include <stack>
#include <exception>
#include <stdexcept>
#include <string>
#include <bitset>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <string.h>

#include <cassert>
#include <limits>
#include <iomanip>

#if !defined(_WIN32)
#include <unistd.h>
#include <sys/mman.h>
#include <sys/time.h>
#endif

#if !defined(AD_DOXYGEN) & defined(_WIN32)

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <windows.h>
#include <algorithm>

#endif

// for stce::timer

#include <fcntl.h>
#include <typeinfo>

// define all DEPRECATED, to generate compiler-warnings,
// if a marked function is instantiated
// (only implemented for gnu compiler and visual C++)
#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#else
// #pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

#ifndef AD_DEFAULT_TAPE_SIZE
#define AD_DEFAULT_TAPE_SIZE (1024*1024*200) // Increased from 10 MB to 200 MB
#endif

#define ADi_INLINEDEF inline

typedef long int AD_TAPE_INT;
#define ADi_BASE_TYPE double

#ifdef AD_DEBUG
# ifndef AD_TAPE_BOUNDS_CHECK
#  define AD_TAPE_BOUNDS_CHECK
# endif
# define CHECK_OVERFLOW(x, a) assert(x < static_cast<AD_TAPE_INT>(std::numeric_limits<AD_TAPE_INT>::max()-a))
#else
# define CHECK_OVERFLOW(x, a)
#endif

// this macro enables a constructor to the type passed as argument
namespace ad {
template <typename T,typename T2=T>
  struct ad_type_constructable_from { const static bool value = false; };
}

#define AD_ENABLE_TYPE_CONSTRUCTION_FROM(T)    \
  namespace ad {                          \
  template <typename T2>                        \
  struct ad_type_constructable_from<T,T2> {    \
    const static bool value = true;             \
    typedef T2 type;                            \
  };                                            \
  }

typedef long int AD_TAPE_INT;
namespace ad {
    namespace operations {
        template<class AD_TAPE_REAL>struct ad_add_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval (const T1 &arg1, const T2 &arg2) {
                return arg1._value() + arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_sub_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval (const T1 &arg1, const T2 &arg2) {
                return arg1._value() - arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return -1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_mul_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval (const T1 &arg1, const T2 &arg2) {
                return arg1._value() * arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return arg1._value();
            }
        };
        template<class AD_TAPE_REAL>struct ad_div_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval (const T1 &arg1, const T2 &arg2) {
                return arg1._value() / arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return 1.0 / arg2._value();
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const T1 &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return -_value / arg2._value();
            }
        };
        template<class AD_TAPE_REAL>struct ad_add_ap {
            template<class T1>static inline const AD_TAPE_REAL eval(const T1 &arg1, const AD_TAPE_REAL &arg2) {
                return arg1._value() + arg2;
            } template<class T1>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const AD_TAPE_REAL &arg2 ) {
                (void)_value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_add_pa {
            template<class T2>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T2 &arg2) {
                return arg1 + arg2._value();
            } template<class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const AD_TAPE_REAL &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_sub_ap {
            template<class T1>static inline const AD_TAPE_REAL eval(const T1 &arg1, const AD_TAPE_REAL &arg2) {
                return arg1._value() - arg2;
            } template<class T1>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const AD_TAPE_REAL &arg2 ) {
                (void)_value;
                (void)arg1;
                (void)arg2;
                return 1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_sub_pa {
            template<class T2>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T2 &arg2) {
                return arg1 - arg2._value();
            } template<class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const AD_TAPE_REAL &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return -1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_mul_ap {
            template<class T1>static inline const AD_TAPE_REAL eval(const T1 &arg1, const AD_TAPE_REAL &arg2) {
                return arg1._value() * arg2;
            } template<class T1>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const AD_TAPE_REAL &arg2 ) {
                (void)_value;
                (void)arg1;
                (void)arg2;
                return arg2;
            }
        };
        template<class AD_TAPE_REAL>struct ad_mul_pa {
            template<class T2>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T2 &arg2) {
                return arg1 * arg2._value();
            } template<class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const AD_TAPE_REAL &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return arg1;
            }
        };
        template<class AD_TAPE_REAL>struct ad_div_ap {
            template<class T1>static inline const AD_TAPE_REAL eval(const T1 &arg1, const AD_TAPE_REAL &arg2) {
                return arg1._value() / arg2;
            } template<class T1>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL &_value, const T1 &arg1, const AD_TAPE_REAL &arg2 ) {
                (void)_value;
                (void)arg1;
                (void)arg2;
                return 1.0 / arg2;
            }
        };
        template<class AD_TAPE_REAL>struct ad_div_pa {
            template<class T2>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T2 &arg2) {
                return arg1 / arg2._value();
            } template<class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL &_value, const AD_TAPE_REAL &arg1, const T2 &arg2 ) {
                (void) _value;
                (void)arg1;
                (void)arg2;
                return -_value / arg2._value();
            }
        };
        using std::sin;
        template<class AD_TAPE_REAL>struct ad_sin {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return sin(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (cos(x._value()));
            }
        };
        using std::cos;
        template<class AD_TAPE_REAL>struct ad_cos {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return cos(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (-sin(x._value()));
            }
        };
        using std::tan;
        template<class AD_TAPE_REAL>struct ad_tan {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return tan(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return ((1.0 + (tan(x._value())*tan(x._value()))));
            }
        };
        using std::cosh;
        template<class AD_TAPE_REAL>struct ad_cosh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return cosh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (sinh(x._value()));
            }
        };
        using std::sinh;
        template<class AD_TAPE_REAL>struct ad_sinh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return sinh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (cosh(x._value()));
            }
        };
        using std::asin;
        template<class AD_TAPE_REAL>struct ad_asin {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return asin(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1 / sqrt(1.0 - x._value()*x._value()));
            }
        };
        using std::acos;
        template<class AD_TAPE_REAL>struct ad_acos {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return acos(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (-1 / sqrt(1.0 - x._value()*x._value()));
            }
        };
        using std::exp;
        template<class AD_TAPE_REAL>struct ad_exp {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return exp(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (exp(x._value()));
            }
        };
        using std::atan;
        template<class AD_TAPE_REAL>struct ad_atan {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return atan(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0 / (1.0 + x._value()*x._value()));
            }
        };
        using std::tanh;
        template<class AD_TAPE_REAL>struct ad_tanh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return tanh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0 - tanh(x._value())*tanh(x._value()));
            }
        };
        using std::sqrt;
        template<class AD_TAPE_REAL>struct ad_sqrt {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return sqrt(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0/(2.0*sqrt(x._value() + (x._value() > 0 ? 0.0 : 1e-8))));
            }
        };
        using std::log;
        template<class AD_TAPE_REAL>struct ad_log {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return log(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0 / x._value());
            }
        };
        using ::erf;
        template<class AD_TAPE_REAL>struct ad_erf {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return erf(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (2.0 / sqrt(3.14159265358979323846264338327950288) * exp(-x._value() * x._value()));
            }
        };
        using ::erfc;
        template<class AD_TAPE_REAL>struct ad_erfc {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return erfc(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (-2.0 / sqrt(3.14159265358979323846264338327950288) * exp(-x._value() * x._value()));
            }
        };
        using ::asinh;
        template<class AD_TAPE_REAL>struct ad_asinh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return asinh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1. / sqrt(1. + (x._value()*x._value())));
            }
        };
        using ::acosh;
        template<class AD_TAPE_REAL>struct ad_acosh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return acosh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1. / sqrt((x._value()*x._value()) - 1.));
            }
        };
        using ::expm1;
        template<class AD_TAPE_REAL>struct ad_expm1 {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return expm1(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (exp(x._value()));
            }
        };
        using ::atanh;
        template<class AD_TAPE_REAL>struct ad_atanh {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return atanh(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1. / (1. - (x._value()*x._value())));
            }
        };
        using ::log1p;
        template<class AD_TAPE_REAL>struct ad_log1p {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return log1p(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0 / (x._value() + 1));
            }
        };
        using ::log10;
        template<class AD_TAPE_REAL>struct ad_log10 {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg) {
                return log10(arg._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &x) {
                (void)_value;
                return (1.0 / (x._value()*log(10)));
            }
        };
        template<class AD_TAPE_REAL>struct ad_minus {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1) {
                return -arg1._value();
            }
            template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &arg1 ) {
                (void)_value;
                (void)arg1;
                return -1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_plus {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1) {
                return arg1._value();
            }
            template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &arg1 ) {
                (void)_value;
                (void)arg1;
                return 1.0;
            }
        };
        using ::fabs;
        template<class AD_TAPE_REAL>struct ad_fabs {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1) {
                return fabs(arg1._value());
            }
            template<class T>static inline const AD_TAPE_REAL calc_partial(const AD_TAPE_REAL &_value, const T &arg1 ) {
                (void) _value;
                if (arg1._value() < 0) return -1.0;
                else return 1.0;
            }
        };
        template<class AD_TAPE_REAL>struct ad_atan2_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval(const T1 &arg1, const T2 &arg2) {
                (void) arg1;
                (void) arg2;
                return atan2(arg1._value(),arg2._value());
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg2._value() / (arg2._value() * arg2._value() + arg1._value() * arg1._value());
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return -arg1._value() / (arg2._value() * arg2._value() + arg1._value() * arg1._value());
            }
        };
        template<class AD_TAPE_REAL>struct ad_atan2_ap {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1, const AD_TAPE_REAL &arg2) {
                (void) arg1;
                (void) arg2;
                return atan2(arg1._value(),arg2);
            } template<class T>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T &arg1, const AD_TAPE_REAL &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg2 / (arg2 * arg2 + arg1._value() * arg1._value());
            }
        };
        template<class AD_TAPE_REAL>struct ad_atan2_pa {
            template<class T>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T &arg2) {
                (void) arg1;
                (void) arg2;
                return atan2(arg1,arg2._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const AD_TAPE_REAL &arg1, const T &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return -arg1 / (arg2._value() * arg2._value() + arg1 *arg1);
            }
        };
        template<class AD_TAPE_REAL>
        struct ad_pow_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval(const T1 &arg1, const T2 &arg2) {
                return pow(arg1._value(), arg2._value());
            }
            template<class T1, class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void)_value;
                return arg2._value() * pow(arg1._value(), arg2._value() - 1.0);
            }
            template<class T1, class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void) _value;
                if (arg1 <= 0)
                    return 0;
                else
                    return log(arg1._value()) * pow(arg1._value(), arg2._value());
            }
        };
        template<class AD_TAPE_REAL>
        struct ad_pow_ap {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1, const AD_TAPE_REAL &arg2) {
                return pow(arg1._value(), arg2);
            }
            template<class T>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T &arg1, const AD_TAPE_REAL &arg2) {
                (void) _value;
                return arg2 * pow(arg1._value(), arg2 - 1.0);
            }
        };
        template<class AD_TAPE_REAL>
        struct ad_pow_pa {
            template<class T>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T &arg2) {
                return pow(arg1, arg2._value());
            }
            template<class T>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const AD_TAPE_REAL &arg1, const T &arg2) {
                (void) _value;
                return log(arg1) * pow(arg1, arg2._value());
            }
        };
        template<class AD_TAPE_REAL>struct ad_hypot_aa {
            template<class T1, class T2>static inline const AD_TAPE_REAL eval(const T1 &arg1, const T2 &arg2) {
                (void) arg1;
                (void) arg2;
                return hypot(arg1._value(),arg2._value());
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg1._value() / _value;
            } template<class T1,class T2>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const T1 &arg1, const T2 &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg2._value() / _value;
            }
        };
        template<class AD_TAPE_REAL>struct ad_hypot_ap {
            template<class T>static inline const AD_TAPE_REAL eval(const T &arg1, const AD_TAPE_REAL &arg2) {
                (void) arg1;
                (void) arg2;
                return hypot(arg1._value(),arg2);
            } template<class T>static inline const AD_TAPE_REAL calc_partial1(const AD_TAPE_REAL _value, const T &arg1, const AD_TAPE_REAL &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg1._value() / _value;
            }
        };
        template<class AD_TAPE_REAL>struct ad_hypot_pa {
            template<class T>static inline const AD_TAPE_REAL eval(const AD_TAPE_REAL &arg1, const T &arg2) {
                (void) arg1;
                (void) arg2;
                return hypot(arg1,arg2._value());
            } template<class T>static inline const AD_TAPE_REAL calc_partial2(const AD_TAPE_REAL _value, const AD_TAPE_REAL &arg1, const T &arg2) {
                (void) _value;
                (void) arg1;
                (void) arg2;
                return arg2._value() / _value;
            }
        };
    }
}
namespace ad {
    namespace helper {
        template <typename T>
        struct type_identity {
            typedef T type;
        };
    }
}
namespace ad {
    template<typename TYPE>
    class reference_constructor_wrapper {
            const TYPE &T;
        public:
            reference_constructor_wrapper(const TYPE &T_) : T(T_) {}
            template <typename TO_CREATE>
            TO_CREATE *create() const {
                return new TO_CREATE(T);
            }
    };
    template<>
    class reference_constructor_wrapper<void *> {
        public:
            reference_constructor_wrapper(void *v) {
                (void) v;
            }
            template <typename TO_CREATE>
            TO_CREATE *create() const {
                return new TO_CREATE();
            }
    };
}
namespace ad {
    namespace internal {
        template<class AD_TAPE_REAL, class AD_ARG, class AD_OPERATION> struct unary_intermediate {
            const AD_TAPE_REAL _value_;
            const AD_ARG &_arg;
            typedef AD_TAPE_REAL VALUE_TYPE;
            typedef typename AD_ARG::DATA_TYPE DATA_TYPE;
            explicit unary_intermediate(const AD_ARG &arg) : _value_(AD_OPERATION::eval(arg)), _arg(arg) {
            }
            inline const AD_TAPE_REAL &_value() const {
                return _value_;
            }
            inline const AD_TAPE_REAL pval() const {
                return AD_OPERATION::calc_partial(_value(), _arg);
            }
        };
        template<class AD_TAPE_REAL, class AD_ARG1, class AD_ARG2, class AD_OPERATION> struct binary_intermediate_aa {
            AD_TAPE_REAL _value_;
            const AD_ARG1 &_arg1;
            const AD_ARG2 &_arg2;
            typedef AD_TAPE_REAL VALUE_TYPE;
            typedef typename AD_ARG1::DATA_TYPE DATA_TYPE;
            binary_intermediate_aa(const AD_ARG1 &arg1, const AD_ARG2 &arg2) :
                _value_(AD_OPERATION::eval(arg1, arg2)),
                _arg1(arg1),
                _arg2(arg2) {
            }
            inline const AD_TAPE_REAL pval1() const {
                return AD_OPERATION::calc_partial1(_value(), _arg1, _arg2);
            }
            inline const AD_TAPE_REAL pval2() const {
                return AD_OPERATION::calc_partial2(_value(), _arg1, _arg2);
            }
            inline const AD_TAPE_REAL &_value() const {
                return _value_;
            }
        };
        template<class AD_TAPE_REAL, class AD_ARG1, class AD_OPERATION> struct binary_intermediate_ap {
            const AD_TAPE_REAL _value_;
            const AD_ARG1 &_arg1;
            const AD_TAPE_REAL _arg2;
            typedef AD_TAPE_REAL VALUE_TYPE;
            typedef typename AD_ARG1::DATA_TYPE DATA_TYPE;
            binary_intermediate_ap(const AD_ARG1 &arg1, const AD_TAPE_REAL &arg2) :
                _value_(AD_OPERATION::eval(arg1, arg2)), _arg1(arg1), _arg2(arg2) {}
            inline const AD_TAPE_REAL &_value() const {
                return _value_;
            }
            inline const AD_TAPE_REAL pval1() const {
                return AD_OPERATION::calc_partial1(_value_, _arg1, _arg2);
            }
        };
        template<class AD_TAPE_REAL, class AD_ARG2, class AD_OPERATION> struct binary_intermediate_pa {
            const AD_TAPE_REAL _value_;
            const AD_TAPE_REAL _arg1;
            const AD_ARG2 &_arg2;
            typedef AD_TAPE_REAL VALUE_TYPE;
            typedef typename AD_ARG2::DATA_TYPE DATA_TYPE;
            binary_intermediate_pa(const AD_TAPE_REAL &arg1, const AD_ARG2 &arg2) :
                _value_(AD_OPERATION::eval(arg1, arg2)), _arg1(arg1), _arg2(arg2) {}
            inline const AD_TAPE_REAL &_value() const {
                return _value_;
            }
            inline const AD_TAPE_REAL pval2() const {
                return AD_OPERATION::calc_partial2(_value_, _arg1, _arg2);
            }
        };
        template <typename T>
        struct passive_value_type_of {
            typedef T TYPE;
        };
        template<class AD_TAPE_REAL, class DATA_HANDLER> struct active_type {
            private:
                AD_TAPE_REAL _value_;
                DATA_HANDLER _data_;
            public:
                typedef AD_TAPE_REAL VALUE_TYPE;
                typedef DATA_HANDLER DATA_TYPE;
                typedef typename passive_value_type_of<AD_TAPE_REAL>::TYPE PASSIVE_VALUE_TYPE;
                inline const AD_TAPE_REAL &_value() const {
                    return _value_;
                }
                inline AD_TAPE_REAL &_value() {
                    return _value_;
                }
                inline const DATA_HANDLER &_data() const {
                    return _data_;
                }
                inline DATA_HANDLER &_data() {
                    return _data_;
                }
                inline void _clear() {
                    _data_.clear();
                }

                inline active_type(const active_type&) = default;
                template<typename TYPE, typename TEST = typename ad_type_constructable_from<TYPE,active_type>::type>
                inline active_type& operator=(const TYPE &val) {
                    _value_ = val;
                    _data_.clear();
                    return *this;
                }
                template<typename TYPE, typename TEST = typename ad_type_constructable_from<TYPE,active_type>::type>
                inline active_type(const TYPE &val) :
                    _value_(val) {}

                inline active_type() : _value_(static_cast<AD_TAPE_REAL>(0.0)) {}
                template<class AD_TAPE_REAL_TMP, class DATA_HANDLER_TMP>
                inline active_type(const active_type<AD_TAPE_REAL_TMP, DATA_HANDLER_TMP> &val) : _value_(val) {}
                inline active_type(const PASSIVE_VALUE_TYPE &val) : _value_(val) {}
                inline active_type &operator = (const active_type &x) {
                    this->_value_ = x._value_;
                    this->_data_ = x._data_;
                    return *this;
                }
                inline active_type &operator =(const PASSIVE_VALUE_TYPE &vneu) {
                    this->_value_ = vneu;
                    this->_data_.clear();
                    return *this;
                }
                template<class AD_TAPE_REAL_TMP, class DATA_HANDLER_TMP>
                inline active_type &operator =(const active_type<AD_TAPE_REAL_TMP, DATA_HANDLER_TMP> &vneu) {
                    this->_value_ = vneu;
                    this->_data_.clear();
                    return *this;
                }
            private:
                template<class A1_T1, class A1_T2, class A1_OP > inline void build_from(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    DATA_HANDLER::handle(x, *this);
                    this->_value_ = x._value_;
                } public:
                template<class A1_T1, class A1_T2, class A1_OP > active_type(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    build_from(x);
                } template<class A1_T1, class A1_T2, class A1_OP > inline active_type& operator=(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    build_from(x);
                    return *this;
                }
            private:
                template<class A1_T1, class A1_OP > inline void build_from(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    DATA_HANDLER::handle(x, *this);
                    this->_value_ = x._value_;
                } public:
                template<class A1_T1, class A1_OP > active_type(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    build_from(x);
                } template<class A1_T1, class A1_OP > inline active_type& operator=(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    build_from(x);
                    return *this;
                }
            private:
                template<class A1_T2, class A1_OP > inline void build_from(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    DATA_HANDLER::handle(x, *this);
                    this->_value_ = x._value_;
                } public:
                template<class A1_T2, class A1_OP > active_type(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    build_from(x);
                } template<class A1_T2, class A1_OP > inline active_type& operator=(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    build_from(x);
                    return *this;
                }
            private:
                template<class A1_T, class A1_OP > inline void build_from(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    DATA_HANDLER::handle(x, *this);
                    this->_value_ = x._value_;
                } public:
                template<class A1_T, class A1_OP > active_type(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    build_from(x);
                } template<class A1_T, class A1_OP > inline active_type& operator=(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    build_from(x);
                    return *this;
                }
                template<class DATA_HANDLER_TMP> inline active_type& operator += (const active_type<AD_TAPE_REAL, DATA_HANDLER_TMP> &x) {
                    *this = *this + x;
                    return *this;
                } template<class A1_T1, class A1_T2, class A1_OP > inline active_type& operator += (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    *this = *this + x;
                    return *this;
                } template<class A1_T1, class A1_OP > inline active_type& operator += (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    *this = *this + x;
                    return *this;
                } template<class A1_T2, class A1_OP > inline active_type& operator += (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    *this = *this + x;
                    return *this;
                } template<class A1_T, class A1_OP > inline active_type& operator += (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    *this = *this + x;
                    return *this;
                } inline active_type& operator += (const AD_TAPE_REAL &x) {
                    this->_value() += x;
                    return *this;
                }
                template<class DATA_HANDLER_TMP> inline active_type& operator -= (const active_type<AD_TAPE_REAL, DATA_HANDLER_TMP> &x) {
                    *this = *this - x;
                    return *this;
                } template<class A1_T1, class A1_T2, class A1_OP > inline active_type& operator -= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    *this = *this - x;
                    return *this;
                } template<class A1_T1, class A1_OP > inline active_type& operator -= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    *this = *this - x;
                    return *this;
                } template<class A1_T2, class A1_OP > inline active_type& operator -= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    *this = *this - x;
                    return *this;
                } template<class A1_T, class A1_OP > inline active_type& operator -= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    *this = *this - x;
                    return *this;
                } inline active_type& operator -= (const AD_TAPE_REAL &x) {
                    this->_value() -= x;
                    return *this;
                }
                template<class DATA_HANDLER_TMP> inline active_type& operator *= (const active_type<AD_TAPE_REAL, DATA_HANDLER_TMP> &x) {
                    *this = *this * x;
                    return *this;
                } template<class A1_T1, class A1_T2, class A1_OP > inline active_type& operator *= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    *this = *this * x;
                    return *this;
                } template<class A1_T1, class A1_OP > inline active_type& operator *= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    *this = *this * x;
                    return *this;
                } template<class A1_T2, class A1_OP > inline active_type& operator *= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    *this = *this * x;
                    return *this;
                } template<class A1_T, class A1_OP > inline active_type& operator *= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    *this = *this * x;
                    return *this;
                } inline active_type& operator *= (const AD_TAPE_REAL &x) {
                    *this = *this * x;
                    return *this;
                }
                template<class DATA_HANDLER_TMP> inline active_type& operator /= (const active_type<AD_TAPE_REAL, DATA_HANDLER_TMP> &x) {
                    *this = *this / x;
                    return *this;
                } template<class A1_T1, class A1_T2, class A1_OP > inline active_type& operator /= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                    *this = *this / x;
                    return *this;
                } template<class A1_T1, class A1_OP > inline active_type& operator /= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                    *this = *this / x;
                    return *this;
                } template<class A1_T2, class A1_OP > inline active_type& operator /= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                    *this = *this / x;
                    return *this;
                } template<class A1_T, class A1_OP > inline active_type& operator /= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                    *this = *this / x;
                    return *this;
                } inline active_type& operator /= (const AD_TAPE_REAL &x) {
                    *this = *this / x;
                    return *this;
                }
                inline active_type &operator ++() {
                    ++this->_value_;
                    return *this;
                }
                inline active_type &operator --() {
                    --this->_value_;
                    return *this;
                }
                inline active_type operator ++(int) {
                    active_type ret(*this);
                    ++this->_value_;
                    return ret;
                }
                inline active_type operator --(int) {
                    active_type ret(*this);
                    --this->_value_;
                    return ret;
                }
        };
        template<class AD_TAPE_REAL, class DATA_HANDLER>
        struct passive_value_type_of<active_type<AD_TAPE_REAL, DATA_HANDLER> > {
            typedef typename passive_value_type_of<AD_TAPE_REAL>::TYPE TYPE;
        };

        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_minus<AD_TAPE_REAL> > operator -( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_minus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_minus<AD_TAPE_REAL> > operator -( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_minus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_minus<AD_TAPE_REAL> > operator -( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_minus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_minus<AD_TAPE_REAL> > operator -( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_minus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_minus<AD_TAPE_REAL> > operator -( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_minus<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_plus<AD_TAPE_REAL> > operator +( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_plus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_plus<AD_TAPE_REAL> > operator +( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_plus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_plus<AD_TAPE_REAL> > operator +( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_plus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_plus<AD_TAPE_REAL> > operator +( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_plus<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_plus<AD_TAPE_REAL> > operator +( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_plus<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sin<AD_TAPE_REAL> > sin( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_sin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sin<AD_TAPE_REAL> > sin( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_sin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sin<AD_TAPE_REAL> > sin( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_sin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sin<AD_TAPE_REAL> > sin( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_sin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sin<AD_TAPE_REAL> > sin( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_sin<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_cos<AD_TAPE_REAL> > cos( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_cos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_cos<AD_TAPE_REAL> > cos( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_cos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_cos<AD_TAPE_REAL> > cos( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_cos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_cos<AD_TAPE_REAL> > cos( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_cos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_cos<AD_TAPE_REAL> > cos( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_cos<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_tan<AD_TAPE_REAL> > tan( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_tan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_tan<AD_TAPE_REAL> > tan( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_tan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_tan<AD_TAPE_REAL> > tan( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_tan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_tan<AD_TAPE_REAL> > tan( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_tan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_tan<AD_TAPE_REAL> > tan( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_tan<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_cosh<AD_TAPE_REAL> > cosh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_cosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_cosh<AD_TAPE_REAL> > cosh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_cosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_cosh<AD_TAPE_REAL> > cosh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_cosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_cosh<AD_TAPE_REAL> > cosh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_cosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_cosh<AD_TAPE_REAL> > cosh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_cosh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sinh<AD_TAPE_REAL> > sinh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_sinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sinh<AD_TAPE_REAL> > sinh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_sinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sinh<AD_TAPE_REAL> > sinh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_sinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sinh<AD_TAPE_REAL> > sinh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_sinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sinh<AD_TAPE_REAL> > sinh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_sinh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_asin<AD_TAPE_REAL> > asin( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_asin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_asin<AD_TAPE_REAL> > asin( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_asin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_asin<AD_TAPE_REAL> > asin( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_asin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_asin<AD_TAPE_REAL> > asin( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_asin<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_asin<AD_TAPE_REAL> > asin( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_asin<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_acos<AD_TAPE_REAL> > acos( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_acos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_acos<AD_TAPE_REAL> > acos( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_acos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_acos<AD_TAPE_REAL> > acos( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_acos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_acos<AD_TAPE_REAL> > acos( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_acos<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_acos<AD_TAPE_REAL> > acos( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_acos<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_exp<AD_TAPE_REAL> > exp( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_exp<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_exp<AD_TAPE_REAL> > exp( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_exp<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_exp<AD_TAPE_REAL> > exp( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_exp<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_exp<AD_TAPE_REAL> > exp( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_exp<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_exp<AD_TAPE_REAL> > exp( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_exp<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan<AD_TAPE_REAL> > atan( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_atan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan<AD_TAPE_REAL> > atan( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_atan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan<AD_TAPE_REAL> > atan( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_atan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan<AD_TAPE_REAL> > atan( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_atan<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan<AD_TAPE_REAL> > atan( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_atan<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_tanh<AD_TAPE_REAL> > tanh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_tanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_tanh<AD_TAPE_REAL> > tanh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_tanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_tanh<AD_TAPE_REAL> > tanh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_tanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_tanh<AD_TAPE_REAL> > tanh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_tanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_tanh<AD_TAPE_REAL> > tanh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_tanh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sqrt<AD_TAPE_REAL> > sqrt( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_sqrt<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sqrt<AD_TAPE_REAL> > sqrt( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_sqrt<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sqrt<AD_TAPE_REAL> > sqrt( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_sqrt<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sqrt<AD_TAPE_REAL> > sqrt( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_sqrt<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sqrt<AD_TAPE_REAL> > sqrt( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_sqrt<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_log<AD_TAPE_REAL> > log( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_log<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_log<AD_TAPE_REAL> > log( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_log<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_log<AD_TAPE_REAL> > log( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_log<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_log<AD_TAPE_REAL> > log( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_log<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_log<AD_TAPE_REAL> > log( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_log<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_erf<AD_TAPE_REAL> > erf( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_erf<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_erf<AD_TAPE_REAL> > erf( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_erf<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_erf<AD_TAPE_REAL> > erf( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_erf<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_erf<AD_TAPE_REAL> > erf( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_erf<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_erf<AD_TAPE_REAL> > erf( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_erf<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_erfc<AD_TAPE_REAL> > erfc( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_erfc<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_erfc<AD_TAPE_REAL> > erfc( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_erfc<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_erfc<AD_TAPE_REAL> > erfc( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_erfc<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_erfc<AD_TAPE_REAL> > erfc( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_erfc<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_erfc<AD_TAPE_REAL> > erfc( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_erfc<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_expm1<AD_TAPE_REAL> > expm1( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_expm1<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_expm1<AD_TAPE_REAL> > expm1( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_expm1<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_expm1<AD_TAPE_REAL> > expm1( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_expm1<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_expm1<AD_TAPE_REAL> > expm1( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_expm1<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_expm1<AD_TAPE_REAL> > expm1( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_expm1<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_asinh<AD_TAPE_REAL> > asinh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_asinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_asinh<AD_TAPE_REAL> > asinh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_asinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_asinh<AD_TAPE_REAL> > asinh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_asinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_asinh<AD_TAPE_REAL> > asinh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_asinh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_asinh<AD_TAPE_REAL> > asinh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_asinh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_acosh<AD_TAPE_REAL> > acosh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_acosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_acosh<AD_TAPE_REAL> > acosh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_acosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_acosh<AD_TAPE_REAL> > acosh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_acosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_acosh<AD_TAPE_REAL> > acosh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_acosh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_acosh<AD_TAPE_REAL> > acosh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_acosh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atanh<AD_TAPE_REAL> > atanh( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_atanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atanh<AD_TAPE_REAL> > atanh( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_atanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atanh<AD_TAPE_REAL> > atanh( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_atanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atanh<AD_TAPE_REAL> > atanh( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_atanh<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atanh<AD_TAPE_REAL> > atanh( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_atanh<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_log1p<AD_TAPE_REAL> > log1p( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_log1p<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_log1p<AD_TAPE_REAL> > log1p( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_log1p<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_log1p<AD_TAPE_REAL> > log1p( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_log1p<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_log1p<AD_TAPE_REAL> > log1p( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_log1p<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_log1p<AD_TAPE_REAL> > log1p( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_log1p<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_log10<AD_TAPE_REAL> > log10( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_log10<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_log10<AD_TAPE_REAL> > log10( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_log10<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_log10<AD_TAPE_REAL> > log10( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_log10<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_log10<AD_TAPE_REAL> > log10( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_log10<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_log10<AD_TAPE_REAL> > log10( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_log10<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_fabs<AD_TAPE_REAL> > fabs( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > fabs( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > fabs( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > fabs( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > fabs( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_fabs<AD_TAPE_REAL> > abs( const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > abs( const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > abs( const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > abs( const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_fabs<AD_TAPE_REAL> > abs( const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1) {
            return ad::internal::unary_intermediate<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::operations::ad_fabs<AD_TAPE_REAL> >(x1);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_add_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_ap<AD_TAPE_REAL> > operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_ap<AD_TAPE_REAL> > >::type operator + (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_pa<AD_TAPE_REAL> > operator + (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_pa<AD_TAPE_REAL> > >::type operator + (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > >::type operator + (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > operator + (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > >::type operator + (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > >::type operator + (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > operator + (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > >::type operator + (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > >::type operator + (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > operator + (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > >::type operator + (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> > >::type operator + (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > operator + (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> > >::type operator + (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_add_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_sub_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > >::type operator - (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > operator - (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > >::type operator - (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > >::type operator - (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > operator - (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > >::type operator - (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > >::type operator - (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > operator - (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > >::type operator - (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > >::type operator - (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > operator - (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > >::type operator - (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> > >::type operator - (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > operator - (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> > >::type operator - (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_sub_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_mul_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > >::type operator * (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > operator * (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > >::type operator * (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > >::type operator * (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > operator * (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > >::type operator * (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > >::type operator * (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > operator * (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > >::type operator * (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > >::type operator * (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > operator * (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > >::type operator * (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> > >::type operator * (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > operator * (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> > >::type operator * (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_mul_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_div_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_ap<AD_TAPE_REAL> > operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_ap<AD_TAPE_REAL> > >::type operator / (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_pa<AD_TAPE_REAL> > operator / (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_pa<AD_TAPE_REAL> > >::type operator / (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > >::type operator / (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > operator / (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > >::type operator / (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > >::type operator / (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > operator / (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > >::type operator / (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > >::type operator / (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > operator / (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > >::type operator / (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> > >::type operator / (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > operator / (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> > >::type operator / (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_div_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_atan2_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > >::type atan2 (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > atan2 (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > >::type atan2 (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > >::type atan2 (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > atan2 (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > >::type atan2 (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > >::type atan2 (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > atan2 (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > >::type atan2 (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > >::type atan2 (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > atan2 (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > >::type atan2 (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> > >::type atan2 (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > atan2 (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> > >::type atan2 (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_atan2_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_pow_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > >::type pow (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > pow (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > >::type pow (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > >::type pow (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > pow (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > >::type pow (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > >::type pow (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > pow (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > >::type pow (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > >::type pow (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > pow (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > >::type pow (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> > >::type pow (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > pow (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> > >::type pow (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_pow_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &x2) {
            return ad::internal::binary_intermediate_aa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>,ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>, ad::operations::ad_hypot_aa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > >::type hypot (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x1, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > hypot (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > typename std::enable_if<!std::is_same<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > >::type hypot (const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::PASSIVE_VALUE_TYPE &x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > >::type hypot (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x1, const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > hypot (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > >::type hypot (const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > >::type hypot (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > hypot (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > >::type hypot (const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > >::type hypot (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x1, const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > hypot (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > >::type hypot (const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> > >::type hypot (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x1, const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x2) {
            return ad::internal::binary_intermediate_ap<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_ap<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > hypot (const typename ad::helper::type_identity<AD_TAPE_REAL>::type &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > typename std::enable_if<!std::is_same<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE>::value, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> > >::type hypot (const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::PASSIVE_VALUE_TYPE &x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x2) {
            return ad::internal::binary_intermediate_pa<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::operations::ad_hypot_pa<AD_TAPE_REAL> >(x1,x2);
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() == x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() == x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator == (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 == x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator == (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() == x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() == x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 == x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() == x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 == x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator == (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() == x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator == (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 == x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() != x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() != x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator != (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 != x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator != (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() != x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() != x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 != x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() != x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 != x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator != (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() != x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator != (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 != x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() < x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() < x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator < (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 < x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator < (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() < x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() < x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 < x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() < x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 < x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator < (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() < x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator < (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 < x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() <= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() <= x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator <= (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator <= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() <= x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() <= x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() <= x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 <= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator <= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() <= x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator <= (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 <= x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() > x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() > x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator > (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 > x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator > (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() > x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() > x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 > x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() > x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 > x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator > (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() > x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator > (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 > x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) ;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP>& x2) {
            return x1._value() >= x2._value();
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x1, const double& x2) {
            return x1._value() >= x2;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2);
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool operator >= (const double& x1, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x2) {
            return x1 >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator >= (const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x1, const double& x2) {
            return x1._value() >= x2;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x2) {
            return x1 >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() >= x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x2) {
            return x1 >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x1, const double& x2) {
            return x1._value() >= x2;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x2) {
            return x1 >= x2._value();
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator >= (const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x1, const double& x2) {
            return x1._value() >= x2;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2);
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool operator >= (const double& x1, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x2) {
            return x1 >= x2._value();
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 >
        static inline void reset_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
            double tmp = 0;
            get(x, tmp);
            x = tmp;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 >
        static inline std::istream &operator >> (std::istream &in, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
            AD_TAPE_REAL &tmp = const_cast<AD_TAPE_REAL &>(x._value());
            in >> tmp;
            return in;
        }
        using std::ceil;
        using std::floor;
        using std::isfinite;
        using std::isnan;
        using std::isinf;
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool isnan(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            return isnan(x._value());
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool isinf(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            return isinf(x._value());
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline double ceil(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            return ceil(x._value());
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline double floor(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            return floor(x._value());
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline bool isfinite(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            return isfinite(x._value());
        }
        using std::ceil;
        using std::floor;
        using std::isfinite;
        using std::isnan;
        using std::isinf;
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool isnan(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            return isnan(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool isinf(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            return isinf(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline double ceil(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            return ceil(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline double floor(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            return floor(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline bool isfinite(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            return isfinite(x._value());
        }
        using std::ceil;
        using std::floor;
        using std::isfinite;
        using std::isnan;
        using std::isinf;
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool isnan(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            return isnan(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool isinf(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            return isinf(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline double ceil(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            return ceil(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline double floor(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            return floor(x._value());
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline bool isfinite(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            return isfinite(x._value());
        }
        using std::ceil;
        using std::floor;
        using std::isfinite;
        using std::isnan;
        using std::isinf;
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool isnan(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            return isnan(x._value());
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool isinf(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            return isinf(x._value());
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline double ceil(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            return ceil(x._value());
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline double floor(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            return floor(x._value());
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline bool isfinite(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            return isfinite(x._value());
        }
        using std::ceil;
        using std::floor;
        using std::isfinite;
        using std::isnan;
        using std::isinf;
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool isnan(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            return isnan(x._value());
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool isinf(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            return isinf(x._value());
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline double ceil(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            return ceil(x._value());
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline double floor(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            return floor(x._value());
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline bool isfinite(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            return isfinite(x._value());
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline std::ostream& operator << (std::ostream& out, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x) {
            out << x._value();
            return out;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline std::ostream& operator << (std::ostream& out, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x) {
            out << x._value();
            return out;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline std::ostream& operator << (std::ostream& out, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x) {
            out << x._value();
            return out;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline std::ostream& operator << (std::ostream& out, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x) {
            out << x._value();
            return out;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline std::ostream& operator << (std::ostream& out, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x) {
            out << x._value();
            return out;
        }
    }
}
namespace ad {
    namespace internal {
        template<class T> struct is_not_ad_type {
            const static bool RET = true;
        };
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct is_not_ad_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
            const static bool RET=false;
        };
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct is_not_ad_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
            const static bool RET=false;
        };
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct is_not_ad_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
            const static bool RET=false;
        };
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct is_not_ad_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
            const static bool RET=false;
        };
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct is_not_ad_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
            const static bool RET=false;
        };
        template <bool condition, class Then, class Else>
        struct IF {
            typedef Then RET;
        };
        template <class Then, class Else>
        struct IF<false, Then, Else> {
            typedef Else RET;
        };
        template <class T> struct
            active_type_of {
            typedef T RET;
        };
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct active_type_of<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
            typedef ad::internal::active_type< typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE, typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::DATA_TYPE> RET;
        };
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct active_type_of<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
            typedef ad::internal::active_type< typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE, typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::DATA_TYPE> RET;
        };
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct active_type_of<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
            typedef ad::internal::active_type< typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::DATA_TYPE> RET;
        };
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct active_type_of<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
            typedef ad::internal::active_type< typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::DATA_TYPE> RET;
        };
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct active_type_of<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
            typedef ad::internal::active_type< typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE, typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::DATA_TYPE> RET;
        };
        template<class AD_ARG1, class AD_ARG2>
        struct min_max_return_type {
            typedef typename IF<is_not_ad_type<AD_ARG1>::RET, AD_ARG1, typename active_type_of<AD_ARG1>::RET >::RET tmp_type;
            typedef typename IF<is_not_ad_type<AD_ARG2>::RET, AD_ARG2, typename active_type_of<AD_ARG2>::RET >::RET tmp_type2;
            typedef typename IF<is_not_ad_type<AD_ARG1>::RET, tmp_type2, tmp_type>::RET ret_type1;
            typedef typename IF<is_not_ad_type<ret_type1>::RET, AD_ARG1, ret_type1>::RET type;
        };
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class DATA_HANDLER_2 > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> >::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_2> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> >::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A2_T, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A2_T1, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T1, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> >::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A2_T1, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP, class A2_T2, class A2_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> >::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A2_T2, A2_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, AD_TAPE_REAL>::type max(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const AD_TAPE_REAL &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>, AD_TAPE_REAL>::type min(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &a, const AD_TAPE_REAL &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, AD_TAPE_REAL>::type max(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline typename min_max_return_type<ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>, AD_TAPE_REAL>::type min(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, AD_TAPE_REAL>::type max(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>, AD_TAPE_REAL>::type min(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, AD_TAPE_REAL>::type max(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>, AD_TAPE_REAL>::type min(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, AD_TAPE_REAL>::type max(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline typename min_max_return_type<ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>, AD_TAPE_REAL>::type min(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &a, const AD_TAPE_REAL &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >::type max(const AD_TAPE_REAL &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class DATA_HANDLER_1 > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >::type min(const AD_TAPE_REAL &a, const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> >::type max(const AD_TAPE_REAL &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> >::type min(const AD_TAPE_REAL &a, const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> >::type max(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> >::type min(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> >::type max(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T1, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> >::type min(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &b) {
            if (a < b) return a;
            else return b;
        }
        template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> >::type max(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &b) {
            if (a > b) return a;
            else return b;
        } template<class AD_TAPE_REAL, class A1_T2, class A1_OP > static inline typename min_max_return_type<AD_TAPE_REAL, ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> >::type min(const AD_TAPE_REAL &a, const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &b) {
            if (a < b) return a;
            else return b;
        }
    }


}


namespace ad {
    namespace internal {
        template<class AD_TAPE_REAL>
        struct ts_data {
            typedef AD_TAPE_REAL DERIVATIVE_T;
            AD_TAPE_REAL tlm;
            ts_data() : tlm(0) {}
            ts_data &operator = (const ts_data &b) {
                tlm = b.tlm;
                return *this;
            }
            inline void clear() {
                tlm = 0;
            }
            inline const AD_TAPE_REAL &_derivative() const
            {
                (void) "adgt1v";
                return tlm;
            }
            inline AD_TAPE_REAL &_derivative()
            {
                (void) "adgt1v";
                return tlm;
            }
            template<class AD_INTERMEDIATE, class AD_ACTIVE_TYPE>
            static inline void handle(const AD_INTERMEDIATE &vneu, AD_ACTIVE_TYPE &target) {
                ts_data &data = const_cast<ts_data &>(target._data());
                data.tlm = get_tlm(vneu, 1.0);
            }
            template<class DATA_HANDLER_1 >
            static inline AD_TAPE_REAL get_tlm(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x, const AD_TAPE_REAL &pval) {
                return x._data().tlm * pval;
            }
            template<class A1_T1, class A1_T2, class A1_OP >
            static inline AD_TAPE_REAL get_tlm(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x, const AD_TAPE_REAL &pval) {
                return get_tlm(x._arg1, x.pval1() * pval) + get_tlm(x._arg2, x.pval2() * pval);
            }
            template<class A1_T, class A1_OP >
            static inline AD_TAPE_REAL get_tlm(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x, const AD_TAPE_REAL &pval) {
                return get_tlm(x._arg, x.pval() * pval);
            }
            template<class A1_T1, class A1_OP >
            static inline AD_TAPE_REAL get_tlm(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x, const AD_TAPE_REAL &pval) {
                return get_tlm(x._arg1, x.pval1() * pval);
            }
            template<class A1_T2, class A1_OP >
            static inline AD_TAPE_REAL get_tlm(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x, const AD_TAPE_REAL &pval) {
                return get_tlm(x._arg2, x.pval2() * pval);
            }
            typedef void DATA_TAPE_TYPE;
            typedef AD_TAPE_INT TAPE_INDEX_TYPE;
            inline void *_tape() const {
                return NULL;
            }
            inline AD_TAPE_INT _tape_index() const {
                return 0;
            }
        };
    }
}
namespace ad {
    namespace internal {
        class tape_options {
            private:
                AD_TAPE_INT _tapesize;
            public:
                tape_options() : _tapesize(AD_DEFAULT_TAPE_SIZE) { }
                inline AD_TAPE_INT &tapesize() {
                    return _tapesize;
                }
                inline const AD_TAPE_INT &tapesize() const {
                    return _tapesize;
                }
        };
    }
}
namespace ad {
    class exception {
        public:
            template <typename std_exception> static std_exception create(std::string error, std::string file = "", int line = 0) {
                std::stringstream S;
                S << "--- ad/c++ --- " << error;
                if (file != "") S << " --- " << file << ":" << line << ".";
                std::cout << "EXCEPTION thrown: " << S.str() << std::endl;
                return std_exception(S.str());
            }
    };
}
namespace ad {
    namespace internal {
        template <typename AD_TAPE_REAL, class AD_ADJOINT_REAL> struct blob_tape;
    }
/*
    namespace helper {
        template<typename AD_TAPE>
        class callback_object_base {
                template<typename AD_TAPE_REAL, typename AD_ADJOINT_REAL> friend struct ad::internal::blob_tape;
            protected:
                inline virtual ~callback_object_base() {
                }
                AD_TAPE *registered_tape;
            public:
                inline void set_tape(AD_TAPE *t) {
                    if (registered_tape != 0) return;
                    else registered_tape = t;
                }
                inline AD_TAPE *get_tape() {
                    return registered_tape;
                }
                inline callback_object_base() : registered_tape(0) { }
                inline virtual double get_memory_size() {
                    return sizeof(AD_TAPE);
                }
        };
    }
}
*/


  namespace helper {
    template<typename AD_TAPE>
    class callback_object_base {
      template<typename AD_TAPE_REAL, typename AD_ADJOINT_REAL> friend struct ad::internal::blob_tape;
    protected:
      inline virtual ~callback_object_base() {
      }
      AD_TAPE *registered_tape;
    public:
      inline void set_tape(AD_TAPE *t) {
        if (registered_tape != 0) return;
        else registered_tape = t;
      }
      inline AD_TAPE *get_tape() {
        return registered_tape;
      }
      inline callback_object_base() : registered_tape(0) { }
      inline virtual double get_memory_size() {
        return sizeof(AD_TAPE);
      }
    };
    template<class AD_TYPE, class AD_TAPE>
    class userdata_object_base : public callback_object_base<AD_TAPE> {
    private:
      class template_base_class {
      public:
        virtual ~template_base_class() {};
        virtual double size() = 0;
      };
      template <typename X>
      class templated_base_class : public template_base_class {
      public:
        const X _data;
        templated_base_class(const X &d) : _data(d) {}
        const X &get_data() const {
          return _data;
        }
      };
      template <typename X>
      class template_class : public templated_base_class<X> {
      public:
        template_class(X data) : templated_base_class<X>(data) {}
        virtual ~template_class() { }
        virtual double size() {
          return sizeof(X);
        }
      };
      template <typename X>
      class template_vector_class : public templated_base_class<X *> {
      public:
        const int _n;
        template_vector_class(const X *data, int n) : templated_base_class<X *>(new X[n]), _n(n) {
          for (int i = 0; i < n; ++i)
            this->_data[i] = data[i];
        }
        template_vector_class(const X *data, const int inc, const int n) : templated_base_class<X *>(new X[n]), _n(n) {
          for (int i = 0, idx = 0; i < n; ++i, idx += inc)
            this->_data[i] = data[idx];
        }
        virtual ~template_vector_class() {
          delete [] this->_data;
        }
        virtual double size() {
          return _n * sizeof(X);
        }
      };
      unsigned int cp_count;
      std::vector<template_base_class *> checkpoint;
    protected:
      virtual ~userdata_object_base() {
        for (unsigned int i = 0; i < checkpoint.size(); i++)
          delete checkpoint[i];
        checkpoint.clear();
      }
    public:
      userdata_object_base(): callback_object_base<AD_TAPE>(), cp_count(0) {}
      inline virtual double get_memory_size() {
        double S = callback_object_base<AD_TAPE>::get_memory_size();
        for (unsigned int i = 0; i < checkpoint.size(); i++)
          S += checkpoint[i]->size();
        return S;
      }
      template<typename X>
      inline void write_data(const X &cp) {
        checkpoint.push_back(new template_class<X>(cp));
      }
      template<typename X>
      inline void write_data(const X *const cp, const int n) {
        checkpoint.push_back(new template_vector_class<X>(cp, n));
      }
      template<typename X>
      inline void write_data(const X *const &cp, const int inc, const int n) {
        checkpoint.push_back(new template_vector_class<X>(cp, inc, n));
      }
      template<typename X>
      inline const X &read_data() {
        const X &cp = static_cast<templated_base_class<X>* >(checkpoint[cp_count])->get_data();
        ++cp_count;
        if (static_cast<size_t>(cp_count) == checkpoint.size()) cp_count = 0;
        return cp;
      }
    };
    template<class AD_TYPE, class AD_TAPE>
    class external_adjoint_object_base : public userdata_object_base<AD_TYPE, AD_TAPE> {
    protected:
      std::vector<AD_TAPE_INT> inputs;
      std::vector<AD_TAPE_INT> outputs;
      AD_TAPE_INT inputs_count;
      AD_TAPE_INT outputs_count;
    public:
      inline size_t get_number_of_registered_inputs() {
        return inputs.size();
      }
      inline size_t get_number_of_registered_outputs() {
        return outputs.size();
      }
    public:
      inline void check_tape(const AD_TYPE &x) {
        if ((x._data()._tape() != 0) && (this->registered_tape != x._data()._tape()))
          throw ad::exception::create<std::runtime_error>("impossible binding tape - wrong tape in variable!");
      }
    protected:
      ~external_adjoint_object_base() { }
    public:
      external_adjoint_object_base(const std::pair<int, int> &a): userdata_object_base<AD_TYPE, AD_TAPE>(),
        inputs_count(0),
        outputs_count(0) {
        inputs.reserve(a.first);
        outputs.reserve(a.second);
      }
      external_adjoint_object_base(): userdata_object_base<AD_TYPE, AD_TAPE>(), inputs_count(0), outputs_count(0) {
      }
      inline typename AD_TYPE::VALUE_TYPE register_input(const AD_TYPE &x) {
        check_tape(x);
        inputs.push_back(x._data()._tape_index());
        return x._value();
      }
      inline void register_input(const AD_TYPE *const x, typename AD_TYPE::VALUE_TYPE *values, const int n) {
        AD_TAPE_INT oldsize = inputs.size();
        inputs.resize(oldsize + n);
        AD_TAPE_INT *pos = &inputs[oldsize];
        int i;
        for (i = 0; i < n; ++i) {
          check_tape(x[i]);
          pos[i] = x[i]._data()._tape_index();
          values[i] = x[i]._value();
        }
      }
      inline void register_input(const std::vector<AD_TYPE> &x, std::vector<typename AD_TYPE::VALUE_TYPE> &values) {
        assert(x.size() == values.size());
        register_input(&(x[0]), &(values[0]), x.size());
      }
      inline std::vector<typename AD_TYPE::VALUE_TYPE> register_input(const std::vector<AD_TYPE> &x) {
        std::vector<typename AD_TYPE::VALUE_TYPE> values(x.size());
        register_input(x, values);
        return values;
      }
      inline void register_output(AD_TYPE *actives, const size_t n) {
        if (this->registered_tape == NULL) {
          throw ad::exception::create<std::runtime_error>("impossible binding output - no tape available");
        } else {
          AD_TAPE_INT oldsize = outputs.size();
          outputs.resize(oldsize + n);
          AD_TAPE_INT *outs = &outputs[ oldsize ];
          size_t startindex = 0;
          typename AD_TAPE::TAPE_ENTRY *ins = this->registered_tape->_get_insert_ptr_range(n, startindex);
          for (size_t i = 0; i < n; ++i) {
            ins[i].arg = 0;
            typename AD_TYPE::DATA_TYPE &data = actives[i]._data();
            data.register_variable(startindex + i, this->registered_tape);
            outs[i] = static_cast<AD_TAPE_INT>(startindex + i);
          }
        }
      }
      inline void register_output(const typename AD_TYPE::VALUE_TYPE *const pvalues, AD_TYPE *actives, const size_t n) {
        if (this->registered_tape == NULL) {
          throw ad::exception::create<std::runtime_error>("impossible binding output - no tape available");
        } else {
          AD_TAPE_INT oldsize = outputs.size();
          outputs.resize(oldsize + n);
          AD_TAPE_INT *outs = &outputs[ oldsize ];
          size_t startindex = 0;
          typename AD_TAPE::TAPE_ENTRY *ins = this->registered_tape->_get_insert_ptr_range(n, startindex);
          for (size_t i = 0; i < n; ++i) {
            ins[i].arg = 0;
            actives[i] = pvalues[i];
            typename AD_TYPE::DATA_TYPE &data = const_cast<typename AD_TYPE::DATA_TYPE &>(actives[i]._data());
            data.register_variable(startindex + i, this->registered_tape);
            outs[i] = static_cast<AD_TAPE_INT>(startindex + i);
          }
        }
      }
      inline void register_output(const std::vector<typename AD_TYPE::VALUE_TYPE> &pvalues, std::vector<AD_TYPE> &actives) {
        assert(pvalues.size() == actives.size());
        register_output(&(pvalues[0]), &(actives[0]), pvalues.size());
      }
      inline std::vector<AD_TYPE> register_output(const std::vector<typename AD_TYPE::VALUE_TYPE> &pvalues) {
        std::vector<AD_TYPE> actives(pvalues.size());
        register_output(pvalues, actives);
        return actives;
      }
      inline void register_output(std::vector<AD_TYPE> &actives) {
        register_output(&(actives[0]), actives.size());
      }
      inline AD_TYPE register_output(const typename AD_TYPE::VALUE_TYPE &py, AD_TAPE *tape = NULL) {
        AD_TYPE y;
        if (tape != NULL) {
          if (this->registered_tape != NULL && this->registered_tape != tape) {
            throw ad::exception::create<std::runtime_error>("impossible binding output in external function (register_output) - tape of inputs and outputs differ!");
          }
          this->registered_tape = tape;
        }
        if (this->registered_tape != NULL) {
          y = py;
          this->registered_tape->register_variable(y);
        } else
          throw ad::exception::create<std::runtime_error>("impossible binding output in external function - no tape available");
        outputs.push_back(y._data()._tape_index());
        return y;
      }
      inline typename AD_TYPE::VALUE_TYPE get_output_adjoint() {
        AD_TAPE_INT idx = outputs_count;
        outputs_count++;
        if (static_cast<size_t>(outputs_count) == outputs.size())
          outputs_count = 0;
        typename AD_TYPE::VALUE_TYPE back = 0;
        back = this->registered_tape->_adjointEx(outputs[static_cast<size_t>(idx)]);
        return back;
      }
      inline void get_output_adjoint(typename AD_TYPE::VALUE_TYPE *buffer, const size_t n) {
        AD_TAPE_INT idx = outputs_count;
        for (size_t i = 0; i < n; ++i) {
          buffer[i] = this->registered_tape->_adjoint(outputs[idx]);
          idx++;
        }
        outputs_count += n;
        if (static_cast<size_t>(outputs_count) == outputs.size())
          outputs_count = 0;
      }
      inline void get_output_adjoint(std::vector<typename AD_TYPE::VALUE_TYPE> &buffer) {
        assert(buffer.size());
        get_output_adjoint(&(buffer[0]), buffer.size());
      }
      inline void increment_input_adjoint(const typename AD_TYPE::VALUE_TYPE *const adj, const int n) {
        for (int i = 0; i < n; ++i) {
          this->registered_tape->_adjoint(inputs[inputs_count + i]) += adj[i];
        }
        inputs_count += n;
        if (static_cast<size_t>(inputs_count) == inputs.size())
          inputs_count = 0;
      }
      inline void increment_input_adjoint(const std::vector<typename AD_TYPE::VALUE_TYPE> &adj) {
        assert(adj.size() != 0);
        increment_input_adjoint(&(adj[0]), adj.size());
      }
      inline bool all_adjoints_written() {
        if (inputs_count == 0) return true;
        else return false;
      }
      inline bool all_adjoints_read() {
        if (outputs_count == 0) return true;
        else return false;
      }
      inline void increment_input_adjoint(const typename AD_TYPE::VALUE_TYPE &adj) {
        AD_TAPE_INT idx = inputs_count;
        inputs_count++;
        if (static_cast<size_t>(inputs_count) == inputs.size())
          inputs_count = 0;
        this->registered_tape->_adjointEx(inputs[static_cast<size_t>(idx)]) += adj;
      }
    };
  }
}

namespace ad {
    template< typename T > struct trait_value {
        typedef T RETURN_TYPE;
        static inline RETURN_TYPE &value(T &value) {
            return value;
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_value <ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_value <const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_value <ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_value <const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_value <ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_value <const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_value <ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_value <const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_value <ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &value) {
            return value._value();
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_value <const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &value) {
            return value._value();
        }
    };
    template <typename T>
    inline typename trait_value<T>::RETURN_TYPE &value(T &x) {
        return trait_value<T>::value(x);
    }
    template <typename T>
    inline typename trait_value<const T>::RETURN_TYPE &value(const T &x) {
        return trait_value<const T>::value(x);
    }
    template< typename T > struct trait_passive_value : trait_value<T> {};
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_passive_value <ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename trait_passive_value<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
            return trait_passive_value<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_passive_value <const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef const typename trait_passive_value<typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
            return trait_passive_value<const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_passive_value <ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef typename trait_passive_value<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
            return trait_passive_value<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_passive_value <const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef const typename trait_passive_value<typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
            return trait_passive_value<const typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_passive_value <ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef typename trait_passive_value<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
            return trait_passive_value<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_passive_value <const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef const typename trait_passive_value<typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
            return trait_passive_value<const typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_passive_value <ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef typename trait_passive_value<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
            return trait_passive_value<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_passive_value <const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef const typename trait_passive_value<typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
            return trait_passive_value<const typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_passive_value <ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef typename trait_passive_value<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
            return trait_passive_value<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_passive_value <const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef const typename trait_passive_value<typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE>::RETURN_TYPE RETURN_TYPE;
        static inline RETURN_TYPE& value(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
            return trait_passive_value<const typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::VALUE_TYPE>::value(ad::value(x));
        }
    };
    template <typename T>
    inline typename trait_passive_value<T>::RETURN_TYPE &passive_value(T &x) {
        return trait_passive_value<T>::value(x);
    }
    template <typename T>
    inline const typename trait_passive_value<T>::RETURN_TYPE &passive_value(const T &x) {
        return trait_passive_value<const T>::value(x);
    }
    template< typename T > struct trait_derivative {
        typedef T RETURN_TYPE;
        static inline RETURN_TYPE value(const T &vvalue) {
            (void) vvalue;
            return RETURN_TYPE();
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_derivative <ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::DATA_TYPE::DERIVATIVE_T& RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return value._data()._derivative();
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_derivative <const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::DATA_TYPE::DERIVATIVE_T& RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return const_cast<RETURN_TYPE>(value._data()._derivative());
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_derivative <ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_derivative <const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_derivative <ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_derivative <const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_derivative <ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_derivative <const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_derivative <ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> >;
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_derivative <const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> >;
    template <typename T>
    inline typename trait_derivative<T>::RETURN_TYPE derivative(T &x) {
        return trait_derivative<T>::value(x);
    }
    template <typename T>
    inline typename trait_derivative<const T>::RETURN_TYPE derivative(const T &x) {
        return trait_derivative<const T>::value(x);
    }
    template< typename T > struct trait_tape_index {
        typedef AD_TAPE_INT RETURN_TYPE;
        static inline RETURN_TYPE value(const T &vvalue) {
            (void) vvalue;
            return 0;
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_tape_index <ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class DATA_HANDLER_1 > struct trait_tape_index <const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > {
        typedef typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_tape_index <ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T, class A1_OP > struct trait_tape_index <const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> > {
        typedef typename ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_tape_index <ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP > struct trait_tape_index <const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_tape_index <ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T1, class A1_OP > struct trait_tape_index <const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_tape_index <ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template<class AD_TAPE_REAL, class A1_T2, class A1_OP > struct trait_tape_index <const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> > {
        typedef typename ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>::DATA_TYPE::TAPE_INDEX_TYPE RETURN_TYPE;
        static inline RETURN_TYPE value(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &value) {
            return value._data()._tape_index();
        }
    };
    template <typename T>
    inline typename trait_tape_index<T>::RETURN_TYPE tape_index(T &x) {
        return trait_tape_index<T>::value(x);
    }
    template <typename T>
    inline typename trait_tape_index<const T>::RETURN_TYPE tape_index(const T &x) {
        return trait_tape_index<const T>::value(x);
    }
}
namespace ad {
    namespace internal {
        template<class AD_TAPE_REAL, class AD_ADJOINT_REAL = AD_TAPE_REAL>
        struct blob_tape {
            public:
                struct TAPE_ENTRY {
                    AD_TAPE_INT arg;
                    AD_TAPE_REAL pval;
                };
                struct position_t {
                        friend struct blob_tape;
                    private:
                        AD_TAPE_INT stackcounter;
                        AD_TAPE_INT progvarcounter;
                        position_t(const AD_TAPE_INT nstackcounter, const AD_TAPE_INT nprogvarcounter) : stackcounter(nstackcounter), progvarcounter(nprogvarcounter) {}
                    public:
                        position_t() : stackcounter(0), progvarcounter(0) {}
                        inline const AD_TAPE_INT &_stackcounter() const {
                            return stackcounter;
                        }
                        inline const AD_TAPE_INT &_progvarcounter() const {
                            return progvarcounter;
                        }
                        inline bool operator== (const position_t &other) const {
                            return _stackcounter() == other._stackcounter();
                        }
                        inline bool operator> (const position_t &other) const {
                            return _stackcounter() > other._stackcounter();
                        }
                        inline bool operator< (const position_t &other) const {
                            return _stackcounter() < other._stackcounter();
                        }
                };
            private:
                TAPE_ENTRY *_topOfStack;
                AD_TAPE_INT _progvarcounter;
                AD_ADJOINT_REAL *_adjoints;
                TAPE_ENTRY *_stack;
                AD_TAPE_INT _stack_size;
                AD_TAPE_INT _adjoint_size;
                bool _isdead;
                bool _isactive;
                int _vecidx;
            private:
                ~blob_tape() {
                    reset();
                    delete [] _stack;
                    delete [] _adjoints;
                }
            public:
                blob_tape() : _topOfStack(0), _progvarcounter(0), _adjoints(0), _stack(0), _stack_size(0), _adjoint_size(0), _isdead(false), _isactive(false), _vecidx(0) {
                }
            public:
                inline void _finish_current_insert_ptr(TAPE_ENTRY *const end) {
                    _topOfStack = end + 1;
                }
                inline TAPE_ENTRY *_get_insert_ptr(const int num_entries2fill, AD_TAPE_INT &new_tape_index) {
                    ;
                    _progvarcounter++;
                    new_tape_index = _progvarcounter;
                    TAPE_ENTRY *ret = _topOfStack;
                    _topOfStack += num_entries2fill;
                    return ret;
                }
                inline TAPE_ENTRY *_get_insert_ptr_range(const int num_entries2fill, AD_TAPE_INT &new_tape_index) {
                    new_tape_index = _progvarcounter + 1;
                    ;
                    _progvarcounter += num_entries2fill;
                    TAPE_ENTRY *ret = _topOfStack;
                    _topOfStack += num_entries2fill;
                    return ret;
                }
            public:
                struct interpretation_settings;
            private:
                void _interpret_chunk(TAPE_ENTRY *start, TAPE_ENTRY *end, AD_TAPE_INT &progvaridx, const interpretation_settings &settings) {
                    TAPE_ENTRY *cur = start;
                    if (settings.zeroadjoints) {
                        while (cur != end) {
                            cur--;
                            AD_ADJOINT_REAL &adj = _adjoints[progvaridx];
                            --progvaridx;
                            const int &edgecount = cur->arg;
                            for (int i = 0; i < edgecount; ++i) {
                                cur--;
                                _adjoints[cur->arg] += adj * cur->pval;
                            }
                            adj = 0;
                        }
                    } else {
                        while (cur != end) {
                            cur--;
                            AD_ADJOINT_REAL &adj = _adjoints[progvaridx];
                            --progvaridx;
                            const int &edgecount = cur->arg;
                            for (int i = 0; i < edgecount; ++i) {
                                cur--;
                                _adjoints[cur->arg] += adj * cur->pval;
                            }
                        }
                    }
                }
                inline void _interpret_adjoint_internal_plain(const position_t &from, const position_t &to, const interpretation_settings &settings) {
                    AD_TAPE_INT progvaridx = from._progvarcounter();
                    TAPE_ENTRY *start = _stack + from._stackcounter();
                    TAPE_ENTRY *end = _stack + to._stackcounter();
                    _interpret_chunk(start, end, progvaridx, settings);
                }
                inline void _zero_adjoints_internal(const position_t &from, const position_t &to) {
                    if (_adjoints == 0) return;
                    for (AD_TAPE_INT i = from._progvarcounter(); i > to._progvarcounter(); --i) {
                        _adjoints[i] = 0;
                    }
                }
                inline void _reset_to_internal(const position_t &to) {
                    _zero_adjoints_internal(get_position(), to);
                    _topOfStack = _stack + to._stackcounter();
                    _progvarcounter = to._progvarcounter();
                }
            public:
                template<class AD_ACTIVE>
                inline void _register_variables_internal(AD_ACTIVE *actives, const typename AD_ACTIVE::VALUE_TYPE *values, int *outs, const int n) {
                    TAPE_ENTRY *entries = _topOfStack;
                    ;
                    AD_TAPE_INT startidx = _progvarcounter + 1;
                    AD_TAPE_INT i;
                    typedef typename AD_ACTIVE::DATA_TYPE AD_DATA;
                    for (i = 0; i < n; ++i) {
                        entries[i].arg = 0;
                        if (values) actives[i] = values[i];
                        if (outs) outs[i] = startidx + i;
                        AD_DATA &data = const_cast<AD_DATA &>(actives[i]._data());
                        data.register_variable(startidx + i, this);
                    }
                    _topOfStack += n;
                    ;
                    _progvarcounter += n;
                }
                inline position_t get_position() const {
                    return position_t(static_cast<AD_TAPE_INT>(_topOfStack - _stack), _progvarcounter);
                }
                inline AD_TAPE_REAL &_adjointEx(const size_t tape_index) {
                    AD_ADJOINT_REAL &adj = _adjoint(tape_index);
                    AD_TAPE_REAL *adjvec;
                    adjvec = reinterpret_cast<AD_TAPE_REAL *>(&adj);
                    return adjvec[_vecidx];
                }
                inline AD_ADJOINT_REAL &_adjoint(const size_t tape_index) {
                    return _adjoints[tape_index];
                }
                inline double get_checkpoint_memory_size() {
                    double checkpoint_size = 0;
                    for (AD_TAPE_INT i = 0; i < tape_callbacks.size(); i++)
                        checkpoint_size += tape_callbacks[i].userdata->get_memory_size();
                    return checkpoint_size / 1024.0 / 1024.0;
                }
                inline size_t get_tape_memory_size() const {
                    return (get_position()._progvarcounter() * sizeof(AD_TAPE_REAL) + get_position()._stackcounter() * sizeof(TAPE_ENTRY));
                }
                inline double get_allocated_tape_memory_size() const {
                    return ((_adjoint_size / 1024.0 / 1024.0) * sizeof(AD_TAPE_REAL) + (_stack_size / 1024.0 / 1024.0) * sizeof(TAPE_ENTRY));
                }
                inline size_t _get_tape_memory() const {
                    return (((get_position()._progvarcounter()) * sizeof(AD_TAPE_REAL)) + ((get_position()._stackcounter()) * sizeof(TAPE_ENTRY)));
                }
                static blob_tape *create(tape_options options = tape_options()) {
                    return create(options.tapesize());
                }
                static blob_tape *create(AD_TAPE_INT size, AD_TAPE_INT progvarcounter = 0) {
                    if (progvarcounter == 0) progvarcounter = size / 2;
                    blob_tape *ret = new blob_tape();
                    ret->_stack = new TAPE_ENTRY[size + 1];
                    ret->_topOfStack = ret->_stack;
                    ret->_adjoints = new AD_ADJOINT_REAL[progvarcounter + 1];
                    for (AD_TAPE_INT i = 0; i <= progvarcounter; i++)
                        ret->_adjoints[i] = 0;
                    ret->_stack_size = size;
                    ret->_adjoint_size = progvarcounter;
                    ret->_isactive = true;
                    for (AD_TAPE_INT i = 0; i <= size; ++i) {
                        ret->_stack[i].arg = 0;
                        ret->_stack[i].pval = 0;
                    }
                    return ret;
                }
                static void remove(blob_tape *&tape) {
                    if (tape == 0) return;
                    delete tape;
                    tape = 0;
                }
                typedef ad::internal::blob_tape<AD_TAPE_REAL, AD_ADJOINT_REAL> AD_TAPE_CLASS;
            public:
                struct interpretation_settings {
                    bool reset;
                    bool zeroadjoints;
                    interpretation_settings() : reset(false), zeroadjoints(false) {}
                };
                typedef ad::helper::callback_object_base<AD_TAPE_CLASS> callback_object_t;

                template <typename EXT_DATA>
                class CALLBACK_DATA_POINTER {
                    public:
                        typedef void (*TAPE_CALLBACK_w_all_base)(AD_TAPE_CLASS &caller, const interpretation_settings &s, EXT_DATA *userdata);
                        typedef void (*TAPE_CALLBACK_w_tape_base)(AD_TAPE_CLASS &caller, EXT_DATA *userdata);
                        typedef void (*TAPE_CALLBACK_plain_base)(EXT_DATA *userdata);
                };
                class CALLBACK_FCN_HANDLER_BASE {
                    public:
                        virtual void run_callback(AD_TAPE_CLASS &caller, const interpretation_settings &s, callback_object_t *userdata) = 0;
                        virtual ~CALLBACK_FCN_HANDLER_BASE() {};
                };
                template <typename EXT_DATA>
                class CALLBACK_FCN_HANDLER : public CALLBACK_FCN_HANDLER_BASE {
                    private:
                        union {
                            typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_plain_base fcn;
                            typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_w_tape_base fcn_w_tape;
                            typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_w_all_base fcn_w_all;
                        } fcn;
                        int fcn_type_id;
                    public:
                        CALLBACK_FCN_HANDLER(typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_plain_base fcn_) : fcn_type_id(0) {
                            fcn.fcn = fcn_;
                        }
                        CALLBACK_FCN_HANDLER(typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_w_tape_base fcn_) : fcn_type_id(1) {
                            fcn.fcn_w_tape = fcn_;
                        }
                        CALLBACK_FCN_HANDLER(typename CALLBACK_DATA_POINTER<EXT_DATA>::TAPE_CALLBACK_w_all_base fcn_) : fcn_type_id(2) {
                            fcn.fcn_w_all = fcn_;
                        }
                        void run_callback(AD_TAPE_CLASS &caller, const interpretation_settings &s, callback_object_t *userdata) {
                            EXT_DATA *casted_userdata = static_cast<EXT_DATA *>(userdata);
                            switch (fcn_type_id) {
                            case 0:
                                fcn.fcn(casted_userdata);
                                break;
                            case 1:
                                fcn.fcn_w_tape(caller, casted_userdata);
                                break;
                            case 2:
                                fcn.fcn_w_all(caller, s, casted_userdata);
                                break;
                            default:
                                throw ad::exception::create<std::runtime_error>("unknown callback to run.");
                                break;
                            }
                        }
                        ~CALLBACK_FCN_HANDLER() {}
                };
                class tape_callback {
                        callback_object_t *userdata;
                        CALLBACK_FCN_HANDLER_BASE *callback_handler;
                        position_t position;
                    public:
                        tape_callback() : userdata(0),
                            callback_handler(0) { }
                        template <typename EXT_DATA, typename FCN_CALLBACK>
                        void set_callback(FCN_CALLBACK fcn_) {
                            if (callback_handler)
                                throw ad::exception::create<std::runtime_error>("currently not supported to insert external_adjoint_object_bases twice.");
                            callback_handler = new CALLBACK_FCN_HANDLER<EXT_DATA>(fcn_);
                        }
                        callback_object_t *&_userdata() {
                            return userdata;
                        }
                        void free_userdata() {
                            delete userdata;
                            if (callback_handler)
                                delete callback_handler;
                        }
                        position_t &_position() {
                            return position;
                        }
                        void run_callback(AD_TAPE_CLASS &caller, const interpretation_settings &s) {
                            if (callback_handler)
                                callback_handler->run_callback(caller, s, userdata);
                        }
                };
                std::vector<tape_callback> tape_callbacks;
                template <class ext_fcn_data_type, typename FCN_PARAMETERS>
                inline ext_fcn_data_type *create_callback_object(const FCN_PARAMETERS &parameters) {
                    const reference_constructor_wrapper<FCN_PARAMETERS> ref_wrapper(parameters);
                    tape_callback external_fcn;
                    tape_callbacks.push_back(external_fcn);
                    tape_callbacks.back()._position() = get_position();
                    ext_fcn_data_type *userdata = ref_wrapper.template create<ext_fcn_data_type>();
                    tape_callbacks.back()._userdata() = userdata;
                    userdata->set_tape(this);
                    return userdata;
                }
                template <class ext_fcn_data_type>
                inline ext_fcn_data_type *create_callback_object() {
                    void *dummy;
                    return create_callback_object<ext_fcn_data_type>(dummy);
                }
                template <class ext_fcn_data_type, typename FCN_CALLBACK>
                inline void insert_callback(FCN_CALLBACK callback_handler, ext_fcn_data_type *D) {
                    if (tape_callbacks.back()._userdata() == D) {
                        tape_callbacks.back().template set_callback<ext_fcn_data_type>(callback_handler);
                        tape_callbacks.back()._position() = get_position();
                        AD_TAPE_INT tmp = 0;
                        _get_insert_ptr(1, tmp)->arg = 0;
                    } else {
                        throw ad::exception::create<std::runtime_error>("please always insert most recently created external function object.");
                    }
                }
            private:
                inline void _reset_tape_callbacks_to(const position_t &to) {
                    for (int i = static_cast<int>(tape_callbacks.size()) - 1; i >= 0; i--) {
                        size_t offset = static_cast<size_t>(i);
                        if (tape_callbacks[offset]._position()._progvarcounter() > to._progvarcounter()) {
                            tape_callbacks[offset].free_userdata();
                            tape_callbacks.pop_back();
                        } else {
                            break;
                        }
                    }
                }
            private:
                inline void _interpret_adjoint_internal(const position_t &from, const position_t &to, const interpretation_settings &settings)
                {
                    if (from > get_position())
                        throw ad::exception::create<std::runtime_error>("you try to use a tape position outside of the current tape. error.");
                    _adjoint(static_cast<AD_TAPE_INT>(from._progvarcounter()));
                    int external_first = -1;
                    int external_count = 0;
                    for (int i = static_cast<int>(tape_callbacks.size()) - 1; i >= 0; --i) {
                        size_t offset = static_cast<size_t>(i);
                        if ((tape_callbacks[offset]._position()._progvarcounter() <= from._progvarcounter()) &&
                                (tape_callbacks[offset]._position()._progvarcounter() >= to._progvarcounter())) {
                            if (external_first == -1) external_first = i;
                            ++external_count;
                        }
                    }
                    position_t myfrom = from;
                    position_t myto;
                    for (int i = external_first; external_count > 0; --external_count) {
                        size_t offset = static_cast<size_t>(i);
                        myto = tape_callbacks[offset]._position();
                        _interpret_adjoint_internal_plain(myfrom, myto, settings);
                        if (settings.reset) {
                            _reset_to_internal(myto);
                            AD_TAPE_INT nti;
                            TAPE_ENTRY *cur = _get_insert_ptr(1, nti);
                            (void) cur;
                        }
                        tape_callbacks[offset].run_callback(*this, settings);
                        if (settings.reset) {
                            _reset_to_internal(myto);
                            tape_callbacks[offset].free_userdata();
                            tape_callbacks.pop_back();
                        }
                        myfrom = myto;
                        --i;
                    }
                    _interpret_adjoint_internal_plain(myfrom, to, settings);
                }
            public:
                inline bool is_active() {
                    return _isactive;
                }
                inline void switch_to_active() {
                    if (!_isactive) _isactive = true;
                }
                inline void switch_to_passive() {
                    if (_isactive) _isactive = false;
                }
                template<class DATA_HANDLER_1 >
                inline void register_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
                    AD_TAPE_INT new_tape_index;
                    TAPE_ENTRY *cur = _get_insert_ptr(1, new_tape_index);
                    (void) cur;
                    cur->arg = 0;
                    DATA_HANDLER_1 &data = const_cast<DATA_HANDLER_1 &>(x._data());
                    data.register_variable(new_tape_index, this);
                }
                template<class DATA_HANDLER_1 >
                inline void register_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE &v) {
                    x = v;
                    register_variable(x);
                }
                template<class DATA_HANDLER_1 >
                inline void register_variable(std::vector<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > &x) {
                    for (size_t i = 0; i < x.size(); i++)
                        register_variable(x[i]);
                }
                template<class DATA_HANDLER_1 >
                inline void register_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> *x, const int n, const typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>::VALUE_TYPE *const v) {
                    if (n == 0) return;
                    for (int i = 0; i < n; ++i)
                        register_variable(x[i], v[i]);
                }
                template<class DATA_HANDLER_1 >
                inline void register_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> *x, const size_t n) {
                    if (n == 0) return;
                    for (size_t i = 0; i < n; ++i)
                        register_variable(x[i]);
                }
                template<class DATA_HANDLER_1 >
                inline void register_output_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
                    x = static_cast<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >(1.0) * x;
                }
                template<class DATA_HANDLER_1 >
                inline void register_output_variable(ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> *x, const size_t n) {
                    if (n == 0) return;
                    for (size_t i = 0; i < n; ++i)
                        x[i] = static_cast<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> >(1.0) * x[i];
                }
                template<class DATA_HANDLER_1 >
                inline void register_output_variable(std::vector<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> > &x) {
                    for (size_t i = 0; i < x.size(); i++)
                        register_output_variable(x[i]);
                }
                inline void reset_to(const position_t &to) {
                    _reset_to_internal(to);
                    _reset_tape_callbacks_to(to);
                }
                inline void reset() {
                    reset_to(position_t());
                }
                inline void interpret_adjoint() {
                    position_t to;
                    interpretation_settings settings;
                    _interpret_adjoint_internal(get_position(), to, settings);
                }
                inline void interpret_adjoint_to(const position_t &to) {
                    interpretation_settings settings;
                    if (to > get_position())
                        throw ad::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "../build_files//../src/ad//ad_tape_interface_inc.hpp", 241);
                    else
                        _interpret_adjoint_internal(get_position(), to, settings);
                }
                inline void interpret_adjoint_from(const position_t &from) {
                    position_t to;
                    assert(!(from < to));
                    interpretation_settings settings;
                    _interpret_adjoint_internal(from, to, settings);
                }
                inline void interpret_adjoint_from_to(const position_t &from, const position_t &to) {
                    interpretation_settings settings;
                    if (to > from)
                        throw ad::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "../build_files//../src/ad//ad_tape_interface_inc.hpp", 264);
                    else
                        _interpret_adjoint_internal(from, to, settings);
                }
                inline void interpret_adjoint_and_reset_to(const position_t &to) {
                    position_t from(get_position());
                    interpretation_settings settings;
                    settings.reset = true;
                    settings.zeroadjoints = true;
                    _interpret_adjoint_internal(from, to, settings);
                    _reset_to_internal(to);
                    _reset_tape_callbacks_to(to);
                }
                inline void interpret_adjoint_and_zero_adjoints_to(const position_t &to) {
                    position_t from(get_position());
                    interpretation_settings settings;
                    settings.reset = false;
                    settings.zeroadjoints = true;
                    _interpret_adjoint_internal(from, to, settings);
                }
                inline void interpret_adjoint_and_zero_adjoints_from_to(const position_t &from, const position_t &to) {
                    interpretation_settings settings;
                    settings.reset = false;
                    settings.zeroadjoints = true;
                    _interpret_adjoint_internal(from, to, settings);
                }
                inline void zero_adjoints() {
                    position_t to;
                    _zero_adjoints_internal(get_position(), to);
                }
                inline void zero_adjoints_to(const position_t &to) {
                    _zero_adjoints_internal(get_position(), to);
                }
                inline void zero_adjoints_from(const position_t &from) {
                    position_t to;
                    _zero_adjoints_internal(from, to);
                }
                inline void zero_adjoints_from_to(const position_t &from, const position_t &to) {
                    _zero_adjoints_internal(from, to);
                }
                struct handler_base {
                    protected:
                        template<class DATA_HANDLER_1 >
                        static inline int _get_edge_count(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x) {
                            return x._data()._edgecount();
                        }
                        template<class A1_T1, class A1_T2, class A1_OP >
                        static inline int _get_edge_count(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x) {
                            return _get_edge_count(x._arg1) + _get_edge_count(x._arg2);
                        }
                        template<class A1_T, class A1_OP >
                        static inline int _get_edge_count(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x) {
                            return _get_edge_count(x._arg);
                        }
                        template<class A1_T1, class A1_OP >
                        static inline int _get_edge_count(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x) {
                            return _get_edge_count(x._arg1);
                        }
                        template<class A1_T2, class A1_OP >
                        static inline int _get_edge_count(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x) {
                            return _get_edge_count(x._arg2);
                        }
                        struct tapehandler {
                            TAPE_ENTRY *_ins_ptr;
                            template<class AD_INTERMEDIATE>
                            tapehandler(TAPE_ENTRY *ins_ptr, const int edgecount, AD_INTERMEDIATE vneu) : _ins_ptr(ins_ptr) {
                                interpret(vneu, 1.0);
                                _ins_ptr->arg = edgecount;
                            }
                            template<class DATA_HANDLER_1 >
                            inline void interpret(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1> &x, const AD_TAPE_REAL &pval) {
                                if (x._data()._tape_index() != 0) {
                                    _ins_ptr->arg = x._data()._tape_index();
                                    _ins_ptr->pval = pval;
                                    ++_ins_ptr;
                                }
                            }
                            template<class A1_T1, class A1_T2, class A1_OP >
                            inline void interpret(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP> &x, const AD_TAPE_REAL &pval) {
                                this->interpret(x._arg1, x.pval1()*pval);
                                this->interpret(x._arg2, x.pval2()*pval);
                            }
                            template<class A1_T, class A1_OP >
                            inline void interpret(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP> &x, const AD_TAPE_REAL &pval) {
                                this->interpret(x._arg, x.pval() * pval);
                            }
                            template<class A1_T1, class A1_OP >
                            inline void interpret(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP> &x, const AD_TAPE_REAL &pval) {
                                this->interpret(x._arg1, x.pval1()*pval);
                            }
                            template<class A1_T2, class A1_OP >
                            inline void interpret(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP> &x, const AD_TAPE_REAL &pval) {
                                this->interpret(x._arg2, x.pval2()*pval);
                            }
                        };
                };
                template<AD_TAPE_CLASS *&global_tape>
                struct single_tape_data : handler_base {
                        typedef AD_TAPE_CLASS DATA_TAPE_TYPE;
                        typedef const AD_TAPE_INT &TAPE_INDEX_TYPE;
                        typedef AD_ADJOINT_REAL DERIVATIVE_T;
                        AD_TAPE_INT _tape_index_;
                    public:
                        single_tape_data() : _tape_index_(0) {}
                        inline void clear() {
                            _tape_index_ = 0;
                        }
                        inline const AD_ADJOINT_REAL &_derivative() const {
                            return _adjoint();
                        }
                        inline AD_ADJOINT_REAL &_derivative() {
                            return _adjoint();
                        }
                        inline AD_ADJOINT_REAL &_adjoint() const {
                            return global_tape->_adjoint(_tape_index_);
                        }
                        inline bool _is_registered() const {
                            return _tape_index_ == 0 ? false : true;
                        }
                        inline AD_TAPE_INT &_tape_index() {
                            return _tape_index_;
                        }
                        inline const AD_TAPE_INT &_tape_index() const {
                            return _tape_index_;
                        }
                        inline int _edgecount() const {
                            return _tape_index_ == 0 ? 0 : 1;
                        }
                        inline void register_variable(AD_TAPE_INT new_tape_index, AD_TAPE_CLASS *tape) {
                            (void) tape;
                            _tape_index_ = new_tape_index;
                        }
                        template<class AD_INTERMEDIATE, class AD_ACTIVE_TYPE>
                        static inline void handle(const AD_INTERMEDIATE &vneu, AD_ACTIVE_TYPE &target) {
                            int edgecount = handler_base::_get_edge_count(vneu);
                            if (edgecount > 0) {
                                if (NULL != global_tape && !global_tape->is_active()) {
                                    single_tape_data &data = const_cast<single_tape_data &>(target._data());
                                    data.clear();
                                    return;
                                }
                                AD_TAPE_INT newTapeIndex;
                                typename AD_TAPE_CLASS::TAPE_ENTRY *ins_ptr = global_tape->_get_insert_ptr(edgecount + 1, newTapeIndex);
                                if (ins_ptr != 0) {
                                    typename handler_base::tapehandler tmp(ins_ptr, edgecount, vneu);
                                }
                                single_tape_data &data = const_cast<single_tape_data &>(target._data());
                                data._tape_index_ = newTapeIndex;
                            } else {
                                single_tape_data &data = const_cast<single_tape_data &>(target._data());
                                data.clear();
                            }
                        }
                        static inline AD_TAPE_CLASS *_tape() {
                            return global_tape;
                        }
                        inline void _set_tape(AD_TAPE_CLASS *T) {
                            (void) T;
                        }
                };
        };
    }
}
namespace ad {
    template <class VALUE_T>
    class gvalue {
        public:
            typedef VALUE_T value_t;
            typedef VALUE_T active_t;
            typedef VALUE_T passive_t;
            typedef VALUE_T type;
            typedef void derivative_t;
            typedef void tape_t;
            static const bool is_ad_type = false;
            static const bool is_adjoint_type = false;
            static const bool is_tangent_type = false;
            static const int order = 0;
    };
    template <class DcoT>
    struct mode : public ad::gvalue<DcoT> {};

    template <class T>
    class gt1s {
        public:
            typedef T value_t;
            typedef T derivative_t;
            typedef typename mode<value_t>::passive_t passive_t;
            typedef ad::internal::ts_data<derivative_t> _data;
            typedef ad::internal::active_type<T, _data> type;
            typedef void tape_t;
            static const bool is_ad_type = true;
            static const bool is_adjoint_type = false;
            static const bool is_tangent_type = true;
            static const int order = ad::mode<T>::order + 1;
    };
    template <class T>
    class ga1s {
        public:
            typedef T value_t;
            typedef T derivative_t;
            typedef typename mode<value_t>::passive_t passive_t;
            typedef ad::internal::blob_tape<derivative_t> tape_t;
            static tape_t *global_tape;
            typedef typename tape_t::template single_tape_data<global_tape> _data;
            typedef ad::internal::active_type<T, _data> type;
            typedef ad::internal::tape_options tape_options_t;
            typedef ad::helper::callback_object_base<tape_t> callback_object_t;
            typedef ad::helper::userdata_object_base<type, tape_t> userdata_object_t;
            typedef ad::helper::external_adjoint_object_base<type, tape_t> external_adjoint_object_t;
            typedef external_adjoint_object_t efo_t;
            static const bool is_ad_type = true;
            static const bool is_adjoint_type = true;
            static const bool is_tangent_type = false;
            static const int order = ad::mode<T>::order + 1;
    };
    template <typename T> typename ga1s<T>::tape_t *ga1s<T>::global_tape;
}
namespace ad {
    template<typename tape_t>
    static size_t size_of(const tape_t *tape) {
        return tape->_get_tape_memory();
    }
    template <typename AD_TAPE_REAL>
    struct mode<ad::internal::active_type<AD_TAPE_REAL, typename ad::gt1s<AD_TAPE_REAL>::_data> > : public ad::gt1s<AD_TAPE_REAL> {};
    template <typename AD_TAPE_REAL>
    struct mode<ad::internal::active_type<AD_TAPE_REAL, typename ad::ga1s<AD_TAPE_REAL>::_data> > : public ad::ga1s<AD_TAPE_REAL> {};
}

#endif 
