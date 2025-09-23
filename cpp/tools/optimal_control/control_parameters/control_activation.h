#pragma once

#include <cmath>
#include <variant>

struct LinearActivation {
    template <typename FP>
    FP operator()(FP x) const { return x; }
};

struct SigmoidActivation {
    template <typename FP>
    FP operator()(FP x) const {
        using std::tanh;
        FP gain = 100.0;
        return 0.5 + 0.5 * tanh(gain * (x - FP(0.5))) / tanh(FP(0.5) * gain);
    }
};

using ActivationVariant = std::variant<LinearActivation, SigmoidActivation>;
