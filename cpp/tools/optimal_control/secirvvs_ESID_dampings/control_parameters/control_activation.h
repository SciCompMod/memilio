#pragma once

#include <cmath>

enum class ControlActivation
{
    Linear,
    Sigmoid
};

struct ControlActivationFunction {
    ControlActivationFunction(ControlActivation activation)
        : m_activation(activation)
    {
    }

    template <typename FP>
    FP operator()(FP x) const
    {
        switch (m_activation) {
        case ControlActivation::Linear:
            return linear(x);
        case ControlActivation::Sigmoid:
            return sigmoid(x);
        }
        return x;
    }

private:
    ControlActivation m_activation;

    // Linear activation: f(x) = x
    template <typename FP>
    static FP linear(FP x)
    {
        return x;
    }

    // Sigmoid activation: f(x) = 0.5 + 0.5 * tanh(gain*(x-0.5)) / tanh(0.5 * gain)
    template <typename FP>
    static FP sigmoid(FP x)
    {
        using std::tanh;
        FP gain = 100.0;
        return 0.5 + 0.5 * tanh(gain * (x - 0.5)) / tanh(0.5 * gain);
    }
};
