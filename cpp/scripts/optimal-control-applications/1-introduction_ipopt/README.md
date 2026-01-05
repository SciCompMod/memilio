# Introduction to IPOPT Optimization Examples

This directory contains various examples demonstrating how to use the [IPOPT](https://coin-or.github.io/Ipopt/) solver to solve optimization problems with different approaches to derivative computation. The examples are intended for users interested in understanding how to apply IPOPT for nonlinear optimization problems, with a focus on the trade-offs between different derivative calculation methods, such as **exact derivatives** and **automatic differentiation (AD)**.

## Goals of the Example

The primary goal of this example is to introduce how to solve nonlinear optimization problems using IPOPT in C++. In particular, the example focuses on providing gradients for the objective function and constraints using different methods:

1. **Exact Derivatives**: The first example uses manually coded exact derivatives.
2. **Forward Mode Automatic Differentiation (AD)**: The second example leverages forward-mode AD to compute exact gradients efficiently for small numbers of input variables or when many CPU cores are available.
3. **Reverse Mode Automatic Differentiation (AD)**: The third example uses reverse-mode AD, which is particularly efficient for large numbers of input variables and small numbers of outputs.

These examples show the different strategies and their impact on optimization performance, along with how to configure IPOPT options to solve the nonlinear programming (NLP) problem.
