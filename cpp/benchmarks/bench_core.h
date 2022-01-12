#ifndef BENCH_CORE_H_
#define BENCH_CORE_H_

#include "bench_model.h"

#include "memilio/utils/logging.h"
#include "memilio/math/adapt_rk_all.h"
#include "memilio/math/vadapt_rk.h"
#include "memilio/math/vadapt_rk_all.h"
#include "memilio/math/vadapt_rk_opt.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/random_number_generator.h"

#include <benchmark/benchmark.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

namespace mio {

namespace benchmark {

    namespace detail {
        inline void reset_time(double& t, double& dt, const double t_init, const double dt_init) {
            t       = t_init;
            dt      = dt_init;
        }

        template <class dtype>
        void read(dtype& data, std::ifstream& stream) {
            std::string reader;
            if (!std::getline(stream, reader)) abort();
            std::stringstream parser;
            parser << reader;
            parser >> data;
        }

        void read(double& data, std::ifstream& stream) {
            std::string reader;
            if (!std::getline(stream, reader)) abort();
            if (reader.compare("max") == 0) {
                data = std::numeric_limits<double>::max();
            }
            else if (reader.compare("min") == 0) {
                data = std::numeric_limits<double>::min();
            }
            else {
                std::stringstream parser;
                parser << reader;
                parser >> data;
            }
        }

        void read(std::vector<double>& data, std::ifstream& stream) {
            std::string reader;
            std::stringstream parser;
            if (!std::getline(stream, reader)) abort();
            parser << reader;
            for (double& d : data) {
                parser >> d;
            }
        }
    }

    template <class Model=model::SecirAgeres, size_t num_agegroups=1>
    struct default_simulation_init {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t, double& t_max, double& dt
        ) {
            t       = 0;
            t_max   = 100;
            dt      = 0.5;
            abs_tol = 1e-10;
            rel_tol = 1e-5;
            dt_min  = std::numeric_limits<double>::min();
            dt_max  = std::numeric_limits<double>::max();
            return Model()(num_agegroups);
        }
    };

    template <char const *filepath, class Model=model::SecirAgeres>
    struct simulation_file_init {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t, double& t_max, double& dt
        ) {
            std::ifstream file(filepath);
            if (!file.good()) abort();
            size_t num_agegroups;
            detail::read(num_agegroups, file);
            detail::read(t, file);
            detail::read(t_max, file);
            detail::read(dt, file);
            detail::read(abs_tol, file);
            detail::read(rel_tol, file);
            detail::read(dt_min, file);
            detail::read(dt_max, file);
            file.close();
            return Model()(num_agegroups);
        }
    };

    template <class Integrator, class Initializer=default_simulation_init<>>
    void simulation(::benchmark::State& state) {
        mio::set_log_level(mio::LogLevel::off);
        double abs_tol, rel_tol, dt_min, dt_max, t, t_max, dt;

        auto model = Initializer()(abs_tol, rel_tol, dt_min, dt_max, t, t_max, dt);

        for (auto _ : state) {
            // This code gets timed
            std::shared_ptr<mio::IntegratorCore> I = std::make_shared<Integrator>(abs_tol, rel_tol, dt_min, dt_max);
            simulate(t, t_max, dt, model, I);
        }
    }

    template <class Model=model::SecirAgeres>
    struct default_integrator_step_init {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t_init, double& dt_init,
            mio::DerivFunction& f, Eigen::VectorXd& yt, Eigen::VectorXd& ytp1
        ) {
            typename Model::type model = Model()(1);
            t_init  = 50;
            dt_init = 1;
            abs_tol = 1e-10;
            rel_tol = 1e-5;
            dt_min  = std::numeric_limits<double>::min();
            dt_max  = std::numeric_limits<double>::max();
            f       = [model](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f) { model.eval_right_hand_side(y_f, y_f, t_f, dydt_f); };
            yt      = Eigen::VectorXd::Zero(model.get_initial_values().size());
            ytp1    = Eigen::VectorXd::Zero(yt.size());
            // set values for yt
            std::vector<double> yt_vals({6377.873644, 35.249156, 30.029611, 182.145865, 66.153059, 79.530621, 3069.383604, 159.634440});
            assert(yt.size() <= yt_vals.size());
            for (int i = 0; i < yt.size(); i++) {
                yt[i] = yt_vals[i];
            }
            return model;
        }
    };

    template <char const *filepath, class Model=model::SecirAgeres>
    struct integrator_step_file_init {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t_init, double& dt_init,
            mio::DerivFunction& f, Eigen::VectorXd& yt, Eigen::VectorXd& ytp1
        ) {
            std::ifstream file(filepath);
            if (!file.good()) abort();
            size_t num_agegroups;
            detail::read(num_agegroups, file);
            typename Model::type model = Model()(num_agegroups);
            detail::read(t_init, file);
            detail::read(dt_init, file);
            detail::read(abs_tol, file);
            detail::read(rel_tol, file);
            detail::read(dt_min, file);
            detail::read(dt_max, file);
            std::vector<double> yt_vals(model.get_initial_values().size());
            detail::read(yt_vals, file);
            file.close();
            f     = [model](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f) { model.eval_right_hand_side(y_f, y_f, t_f, dydt_f); };
            yt    = Eigen::VectorXd::Zero(yt_vals.size());
            ytp1  = Eigen::VectorXd::Zero(yt_vals.size());
            for (int i = 0; i < yt.size(); i++) {
                yt[i] = yt_vals[i];
            }
            return model;
        }
    };

    template <class Integrator, class Initializer=default_integrator_step_init<>>
    void integrator_step(::benchmark::State& state) {
        mio::set_log_level(mio::LogLevel::off);
        double abs_tol, rel_tol, dt_min, dt_max, t_init, dt_init;
        mio::DerivFunction f;
        Eigen::VectorXd yt, ytp1;

        auto model = Initializer()(abs_tol, rel_tol, dt_min, dt_max, t_init, dt_init, f, yt, ytp1);

        auto I = Integrator(abs_tol, rel_tol, dt_min, dt_max);

        double t, dt;
        for (auto _ : state) {
            // This code gets timed
            detail::reset_time(t, dt, t_init, dt_init);
            I.step(f, yt, t, dt, ytp1);
        }
    }
} // namespace benchmark

} // namespace mio

#endif