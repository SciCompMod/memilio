#include "abm/simulation.h"
#include "memilio/utils/stl_util.h"
#include "benchmark/benchmark.h"

mio::abm::Simulation make_simulation(size_t num_persons, std::initializer_list<uint32_t> seeds)
{
    auto rng = mio::RandomNumberGenerator();
    rng.seed(seeds);
    auto world      = mio::abm::World(5);
    world.get_rng() = rng;

    //create persons at home
    const auto mean_home_size    = 5.0;
    const auto min_home_size     = 1;
    auto& home_size_distribution = mio::PoissonDistribution<int>::get_instance();
    auto home                    = world.add_location(mio::abm::LocationType::Home);
    auto planned_home_size       = home_size_distribution(world.get_rng(), mean_home_size);
    auto home_size               = 0;
    for (size_t i = 0; i < num_persons; ++i) {
        if (home_size >= std::max(min_home_size, planned_home_size)) {
            home              = world.add_location(mio::abm::LocationType::Home);
            planned_home_size = home_size_distribution(world.get_rng(), mean_home_size);
            home_size         = 0;
        }

        auto age     = mio::AgeGroup(mio::UniformIntDistribution<size_t>::get_instance()(
            world.get_rng(), size_t(0), world.parameters.get_num_groups() - 1));
        auto& person = world.add_person(home, age);
        person.set_assigned_location(home);
        home_size++;
    }

    //create other locations
    for (auto loc_type : mio::enum_members<mio::abm::LocationType>()) {
        if (loc_type == mio::abm::LocationType::Home) {
            continue; //homes already created
        }

        const auto num_locs = std::max(size_t(1), num_persons / 5'000);
        std::vector<mio::abm::LocationId> locs(num_locs);
        std::generate(locs.begin(), locs.end(), [&] {
            return world.add_location(loc_type);
        });
        for (auto& person : world.get_persons()) {
            auto loc_idx =
                mio::UniformIntDistribution<size_t>::get_instance()(world.get_rng(), size_t(0), num_locs - 1);
            person.set_assigned_location(locs[loc_idx]);
        }
    }

    //infections
    for (auto& person : world.get_persons()) {
        auto prng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        //~0.5% of people are infected, large enough to have some infection activity without everyone dying
        if (mio::UniformDistribution<double>::get_instance()(prng, 0.0, 1.0) < 0.005) {
            auto state = mio::abm::InfectionState(
                mio::UniformIntDistribution<int>::get_instance()(prng, 1, int(mio::abm::InfectionState::Count) - 1));
            auto infection =
                mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                    world.parameters, mio::abm::TimePoint(0), state);
            person.add_new_infection(std::move(infection));
        }
    }

    //testing schemes
    auto sample = [&](auto v, size_t n) { //selects n elements from list v
        std::shuffle(v.begin(), v.end(), world.get_rng());
        return std::vector<typename decltype(v)::value_type>(v.begin(), v.begin() + n);
    };
    std::vector<mio::AgeGroup> ages;
    std::generate_n(std::back_inserter(ages), world.parameters.get_num_groups(), [a = 0]() mutable {
        return mio::AgeGroup(a++);
    });
    auto random_criteria = [&]() {
        auto random_ages   = sample(ages, 2);
        auto random_states = std::vector<mio::abm::InfectionState>(0);
        return mio::abm::TestingCriteria(random_ages, random_states);
    };

    world.get_testing_strategy().add_testing_scheme(
        mio::abm::LocationType::School,
        mio::abm::TestingScheme(random_criteria(), mio::abm::days(3), mio::abm::TimePoint(0),
                                mio::abm::TimePoint(0) + mio::abm::days(10), {}, 0.5));
    world.get_testing_strategy().add_testing_scheme(
        mio::abm::LocationType::Work,
        mio::abm::TestingScheme(random_criteria(), mio::abm::days(3), mio::abm::TimePoint(0),
                                mio::abm::TimePoint(0) + mio::abm::days(10), {}, 0.5));
    world.get_testing_strategy().add_testing_scheme(
        mio::abm::LocationType::Home,
        mio::abm::TestingScheme(random_criteria(), mio::abm::days(3), mio::abm::TimePoint(0),
                                mio::abm::TimePoint(0) + mio::abm::days(10), {}, 0.5));
    world.get_testing_strategy().add_testing_scheme(
        mio::abm::LocationType::SocialEvent,
        mio::abm::TestingScheme(random_criteria(), mio::abm::days(3), mio::abm::TimePoint(0),
                                mio::abm::TimePoint(0) + mio::abm::days(10), {}, 0.5));

    return mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));
}

/**
 * Benchmark for the ABM simulation.
 * @param num_persons Number of persons in the simulation.
 * @param seeds Seeds for the random number generator. 
 */
void abm_benchmark(benchmark::State& state, size_t num_persons, std::initializer_list<uint32_t> seeds)
{
    mio::set_log_level(mio::LogLevel::warn);

    for (auto&& _ : state) {
        state.PauseTiming(); //exclude the setup from the benchmark
        auto sim = make_simulation(num_persons, seeds);
        state.ResumeTiming();
        sim.advance(sim.get_time() + mio::abm::days(10));
    }
}

//Measure ABM simulation run time with different sizes and different seeds.
//Fixed RNG seeds to make runs comparable. When there are code changes, the simulation will still
//run differently due to different sequence of random numbers being drawn. But for large enough sizes
//RNG should average out, so runs should be comparable even with code changes.
//We run a few different benchmarks to hopefully catch abnormal cases. Then seeds may
//have to be adjusted to get the benchmark back to normal.
//For small sizes (e.g. 10k) extreme cases are too likely, i.e. infections die out
//or overwhelm everything, so we don't benchmark these. Results should be mostly transferrable.
BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_50k, 50000, {14159265u, 35897932u})->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_100k, 100000, {38462643u, 38327950u})->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_200k, 200000, {28841971u, 69399375u})->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
