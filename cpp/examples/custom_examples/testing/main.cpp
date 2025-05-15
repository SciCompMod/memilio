#include "memilio/epidemiology/adoption_rate.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/epidemiology/lct_populations.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"

// #include <iostream>
// #include <vector>

// using namespace mio;

// int main() {
//     // Create AgeGroup indices
//     AgeGroup age0(0);
//     AgeGroup age1(1);
//     AgeGroup age2(2);

//     // Test comparison
//     std::cout << "age0 == age1? " << (age0 == age1) << "\n"; // 0 (false)
//     std::cout << "age2 > age1? " << (age2 > age1) << "\n";   // 1 (true)

//     // Arithmetic operations (from OperatorAdditionSubtraction, etc.)
//     AgeGroup age_sum = AgeGroup(static_cast<size_t>(age1 + AgeGroup(2)));
//     std::cout << "age1 + 2 = " << static_cast<size_t>(age_sum) << "\n";

//     // Test multidimensional Index (e.g., AgeGroup + dummy category)
//     struct DummyCategory : public Index<DummyCategory> {
//         DummyCategory(size_t val) : Index<DummyCategory>(val) {}
//     };

//     using MultiIndex = Index<AgeGroup, DummyCategory>;
//     MultiIndex mi(AgeGroup(3), DummyCategory(5));

//     std::cout << "MultiIndex<AgeGroup, DummyCategory>: Age=" 
//               << static_cast<size_t>(get<AgeGroup>(mi)) 
//               << ", Dummy=" << static_cast<size_t>(get<DummyCategory>(mi)) << "\n";

//     // Reduce/Extend Index
//     Index<AgeGroup, DummyCategory> super_index(AgeGroup(4), DummyCategory(7));
//     Index<AgeGroup> sub_index = reduce_index<Index<AgeGroup>>(super_index);
//     std::cout << "Reduced index (AgeGroup only): " << static_cast<size_t>(get<AgeGroup>(sub_index)) << "\n";

//     auto extended_index = extend_index<Index<AgeGroup, DummyCategory>>(sub_index, 99);
//     std::cout << "Extended index: AgeGroup=" << static_cast<size_t>(get<AgeGroup>(extended_index)) 
//               << ", Dummy=" << static_cast<size_t>(get<DummyCategory>(extended_index)) << "\n";

//     // Demonstrate serialization interface (stubbed as pseudo-code)
//     // IOContext io;
//     // age0.serialize(io);
//     // auto result = AgeGroup::deserialize(io);

//     return 0;
// }








#include "memilio/math/eigen.h"
#include "memilio/utils/type_safe.h"
#include "memilio/math/matrix_shape.h"
#include "memilio/math/smoother.h"
#include "memilio/math/matrix_shape.h"

#include <iostream>

// using namespace mio;

// void print_matrix(const std::string& label, const Eigen::MatrixXd& m) {
//     std::cout << label << ":\n" << m << "\n\n";
// }

// int main() {
//     // Example 1: Two dampings at the SAME LEVEL (additive combination)
//     {
//         SquareDampings school_closure(1); // 1x1 matrix (single group)
        
//         // Add 50% reduction at t=0 days (level 0)
//         school_closure.add(0.5, DampingLevel(0), DampingType(0), SimulationTime(0.0));
        
//         // Add 30% reduction at t=10 days (same level 0)
//         school_closure.add(0.3, DampingLevel(0), DampingType(0), SimulationTime(10.0));

//         std::cout << "=== Example 1: Same Level Additive Combination ===\n";
//         print_matrix("Immediate effect (t=0)", school_closure.get_matrix_at(0.0));
//         print_matrix("After second closure (t=10)", school_closure.get_matrix_at(10.0));
//         print_matrix("During transition (t=5)", school_closure.get_matrix_at(5.0));
//     }

//     // Example 2: Two dampings at DIFFERENT LEVELS (multiplicative combination)
//     {
//         SquareDampings pandemic_response(1);
        
//         // School closure: 60% reduction (level 0)
//         pandemic_response.add(0.6, DampingLevel(0), DampingType(0), SimulationTime(0.0));
        
//         // Vaccination: 30% reduction (level 1)
//         pandemic_response.add(0.3, DampingLevel(1), DampingType(1), SimulationTime(0.0));

//         std::cout << "\n=== Example 2: Different Level Multiplicative Combination ===\n";
//         print_matrix("Combined effect (t=0)", pandemic_response.get_matrix_at(0.0));
//     }

//     // Example 3: Complex Scenario with Smooth Transitions
//     {
//         SquareDampings covid_policies(1);
        
//         // Initial lockdown: 70% reduction at t=0
//         covid_policies.add(0.7, DampingLevel(0), DampingType(0), SimulationTime(0.0));
        
//         // Partial reopening: 40% reduction at t=14 days
//         covid_policies.add(0.4, DampingLevel(0), DampingType(1), SimulationTime(14.0));
        
//         // New variant response: 20% reduction at t=21 (different level)
//         covid_policies.add(0.2, DampingLevel(1), DampingType(1), SimulationTime(21.0));

//         // New variant response: 20% reduction at t=21 (different level)
//         covid_policies.add(0.5, DampingLevel(2), DampingType(1), SimulationTime(28.0));

//         std::cout << "\n=== Example 3: Time-Based Transitions ===";
//         for(double t = 0.0; t <= 28.0; t += 7.0) {
//             std::cout << "\n\nDay " << t << ":\n";
//             print_matrix("Effective damping", covid_policies.get_matrix_at(t));
//         }
//     }

//     return 0;
// }


// Dampings vom gleichen Typ werden addiert.
// Dampings von verschiedenen Leveln werden multipliziert.
// Hinzufügen eines vorherigen Dampings, verändert nur den Wert.

// damping.add(0.7, DampingLevel(0), DampingType(0), SimulationTime(0.0));
// damping.add(0.4, DampingLevel(0), DampingType(1), SimulationTime(14.0));
// damping.add(0.2, DampingLevel(1), DampingType(1), SimulationTime(21.0));
// damping.add(0.5, DampingLevel(2), DampingType(1), SimulationTime(28.0));
// damping.add(0.3, DampingLevel(2), DampingType(1), SimulationTime(28.0));

// damping.get_matrix_at(28.0) = 1 - (1 - (0.7 + 0.4)) * (1 - 0.2) * (1 - 0.3)









#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/utils/random_number_generator.h"
#include <iostream>

using namespace mio;

void print_matrix(const std::string& label, const Eigen::MatrixXd& m) {
    std::cout << label << ":\n" << m << "\n\n";
}

int main() {
    // // Example 1: Age-Specific Contact Reduction
    // {
    //     // 3 age groups: children, adults, elderly
    //     Eigen::VectorXd groups(3);
    //     groups << 1.0, 0.5, 0.2; // Full effect on children, half on adults, 20% on elderly
        
    // DampingSampling school_closure(
    //     UniformDistribution<double>(std::uniform_real_distribution<double>(0.4, 0.6)),
    //     DampingLevel(0), DampingType(0), SimulationTime(0.0),
    //     {0}, groups);

    // std::vector<SquareDampings> contact_matrices(1, SquareDampings(3));

    // auto matrix_func = static_cast<decltype(make_contact_damping_matrix)*>(&make_contact_damping_matrix);
    // apply_dampings(contact_matrices, std::vector<DampingSampling>{school_closure}, matrix_func);

    //     DampingSampling school_closure(
    //         UniformDistribution<double>(0.4, 0.6), // 40-60% effectiveness
    //         DampingLevel(0), DampingType(0), SimulationTime(0.0),
    //         {0}, // Apply to first contact matrix
    //         groups
    //     );

    //     std::vector<SquareDampings> contact_matrices(1, SquareDampings(3));
        
    //     // Apply with mean value (50%)
    //     apply_dampings(contact_matrices, {school_closure}, make_contact_damping_matrix);
        
    //     std::cout << "=== Example 1: Age-Specific Contact Reduction ===\n";
    //     std::cout << "Mean estimate (50% base reduction):\n";
    //     print_matrix("School contacts", contact_matrices[0].get_matrix_at(0.0));

    //     // Draw actual sample (e.g., 55%)
    //     school_closure.draw_sample();
    //     contact_matrices[0].clear();
    //     apply_dampings(contact_matrices, {school_closure}, make_contact_damping_matrix);
        
    //     std::cout << "Sampled value (" << (100*school_closure.get_value()) << "% effectiveness):\n";
    //     print_matrix("Adjusted school contacts", contact_matrices[0].get_matrix_at(0.0));
    // }

    // // Example 2: Combined Mobility Restrictions
    // {
    //     // Fixed travel restriction: 30% reduction
    //     SquareDampings fixed_travel_ban(3);
    //     fixed_travel_ban.add(0.3, DampingLevel(0), DampingType(1), SimulationTime(0.0));

    //     // Random mobility reduction: 0-20% additional reduction
    //     DampingSampling mobility_sampling(
    //         UniformDistribution<double>(0.0, 0.2),
    //         DampingLevel(1), DampingType(1), SimulationTime(7.0),
    //         {0}, // Same matrix index
    //         Eigen::VectorXd::Constant(3, 1.0)
    //     );

    //     std::vector<SquareDampings> mobility_matrices(1, SquareDampings(3));
    //     mobility_matrices[0] = fixed_travel_ban;

    //     // Apply random component
    //     mobility_sampling.draw_sample();
    //     apply_dampings(mobility_matrices, {mobility_sampling}, make_contact_damping_matrix);

    //     std::cout << "\n=== Example 2: Layered Mobility Reductions ===\n";
    //     std::cout << "Fixed reduction: 30%\n";
    //     std::cout << "Sampled additional reduction: " << (100*mobility_sampling.get_value()) << "%\n";
    //     print_matrix("Total mobility damping", mobility_matrices[0].get_matrix_at(7.0));
    // }

    // // Example 3: Time-Dependent Workplace Restrictions
    // {
    //     DampingSampling work_schedule(
    //         NormalDistribution<double>(0.4, 0.1), // 40% ±10% reduction
    //         DampingLevel(0), DampingType(2), SimulationTime(14.0),
    //         {1}, // Second contact matrix (workplaces)
    //         (Eigen::VectorXd(3) << 0.1, 0.8, 0.4).finished() // Strong effect on adults
    //     );

    //     std::vector<SquareDampings> workplaces(2, SquareDampings(3)); // [0]=schools, [1]=workplaces
        
    //     // Sample and apply
    //     work_schedule.draw_sample();
    //     apply_dampings(workplaces, {work_schedule}, make_contact_damping_matrix);

    //     std::cout << "\n=== Example 3: Workplace Restrictions ===\n";
    //     std::cout << "Sampled workplace reduction: " << (100*work_schedule.get_value()) << "%\n";
    //     print_matrix("Work contact matrix", workplaces[1].get_matrix_at(14.0));
    //     print_matrix("School contact matrix (unaffected)", workplaces[0].get_matrix_at(14.0));
    // }

    return 0;
}