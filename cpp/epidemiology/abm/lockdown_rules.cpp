#include "epidemiology/abm/lockdown_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/time.h"


namespace epi
{

void set_home_office(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant(1, p);
    params.get<WorkRatio>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}

void set_school_closure(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant(1, p);
    params.get<SchoolRatio>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}

void close_social_events(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant((size_t)AbmAgeGroup::Count, p);
    params.get<SocialEventRate>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}


} //namespace epi