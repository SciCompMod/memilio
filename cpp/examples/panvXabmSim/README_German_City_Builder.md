# German City Builder for Agent-Based Modeling

This enhanced city builder creates representative German cities for epidemiological agent-based modeling simulations. The builder uses authentic German demographic and infrastructure data to generate realistic city structures at any scale.

## Features

- **Realistic Demographics**: Uses 2023 German census data for age distribution
- **Authentic Infrastructure**: Based on German federal statistics for workplaces, schools, healthcare, and retail
- **Scalable Design**: Works from small towns (1,000 people) to large cities (1,000,000+)
- **Evidence-Based**: All parameters sourced from official German statistics

## German Demographic Data Sources

### Population and Demographics

- **German Federal Statistical Office (Destatis) 2023**: Population by age groups
- **Age Distribution**: 0-4 (4.7%), 5-14 (9.6%), 15-34 (24.5%), 35-59 (33.8%), 60-79 (21.4%), 80+ (6.0%)

### Households

- **Destatis Mikrozensus 2023**: Household size distribution
- **Average Household Size**: 1.95 people per household
- **Distribution**: 43.2% single-person, 33.6% two-person, 11.7% three-person, 9.5% four-person, 2.0% five+ person

### Employment and Workplaces

- **Destatis Labor Force Statistics 2023**: Employment rates by age group
- **Overall Employment Rate**: 63.8% of population aged 15-64
- **People per Workplace**: 25 (average employees per workplace)
- **Age-Specific Employment**: 15-34 (78%), 35-59 (85%), 60-79 (32%)

### Education

- **Destatis Education Statistics 2023**: School attendance and infrastructure
- **Students per Elementary School**: 180
- **Students per Secondary School**: 450
- **Elementary to Secondary Ratio**: 3.5:1
- **School Attendance**: 5-14 years (98%), 15-34 years (45% in education/training)

### Healthcare

- **Destatis Healthcare Statistics 2023**: Hospital and ICU bed capacity
- **Hospital Beds**: 8 per 1,000 people
- **ICU Beds**: 0.8 per 1,000 people

### Retail and Services

- **German Trade Association (HDE) 2023**: Retail structure
- **Grocery Stores**: 1 per 2,000 people
- **Pharmacies**: 1 per 3,500 people
- **General Stores**: 1 per 1,500 people

### Social and Event Venues

- **German Hotel and Restaurant Association (DEHOGA) 2023**: Hospitality statistics
- **Restaurants**: 1 per 400 people
- **Bars/Pubs**: 1 per 800 people
- **Large Event Venues**: 1 per 50,000 people (concert halls, stadiums)
- **Small Event Venues**: 1 per 2,000 people (community centers, clubs)

## Usage Examples

### Small Town (10,000 people)

```cpp
CityConfig small_town;
small_town.total_population = 10000;
auto result = CityBuilder::build_world(small_town);
```

**Generated Infrastructure:**

- 5,128 households (1.95 people/household)
- 255 workplaces
- 28 elementary schools, 8 secondary schools
- 1 hospital, 1 ICU
- 5 grocery stores, 3 pharmacies, 7 general stores
- 25 restaurants, 13 bars
- 1 large event venue, 5 small event venues

### Medium City (100,000 people)

```cpp
CityConfig medium_city;
medium_city.total_population = 100000;
auto result = CityBuilder::build_world(medium_city);
```

**Generated Infrastructure:**

- 51,282 households
- 2,552 workplaces
- 277 elementary schools, 79 secondary schools
- 1 hospital, 1 ICU
- 50 grocery stores, 29 pharmacies, 67 general stores
- 250 restaurants, 125 bars
- 2 large event venues, 50 small event venues

### Large City (1,000,000 people)

```cpp
CityConfig large_city;
large_city.total_population = 1000000;
auto result = CityBuilder::build_world(large_city);
```

**Generated Infrastructure:**

- 512,821 households
- 25,520 workplaces
- 2,772 elementary schools, 792 secondary schools
- 8 hospitals, 1 ICU
- 500 grocery stores, 286 pharmacies, 667 general stores
- 2,500 restaurants, 1,250 bars
- 20 large event venues, 500 small event venues

## Technical Implementation

### Key Classes

- **CityConfig**: Configuration with population size
- **CityParameters**: German demographic and infrastructure constants
- **CityBuilder**: Main builder class with static methods
- **CityInfrastructure**: Calculated infrastructure requirements

### Age Group Assignment

Uses `std::discrete_distribution` with German age distribution probabilities for realistic demographic assignment.

### Location Assignment

- **Households**: Distributed based on German household size statistics
- **Schools**: Assigned based on age-specific attendance rates
- **Workplaces**: Assigned based on age-specific employment rates
- **Healthcare**: Randomly distributed among available facilities
- **Retail/Social**: Accessible to all residents

### Scalability

The system scales linearly with population size while maintaining realistic ratios. For extremely large populations (80M+), consider distributed simulation approaches.

## Files

- `city_parameters.h`: German demographic and infrastructure constants
- `city_builder.h`: Enhanced city builder class declaration
- `city_builder.cpp`: Implementation with German demographic logic
- `german_city_example.cpp`: Usage examples and demonstration

## Future Enhancements

- Regional variations within Germany (North vs South, urban vs rural)
- Seasonal population variations
- Economic sector-specific workplace distribution
- Transportation infrastructure integration
- Multi-generational household modeling
