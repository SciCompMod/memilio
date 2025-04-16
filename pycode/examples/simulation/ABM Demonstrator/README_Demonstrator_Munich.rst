MEmilio INSIDe Demonstrator for Munich
=======================================

Input
------

Inputs for the Munich demonstrator are:

- person.csv: A file containing all persons that should be simulated. The file needs to have the following columns:
    * puid: PersonId
    * age: Age of the person in years.
    * home_zone: Traffic cell/zone of person's Home location.
    * home_id: Id of person's Home location.
    * shop_zone: Traffic cell/zone of person's Shop location.
    * shop_id: Id of person's Shop location.
    * event_zone: Traffic cell/zone of person's Event location.
    * event_id: Id of person's Event location.
    * work_zone: Traffic cell/zone of person's Work location. Should be set to -2 for non-work-aged persons.
    * work_id: Id of person's Work location. Should be set to -2 for non-work-aged persons.
    * school_zone: Traffic cell/zone of person's School location. Should be set to -2 for non-work-aged persons.
    * school_id: Id of person's School location. Should be set to -2 for non-work-aged persons.

- hospitals.csv : A file containing all hospitals used in the simulation. The file needs to have the following columns:
    * hospital_id: Id of the hospital
    * hospital_zone: Traffic cell/zone the hospital is in
    * beds: Number of beds (capacity) of the hospital
    * icu_beds: Number of ICU beds (capacity) of the hospital. Can be 0.

- parameter_table.csv: A file containing all parameters associated to the transmission model stratified by age group.

- Verschnitt_DLR_TAN_Rep.shp: Shape file containing the traffic zones and the wastewater areas. The file needs to have the following columns:
    * id_n: Traffic zone id
    * ID_TAN: Wastewater area id the traffic zone with id id_n is in

Output
------

The demonstrator offers the opportunity to return various files as output.

Mapping files:

- {sim_num}_mapping.txt: Contains a table that maps every ABM location to the corresponding traffic cell/zone. 

- {sim_num}_mapping_tan.txt: Contains a table that maps every ABM location to the corresponding wastewater area. 

- {sim_num}_mapping_tan_locs.txt: List of pairs (ABM location, Wastewater area). Is used to set the wastewater area attribute in the Location object during the initialization.

The ABM locations have the form xxyyy where the first two digits xx specify the location's type and the following number yyy specifies the location's index which is a consecutive number for all locations.
The digits xx for the location types can be the following:

    * 00 - Home
    * 01 - School
    * 02 - Work
    * 03 - SocialEvent
    * 04 - BasicsShop
    * 05 - Hospital
    * 06 - ICU
    * 10 - Cemetery

Cumulative output:

- {sim_num}_comps.csv: Output is a table with the number of agents per timestep and infection state.
The timesteps are given in hours. The infection states are

- Susceptible
- Exposed
- InfectedNoSymptoms (Carrier)
- InfectedSymptoms
- InfectedSevere
- InfectedCritical
- Recovered
- Dead

Agent-based output:

- {sim_num}_infection_paths.txt: Output is a table with the number of time steps spent in every infection state during the simulation for every agent.

There is also the option to output various HDF5 files that contain time-resolved information on transmission and location for every agent. The files that can be output are:

- {sim_num}_output_v1.h5: Has an own HDF5 group per agent. Every group (agent) contains two datasets. The first one holds the ABM Location in the format xxyyy for every time point and the second holds the time since transmission for every time point. The time since transmission is - like the timesteps - given in hours. If an agent is susceptible, recovered or dead, its time since transmission is set to -1.

- {sim_num}_output_v2.h5: Also has an own HDF5 group per agent. The first dataset is a vector where the even indices are time points and the odd indices are the time since transmission at change points. The first entries are time point and time since transmission at simulation start. If the agent has experienced an infection during the simulation, the following entries are the time point of transmission and recovery/death (if the agent recovered/died during the simulation) and the corresponding time since transmissions. The second dataset is a vector with the time points of location changes and the third dataset holds the corresponding (new) ABM Locations in the format xxyyy.

- {sim_num}_output_v3.h5: This file only has one HDF5 group 'data' which contains two datasets. The first dataset is a matrix of size #agents x 2 where entry (i, 0) is the time point of transmission and entry (i,1) the time points of recovery or death of agent i. If an agent has not been infected the values are set to 100,000. The second dataset is a matrix of size #agent x #timepoints where entry (i,t) is the ABM Location in the format xxyyy at time point t of agent i.

- {sim_num}_output_v4.h5: This file only has one HDF5 group 'data' which contains two datasets. The first dataset is a matrix of size #agents x 2 where entry (i, 0) is the time point of transmission and entry (i,1) the time points of recovery or death of agent i. If an agent has not been infected the values are set to 100,000. The second dataset is a matrix of size #agent x #timepoints where entry (i,t) is the wastewater area (integer>0) at time point t of agent i. If an agent was at a location that is not in any wastewater area, the value is set to 0.

- {sim_num}_output_v5.h5: This file only has one HDF5 group 'data' which contains two datasets. The first dataset is a matrix of size #agents x 2 where entry (i, 0) is the time point of transmission and entry (i,1) the time points of recovery or death of agent i. If an agent has not been infected the values are set to 100,000. The second dataset is a matrix of size #agent x #timepoints where entry (i,t) is the ABM LocationType (as integer) at time point t of agent i


