MEmilio INSIDe Demonstrator
===========================

Input
------

As input for the demonstrator an area list should be provided as txt file. The txt file should contain the following columns

- id: the area ids
- inhabitants: inhabitants per area
- type: area type

Possible area types are 

- recreational
- shopping_business
- university
- residential
- mixed

Output
------

The demonstrator returns two txt files as output. 

The text file location_mapping.txt contains a table that maps the input area ids to the corresponding abm location ids. 
The abm location ids have the form xxyyy where the first two digits xx specify the location's type and the following number yyy 
specifies the location's index which is a consecutive number for all locations.
The digits xx for the location types can be the following

- 00 - Home
- 01 - School
- 02 - Work
- 03 - SocialEvent
- 04 - BasicsShop
- 05 - Hospital
- 06 - ICU
- 10 - Cemetery

It should be remarked that there is no one-to-one mapping from the input areas to the abm locations and one input area can be mapped to multiple abm locations.

The text file output.txt contains an agent level output for every location and time point. One row represents one location and the table has the following columns

- LocationId
- Number of output timesteps (tn)
- Timestep 1 (t1)
- Number of agents at location at t1 (am-t1)
- Agent id 1 (a1-t1)
- Time since transmission for a1-t1
- Agent id 2 (a2-t1)
- Time since transmission for a2-t1
- ...
- Agent id m (am-t1)
- Time since transmission for am-t1
- Timestep 2 (t2)
- Number of agents at location at t2 (am-t2)
- Agent id 1 (a1-t2)
- Time since transmission for a1-t2
- Agent id 2 (a2-t2)
- Time since transmission for a2-t2
- ...
- Agent id m (am-t2)
- Time since transmission for am-t2
- ...
- Timestep n (tn)
- Number of agents at location at tn (am-tn)
- Agent id 1 (a1-tn)
- Time since transmission for a1-tn
- Agent id 2 (a2-tn)
- Time since transmission for a2-tn
- ...
- Agent id m (am-tn)
- Time since transmission for am-tn

The time since transmission as well as the timesteps are given in hours. If an agent is susceptible, recovered or dead, his time since transmission is set to -1.

Additionally to the txt output, there is a console output. The console output is a table with the number of agents per timestep and infection state.
The timesteps are given in days here instead of hours. The infection states are

- Susceptible
- Exposed
- InfectedNoSymptoms (Carrier)
- InfectedSymptoms
- InfectedSevere
- InfectedCritical
- Recovered
- Dead
