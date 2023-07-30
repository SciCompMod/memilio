# ABM Object 23-07-12

#Revisit  
Participants: [[Sascha Korf (DLR-SC-HPC)]]
Tags: [[PANDEMOS]]

---
- no progress due to vacations
  - will be worked on this week
- probably hdf5 files (but open in theory)

Data contained in Output (for single run):
- Locations Lookup
  | Field      | Description                             |
  |------------|-----------------------------------------|
  | LocationID | ID used in the movement vectors         |
  | Latitude   | Latitude of the location's coordinates  |
  | Longitude  | Longitude of the location's coordinates |

- Agent Lookup(?)
  | Field    | Description                                                                               |
  |----------|-------------------------------------------------------------------------------------------|
  | AgentID  | ID of the Agent                                                                           |
  | (HomeID) | Location ID of the Agents home/household                                                  |
  | Groups   | List of stratification groups (socio-economic) the Agent belongs to (i.e. age16-30, male) | 

- Movement Data
  | Field          | Description                                                                           |
  |----------------|---------------------------------------------------------------------------------------|
  | TripID         | ID of the trip/movement vector                                                        |
  | AgentID        | ID of the Agent doing the trip                                                        |
  | (HomeID)       | ID of the Household (or Location ID of the agents home?)                              |
  |                |                                                                                       |
  | StartLocation  | Location ID of the start                                                              |
  | EndLocation    | Location ID of the end                                                                |
  |                |                                                                                       |
  | StartTime      | Timestamp of the start                                                                |
  | EndTime        | Timestamp of the end                                                                  |
  |                |                                                                                       |
  | TransportMode  | mode of transportation                                                                |
  | Activity       | Activity at the end (reason for trip)                                                 |
  | InfectionState | Compartment(s) the Agent belongs to while on trip (infected, naive vaccination, etc.) |

---

Related:
