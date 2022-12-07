This folder provides three spatially resolved simulations. There are two graph-ODE (or metapopulation) models using one ODE model for each county and realizing inter-county mobility via a graph approach. An agent-based model (abm) using German statistical data is also included.

- 2020_npis_wildtype: Focus on a SECIR model using parameters for Sars-CoV-2 wild type variant and
implementing static nonpharmaceutical interventions (NPIs) as well as dynamic NPIs. Dynamic NPIs
get into play once predefined incidence thresholds (50 and 200) are exceeded.

- 2021_vaccination_delta: Extending the model of 2020_npis_wildtype by vaccination and reinfection and
considering the effect of vaccination in combination with the lifting of NPIs during the arrival of Delta.

- abm: An agent-based model using statistical data taken from Destatis (https://www-genesis.destatis.de/genesis/online?operation=statistic&levelindex=0&levelid=1627908577036&code=12211#abreadcrumb).