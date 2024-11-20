# TODO
#
# for each scenario we get one scenario.json:
# - read in scenario.json
# - read in model.json (sevirvvs.json)
# - set parameters from scenario.json
#      - from model.json: for each group, check if group is of category AgeGroup and accumulate to find out number of AgeGroups
#      - for each parameterID check corresponding model parameter
#      - then check groupID for corresponding model group
#      - set value(s)
# - set NPIs from scenario.json
#      - for each interventionID check intervention name and description
#      - apply intervention to the model with values
# - create graph model
# - initialize model as usual with data from RKI or other data
# - run simulation
# - post results to ESID backend with new API
#
