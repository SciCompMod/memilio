from memilio.simulation import abm

if __name__ == "__main__":
    # Path to simulation results folder. Has to have a "/" at the end e.g. /Documents/outputs/
    sim_result_file = ""
    # Folder the aggregated output should be saved to
    save_folder = sim_result_file
    # Path to person.csv file
    person_file = ""
    num_sims = 1
    start_sim = 0  # Index the first simulation starts with
    abm.calculate_infections_per_quantity(
        sim_result_file, save_folder, num_sims, start_sim, person_file)
    abm.calculate_agents_per_quantity_age_groups(
        sim_result_file, save_folder, num_sims, start_sim, person_file)
    abm.calculate_agents_per_quantity(
        sim_result_file, save_folder, num_sims, start_sim, person_file)
