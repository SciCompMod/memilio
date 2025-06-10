from memilio.simulation import abm

if __name__ == "__main__":
    sim_result_file = ""
    save_folder = sim_result_file
    num_sims = 1
    abm.calculate_infections_per_quantity(
        sim_result_file, save_folder, num_sims)
