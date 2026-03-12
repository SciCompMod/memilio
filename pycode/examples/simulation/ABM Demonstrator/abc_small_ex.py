import os
import pyabc.external
import random
import pandas as pd
import numpy as np

import abm_dem_munich_abc_small_ex

if __name__ == "__main__":

    discrete_domain_1 = np.arange(1, 50)

    prior = pyabc.Distribution(kappa=pyabc.RV("uniform", 0, 10),
                               init_e=pyabc.RV("uniform", 0, 0.01),
                               damp_time=pyabc.RV("rv_discrete", values=(
                                   discrete_domain_1, [1 / 49] * 49)),
                               damp_lvl=pyabc.RV("uniform", 0, 1))

    def model(pars):
        curr_seed = random.randint(100000000, 999999999)
        abm_dem_munich_abc_small_ex.run_abm_simulation(curr_seed, pars["kappa"], pars["init_e"],
                                                       pars["damp_time"], pars["damp_lvl"])

        try:
            df = pd.read_csv(
                f"output_small_ex/{curr_seed}_comps.csv", delimiter=" ", index_col=False)
            df['prev'] = df['Ins'] + df['Isy'] + df['Isev'] + df['Icri']
            df = df[[(x <= 1343) for x in df["t"]]]
            df['week'] = np.floor(df['t'] / 24 / 7)

            os.remove(f'output_small_ex/{curr_seed}_mapping.txt')

            result = df.groupby(['week'])['prev'].mean().to_list()
        except pd.errors.EmptyDataError:
            result = [0] * 8

        return {"prev": result, "seed": curr_seed}

    def distance(x, x0):
        return np.linalg.norm(np.array(x["prev"]) - np.array(x0["prev"]))

    transition = pyabc.AggregatedTransition(
        mapping={
            'kappa': pyabc.MultivariateNormalTransition(),
            'init_e': pyabc.MultivariateNormalTransition(),
            'damp_time': pyabc.DiscreteJumpTransition(domain=discrete_domain_1, p_stay=0.7),
            'damp_lvl': pyabc.MultivariateNormalTransition()
        }
    )

    sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=48)

    abc = pyabc.ABCSMC(model, prior, distance, transitions=transition,
                       population_size=300, sampler=sampler)

    db_path = os.path.join('output_small_ex', "abm.db")
    # generated using seed = 123456789, kappa = 0.5, init_e = 0.002, damp_time = 15, damp_lvl = 0.2
    observed_prev = [8.70, 30.36, 79.89, 57.60, 23.39, 12.05, 8.75, 4.14]
    abc.new("sqlite:///" + db_path, {"prev": observed_prev})

    history = abc.run(minimum_epsilon=0, max_nr_populations=20)
