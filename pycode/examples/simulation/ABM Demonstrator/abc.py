import os
import pyabc.external
import random
import pandas as pd
import numpy as np

import abm_dem_munich_abc

if __name__ == "__main__":

    discrete_domain_1 = np.arange(1, 50)

    prior = pyabc.Distribution(kappa=pyabc.RV("uniform", 0, 10),
                               init_e=pyabc.RV("uniform", 0, 0.00013),
                               damp_time=pyabc.RV("rv_discrete", values=(
                                   discrete_domain_1, [1 / 49] * 49)),
                               damp_lvl=pyabc.RV("uniform", 0, 1))

    def model(pars):
        curr_seed = random.randint(100000000, 999999999)
        abm_dem_munich_abc.run_abm_simulation(curr_seed, pars["kappa"], pars["init_e"],
                                              pars["damp_time"], pars["damp_lvl"])

        try:
            df = pd.read_csv(
                f"output/{curr_seed}_comps.csv", delimiter=" ", index_col=False)
            df['prev'] = df['Ins'] + df['Isy'] + df['Isev'] + df['Icri']
            df = df[[(x <= 1343) for x in df["t"]]]
            df['week'] = np.floor(df['t'] / 24 / 7)

            os.remove(f'output/{curr_seed}_mapping.txt')

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

    sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=12)

    abc = pyabc.ABCSMC(model, prior, distance, transitions=transition,
                       population_size=300, sampler=sampler)

    db_path = os.path.join('output', "abm.db")
    observed_prev = [329, 2350, 6136, 6279, 4679,
                     3357, 2421, 1757]  # from Contento et al
    abc.new("sqlite:///" + db_path, {"prev": observed_prev})

    history = abc.run(minimum_epsilon=0, max_nr_populations=15)
