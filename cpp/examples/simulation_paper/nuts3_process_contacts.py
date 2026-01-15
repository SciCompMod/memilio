import pandas as pd
import os
import numpy as np

cwd = os.path.dirname(__file__)
result_dir = os.path.join(cwd, "results")

# print(df.head())
# print(df.tail())
# print(df.isna().sum())

TIMESTEPS_PER_REGION = 5
TIMESTEP = 0.1
REGIONS = 400
TMAX = 60

all_times = np.arange(0, TMAX, TIMESTEP)

def sanity_check(df, strategy_name):
    times = df["Time"].values
    time_blocks = times[::TIMESTEPS_PER_REGION * REGIONS]

    if not np.array_equal(np.array(time_blocks), all_times[::TIMESTEPS_PER_REGION]):
        raise ValueError("Invalid time values detected")
    if times.size != REGIONS*TMAX*(1/TIMESTEP):
        raise ValueError("Invalid time values detected")

    idx = 0
    for time_block_idx in range(0, TMAX * 2, 1):
        for region_idx in range(0, REGIONS, 1):
            for step_in_block in range(0, TIMESTEPS_PER_REGION, 1):
                current_time = round(time_block_idx * 0.5 + step_in_block * 0.1, 1)
                if current_time != df["Time"].iat[idx]:
                    raise ValueError("Invalid time values detected")
                idx += 1

# Maybe add a check on the population



# Handle the data
def process_data(strategy_name):
    df = pd.read_csv(
        os.path.join(result_dir, f"Output_{strategy_name}.txt"),
        sep='\s+',
        header=0 
    )
    
    sanity_check(df, strategy_name)

    # infection states numbers
    total_newly_infected  = df["NewlyInfected"].sum() * TIMESTEP
    total_newly_severe    = df["NewlyInfectedSevere"].sum() * TIMESTEP
    total_newly_critical  = df["NewlyInfectedCritical"].sum() * TIMESTEP
    total_newly_dead      = df["NewlyDead"].sum() * TIMESTEP
    total_newly_recovered = df["NewlyRecovered"].sum() * TIMESTEP

    # contacts
    total_contacts = (
        df["ContactFrequency"]
        * df["Population"]
        * TIMESTEP
    ).sum()

    # population without dead
    total_population = df[df["Time"] == 0]["Population"].sum()

    total_contacts_per_day = total_contacts / TMAX
    mean_contacts_per_day = total_contacts_per_day / total_population

    # print(f"{'='*40}")
    # print(f"{strategy_name.upper():^50}")
    # print(f"{'Simulation Summary':^40}")
    # print(f"{'='*40}")
    # print(f"Total newly infected (all regions, all days):  {total_newly_infected:,.0f}")
    # print(f"Total newly severe (all regions, all days):    {total_newly_severe:,.0f}")
    # print(f"Total newly critical (all regions, all days):  {total_newly_critical:,.0f}")
    # print(f"Total newly revovered (all regions, all days): {total_newly_recovered:,.0f}")
    # print(f"Total newly dead (all regions, all days):      {total_newly_dead:,.0f}")
    # print(f"Total contacts (all regions, all days):        {total_contacts:,.0f}")
    # print(f"Total contacts per day (all regions):          {total_contacts_per_day:,.0f}")
    # print(f"Mean contacts per person per day:              {mean_contacts_per_day:.2f}")
    # print(f"Total population (without dead):               {total_population:,.0f}")
    # print(f"{'='*40}")

    return {
        "strategy": strategy_name,
        "contacts": mean_contacts_per_day,
        "infected": total_newly_infected,
        "severe": total_newly_severe,
        "critical": total_newly_critical,
        "dead": total_newly_dead,
    }

def print_latex_table(results):
    print(r"\begin{tabular}{|lrrrrl|}")
    print(r"\hline")
    print(r"\textbf{Strategy} & "
          r"\textbf{Contacts} & "
          r"\textbf{Infected} & "
          r"\textbf{Severe} & "
          r"\textbf{Critical} & "
          r"\textbf{Dead} \\")
    print(r"\hline")

    for r in results:
        print(
            f"{r['strategy']} & "
            f"{r['contacts']:.2f} & "
            f"{r['infected']:,.0f} & "
            f"{r['severe']:,.0f} & "
            f"{r['critical']:,.0f} & "
            f"{r['dead']:,.0f} \\\\"
        )

    print(r"\hline")
    print(r"\end{tabular}")

strategies = {
    "No NPIs": "open",
    "NPIs continued": "same",
    "Strict NPIs": "lockdown",
    "Dynamic NPIs": "dynamic",
}

results = []

for label, strategy in strategies.items():
    data = process_data(strategy)
    data["strategy"] = label  # overwrite with LaTeX-friendly name
    results.append(data)

print_latex_table(results)
