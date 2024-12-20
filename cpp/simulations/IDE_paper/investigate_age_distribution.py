import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def get_df_daily(data_dir):
    """ Get a dataframe that contains daily new confirmed cases and daily deaths from RKI data.

    @param[in] data_dir Directory where file with RKI data is stored.
    @returns Dataframe with daily data on confirmed cases and deaths.
    """
    # Read file.
    datafile = os.path.join(data_dir, "pydata", "Germany", "cases_all_age_all_dates.json")
    df = pd.read_json(datafile)

    # Create df_daily, where daily confirmed (and daily deaths) will be stored.
    df_daily = df.copy()
    df_daily = df_daily.drop(columns=["Confirmed", "Deaths", "Recovered"])
    df_daily = df_daily[df_daily.Date != "2020-01-01"]
    df_daily.insert(2, "DailyConfirmed", value=0.)
    df_daily.insert(3, "DailyDeaths", value=0.)

    # Compute daily new confirmed cases from df (where cumulative confirmed cases are stored) for each agegroup. (Same for Deaths.)
    agegroups = df["Age_RKI"].unique()
    for agegroup in agegroups:
        df_daily.loc[df_daily["Age_RKI"] == agegroup,
                     "DailyConfirmed"] = df.loc[df["Age_RKI"] == agegroup, "Confirmed"].diff().dropna()
        df_daily.loc[df_daily["Age_RKI"] == agegroup,
                     "DailyDeaths"] = df.loc[df["Age_RKI"] == agegroup, "Deaths"].diff().dropna()

    return df_daily


def get_proportional_population_per_agegroup():
    """
    Computed the proportion of each age group compared to the total population.
    """
    # Population from Table 12411-04-02-4-B from regionalstatistik.de, data from 31.12.2020.
    population_per_agegroup = np.array(
        [3969138, 7508662, 18921292, 28666166, 18153339, 5936434])
    total_population = population_per_agegroup.sum()
    population_per_agegroup = population_per_agegroup/total_population

    return population_per_agegroup


def get_relevant_confirmed_cases(data_dir, start_date, T_IH, T_HU, T_U):
    """Gets the confirmed cases per age group for the relevant time frame that is determined by T_IH, T_HU and T_U.

    @param[in] data_dir Directory where file with RKI data is stored.
    @param[in] start_date Start date of interest. 
    @param[in] T_IH Mean stay time in Infected compartment before transitioning to Hospitalized. 
    @param[in] T_HU Mean stay time in Hospitalized compartment before transitioning to ICU.
    @param[in] T_U Mean stay time in ICU.
    """
    # Get dataframe with daily confirmed cases.
    df = get_df_daily(data_dir)

    # Considered age groups.
    agegroups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']

    # Extract relevant dates to be considered.
    df_date = df[(df["Date"] >= pd.Timestamp(start_date)-pd.Timedelta(days=T_IH+T_HU+T_U))
                 & (df["Date"] <= pd.Timestamp(start_date)-pd.Timedelta(days=T_IH+T_HU))]

    # Get total confirmed cases in considered time frame.
    totaldailyconfirmed = df_date.DailyConfirmed.sum()

    # Check the proportion of unknown age group.
    unknown_proportion = df_date.loc[df["Age_RKI"] ==
                                     "unknown", "DailyConfirmed"].sum()/totaldailyconfirmed
    if unknown_proportion > 0.001:
        print(
            f"A proportion of {unknown_proportion} of the considered cases has unknown age group.")

    # Get total daily confirmed cases in considered time frame per age group.
    dailyconfirmed_per_agegroup = []
    for agegroup in agegroups:
        dailyconfirmed_per_agegroup.append(
            df_date.loc[df["Age_RKI"] == agegroup, "DailyConfirmed"].sum())

    # Compute proportion of totaldailyconfirmed.
    dailyconfirmed_per_agegroup = dailyconfirmed_per_agegroup/totaldailyconfirmed

    return dailyconfirmed_per_agegroup


def plot(data_dir, start_dates, T_IH, T_HU, T_U, save_dir=""):
    """ Plots proportions per age groups in confirmed cases that we expect in ICU compartment at the start_dates as 
    well as the proportion per age group of the total population of Germany. 

    @param[in] data_dir Directory where file with RKI data is stored.
    @param[in] start_dates List of start dates.
    @param[in] T_IH Mean stay time in Infected compartment before transitioning to Hospitalized. 
    @param[in] T_HU Mean stay time in Hospitalized compartment before transitioning to ICU.
    @param[in] T_U Mean stay time in ICU.
    @param[in] save_dir Directory where plot will be stored.
    """

    agegroups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']

    # Get proportions per age group for June and Ocotber scenario and corresponding share of total population.
    proportions_june = get_relevant_confirmed_cases(data_dir, 
        start_dates[0], T_IH, T_HU, T_U)
    proportions_october = get_relevant_confirmed_cases(data_dir, 
        start_dates[1], T_IH, T_HU, T_U)
    population_per_agegroup = get_proportional_population_per_agegroup()

    proportions = [proportions_june,
                   proportions_october, population_per_agegroup]
    labels = ["June", "October", "Population"]
    colors = [plt.cm.viridis(x) for x in np.linspace(0., 0.9, 3)]

    # Plot.
    fig, ax = plt.subplots()
    x_axis = np.arange(len(agegroups))
    width = 0.2
    shift = -width
    for i in range(len(proportions)):
        ax.bar(x_axis+shift, proportions[i], width=width,
               label=f"{labels[i]}", color=colors[i])
        shift += width

    fig.supxlabel('Age groups')
    fig.supylabel('Proportion')
    plt.legend()
    plt.xticks(x_axis, agegroups)

    plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + "proportions_per_agegroup.png",
                    bbox_inches='tight', dpi=500)


def main():
    # Path where file with RKI case data is stored. 
    data_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/")
    
    # Path where plots will be stored. 
    save_dir =  os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/")

    # Mean stay times according to Covasim paper. 
    T_IH = 6.6
    T_HU = 1.5
    T_U = 15.230258

    start_dates = ["2020-06-01", "2020-10-01"]

    plot(data_dir, start_dates, T_IH, T_HU, T_U, save_dir)

if __name__ == "__main__":
    main()
