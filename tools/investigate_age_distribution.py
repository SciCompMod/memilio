import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def get_df_daily():
    # Read file.
    datafile = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_age_all_dates.json")
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


def get_population_per_agegroup():
    # Population from Table 12411-04-02-4-B from regionalstatistik.de, data from 31.12.2020.
    population_per_agegroup = np.array(
        [3969138, 7508662, 18921292, 28666166, 18153339, 5936434])
    total_population = population_per_agegroup.sum()
    population_per_agegroup = population_per_agegroup/total_population

    return population_per_agegroup


def get_relevant_confirmed_cases(start_date, T_IH, T_HU):
    # Get dataframe with daily confirmed cases.
    df = get_df_daily()

    # Considered age groups.
    agegroups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
    T_U = 16.49
    # Extract relevant dates to be considered.
    df_date = df[(df["Date"] >= pd.Timestamp(start_date)-pd.Timedelta(days=T_IH+T_HU+T_U))
                 & (df["Date"] <= pd.Timestamp(start_date)+pd.Timedelta(days=45-(T_IH+T_HU+T_U)))]  # T_IH+T_HU

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


def plot(start_dates, T_IH, T_HU):

    agegroups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']

    # Get proportions per age group for June and Ocotber scenario and corresponding share of total population.
    proportions_june = get_relevant_confirmed_cases(start_dates[0], T_IH, T_HU)
    proportions_october = get_relevant_confirmed_cases(
        start_dates[1], T_IH, T_HU)
    population_per_agegroup = get_population_per_agegroup()

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

    plt.savefig(f"plots/proportions_per_agegroup.png", dpi=500)


def covasim_to_rki_agegroups(mu_covasim):
    mu_rki = np.zeros(6)
    # Covasim age groups are split every 10 years.
    # A00-A04
    mu_rki[0] = (0.5*mu_covasim[0])/0.5
    # A05-A14
    mu_rki[1] = (0.5*mu_covasim[0] + 0.5*mu_covasim[1])/1.
    # A15-A34
    mu_rki[2] = (0.5*mu_covasim[1] + mu_covasim[2] + 0.5 * mu_covasim[3])/2.
    # A35-A59
    mu_rki[3] = (0.5 * mu_covasim[3] + mu_covasim[4] + mu_covasim[5])/2.5
    # A60-A79
    mu_rki[4] = (mu_covasim[6] + mu_covasim[7])/2.
    # A80+
    mu_rki[5] = (mu_covasim[8] + mu_covasim[9])/2.

    return mu_rki


def compute_covasim_probs_per_rki_agegroup():
    # Vectors with probabilities are directly from Covasim Paper. These are adjusted accordingly to our probabilities by division.
    # In Covasim, there are probabilities given from C->I, C->H, C->U and C->D.
    # We calculate our needed transition probabilities by dividing accordingly.

    mu_CI_covasim = np.array(
        [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.90])

    mu_CH_covasim = np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430,
                              0.07650, 0.13280, 0.20655, 0.24570, 0.24570])
    mu_IH_covasim = mu_CH_covasim/mu_CI_covasim

    mu_CU_covasim = np.array([0.00003, 0.00008, 0.00036, 0.00104, 0.00216,
                              0.00933, 0.03639, 0.08923, 0.17420, 0.17420])
    mu_HU_covasim = mu_CU_covasim/(mu_CH_covasim)

    mu_CD_covasim = np.array([0.00002, 0.00002, 0.00010, 0.00032, 0.00098,
                              0.00265, 0.00766, 0.02439, 0.08292, 0.16190])
    mu_UD_covasim = mu_CD_covasim/mu_CU_covasim

    # # Compute average by population just to test.
    # agegroups_covasim = np.array([7752706.0, 7581868,  9483430, 10871964, 10070748,
    #                               13304542,  10717241, 7436098, 5092743, 843691])
    # total_pop = agegroups_covasim.sum()

    # mu_CI_average = 0
    # mu_IH_average = 0
    # mu_HU_average = 0
    # mu_UD_average = 0
    # for i in range(len(agegroups_covasim)):
    #     mu_CI_average += mu_CI_covasim[i]*agegroups_covasim[i]/total_pop
    #     mu_IH_average += mu_IH_covasim[i]*agegroups_covasim[i]/total_pop
    #     mu_HU_average += mu_HU_covasim[i]*agegroups_covasim[i]/total_pop
    #     mu_UD_average += mu_UD_covasim[i]*agegroups_covasim[i]/total_pop

    # print("Covasim probs by pop: ", mu_CI_average,
    #       mu_IH_average, mu_HU_average, mu_UD_average)

    # Convert from 10 agegroups from Covasim Paper to 6 age groups according to RKI data.
    mu_CI_rki = covasim_to_rki_agegroups(mu_CI_covasim)
    mu_IH_rki = covasim_to_rki_agegroups(mu_IH_covasim)
    mu_HU_rki = covasim_to_rki_agegroups(mu_HU_covasim)
    mu_UD_rki = covasim_to_rki_agegroups(mu_UD_covasim)

    return mu_CI_rki, mu_IH_rki, mu_HU_rki, mu_UD_rki


def compute_adapted_mu(start_date, T_IH, T_HU):
    mu_CI_age, mu_IH_age, mu_HU_age, mu_UD_age = compute_covasim_probs_per_rki_agegroup()

    population_share = get_relevant_confirmed_cases(start_date, T_IH, T_HU)

    mu_CI = 0
    mu_IH = 0
    mu_HU = 0
    mu_UD = 0
    for i in range(len(population_share)):
        mu_CI += mu_CI_age[i] * population_share[i]
        mu_IH += mu_IH_age[i] * population_share[i]
        mu_HU += mu_HU_age[i] * population_share[i]
        mu_UD += mu_UD_age[i] * population_share[i]

    return mu_CI, mu_IH, mu_HU, mu_UD


def compute_mu_by_population():
    mu_CI_age, mu_IH_age, mu_HU_age, mu_UD_age = compute_covasim_probs_per_rki_agegroup()
    population_per_agegroup = get_population_per_agegroup()

    mu_CI = 0
    mu_IH = 0
    mu_HU = 0
    mu_UD = 0
    for i in range(len(population_per_agegroup)):
        mu_CI += mu_CI_age[i] * population_per_agegroup[i]
        mu_IH += mu_IH_age[i] * population_per_agegroup[i]
        mu_HU += mu_HU_age[i] * population_per_agegroup[i]
        mu_UD += mu_UD_age[i] * population_per_agegroup[i]

    return mu_CI, mu_IH, mu_HU, mu_UD


def mu_assessment_by_cases(start_date, T_IH, T_HU):
    population_share = get_relevant_confirmed_cases(
        start_date, T_IH, T_HU)

    mu_CI_assessment = np.array([0.75, 0.75, 0.8, 0.8, 0.8, 0.8])
    mu_IH_assessment = np.array([0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225])
    mu_HU_assessment = np.array([0.075, 0.075, 0.15, 0.15, 0.3, 0.4])
    mu_UD_assessment = np.array([0.05, 0.05, 0.14, 0.14, 0.4, 0.6])

    mu_CI = 0
    mu_IH = 0
    mu_HU = 0
    mu_UD = 0
    for i in range(len(population_share)):
        mu_CI += mu_CI_assessment[i] * population_share[i]
        mu_IH += mu_IH_assessment[i] * population_share[i]
        mu_HU += mu_HU_assessment[i] * population_share[i]
        mu_UD += mu_UD_assessment[i] * population_share[i]

    return mu_CI, mu_IH, mu_HU, mu_UD


def main():

    df_daily = get_df_daily()

    T_IH = 6.6
    T_HU = 1.5
    T_UD = 10.7

    # # from Assessment Paper
    # mu_IH = 0.0786429
    # mu_HU = 0.173176

    # # from Covasim
    # mu_IH = 0.104907
    # mu_HU = 0.369201

    start_dates = ["2020-06-01", "2020-10-01"]

    plot(start_dates, T_IH, T_HU)

    mu_CI, mu_IH, mu_HU, mu_UD = compute_adapted_mu(start_dates[0], T_IH, T_HU)
    print(f"mu {start_dates[0]}: {mu_CI}, {mu_IH}, {mu_HU}, {mu_UD}")

    mu_CI, mu_IH, mu_HU, mu_UD = compute_adapted_mu(start_dates[1], T_IH, T_HU)
    print(f"mu {start_dates[1]}: {mu_CI}, {mu_IH}, {mu_HU}, {mu_UD}")

    mu_CI, mu_IH, mu_HU, mu_UD = compute_mu_by_population()
    print(f"mu by population: {mu_CI}, {mu_IH}, {mu_HU}, {mu_UD}")

    # mu_CI, mu_IH, mu_HU, mu_UD = mu_assessment_by_cases(start_dates[0], T_IH, T_HU)
    # print(f"Assessment mu by cases June: {mu_CI}, {mu_IH}, {mu_HU}, {mu_UD}")

    # mu_CI, mu_IH, mu_HU, mu_UD = mu_assessment_by_cases(start_dates[1], T_IH, T_HU)
    # print(f"Assessment mu by cases October: {mu_CI}, {mu_IH}, {mu_HU}, {mu_UD}")


if __name__ == "__main__":
    main()
