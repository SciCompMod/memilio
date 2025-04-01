import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Set global plot parameters for publication quality
plt.rcParams['font.size'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 20
plt.rcParams['figure.figsize'] = (10, 8)

# Define paths
brunswick_data_path = '/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/data/mobility/braunschweig_result_ffa8.csv'
munich_data_path = '/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/data/mobility/muenchen_microscopic_result.csv'

# Define output path for figures
output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mobility_comparison_figures')
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Dictionary for leisure activities and vehicle choice
dict_leisure = {1: 'work', 2: 'education', 3: 'shopping', 4: 'free time',
                5: 'private matters', 6: 'others', 7: 'home', 0: 'not specified'}
dict_vehicle = {1: 'bicyle', 2: 'car_driver',
                3: 'car_codriver', 4: 'public transport', 5: 'walk'}

# Define age groups as in C++ code
age_group_0_to_4 = 0
age_group_5_to_14 = 1
age_group_15_to_34 = 2
age_group_35_to_59 = 3
age_group_60_to_79 = 4
age_group_80_plus = 5

# Map age cohorts to readable labels
age_group_labels = {
    age_group_0_to_4: '0-4',
    age_group_5_to_14: '5-14',
    age_group_15_to_34: '15-34',
    age_group_35_to_59: '35-59',
    age_group_60_to_79: '60-79',
    age_group_80_plus: '80+'
}

# Load data
print("Loading Brunswick data...")
bs = pd.read_csv(brunswick_data_path, header=None, skiprows=1)
bs.rename(
    columns={0: 'personID', 1: 'startZone', 2: 'destZone', 3: 'loc_id_start', 4: 'loc_id_end',
             5: 'countyStart', 6: 'countyEnd', 7: 'hhID', 8: 'tripDistance', 9: 'startTime', 
             10: 'travelTime', 11: 'loCs', 12: 'laCs', 13: 'loCe', 14: 'laCe', 15: 'vehicleChoice', 
             16: 'ActivityBefore', 17: 'ActivityAfter', 18: 'age', 19: 'home_in_bs', 
             20: 'location_type'},
    inplace=True)

print("Loading Munich data...")
mu = pd.read_csv(munich_data_path, header=None, skiprows=1)
mu.rename(
    columns={0: 'personID', 1: 'startZone', 2: 'destZone', 3: 'loc_id_start', 4: 'loc_id_end',
             5: 'countyStart', 6: 'countyEnd', 7: 'hhID', 8: 'tripDistance', 9: 'startTime', 
             10: 'travelTime', 11: 'loCs', 12: 'laCs', 13: 'loCe', 14: 'laCe', 15: 'vehicleChoice', 
             16: 'ActivityBefore', 17: 'ActivityAfter', 18: 'age', 19: 'home_in_mu', 
             20: 'location_type'},
    inplace=True)

# Replace invalid ages with NaN
bs.loc[bs['age'] == 999, 'age'] = np.nan
mu.loc[mu['age'] == 999, 'age'] = np.nan

print(f"Brunswick data loaded: {len(bs)} trips")
print(f"Munich data loaded: {len(mu)} trips")

# Function to create side-by-side plots
def create_comparison_plot(title, filename, figsize=(16, 7)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    fig.suptitle(title, fontsize=16)
    return fig, ax1, ax2

# Map ages to our defined age groups
def map_to_age_group(age):
    if pd.isna(age):
        return np.nan
    elif age <= 4:
        return age_group_0_to_4
    elif age <= 14:
        return age_group_5_to_14
    elif age <= 34:
        return age_group_15_to_34
    elif age <= 59:
        return age_group_35_to_59
    elif age <= 79:
        return age_group_60_to_79
    else:
        return age_group_80_plus

# Apply age group mapping to both datasets
bs['AgeGroup'] = bs['age'].apply(map_to_age_group)
mu['AgeGroup'] = mu['age'].apply(map_to_age_group)

# 1. Age Distribution Comparison
print("Creating age distribution comparison...")
fig, ax_bs, ax_mu = create_comparison_plot("Age Distribution Comparison", "age_distribution.png")

# Brunswick
bs_persons_ages = bs[['personID', 'AgeGroup']].drop_duplicates()
bs_persons_ages = bs_persons_ages[bs_persons_ages['AgeGroup'].notna()]
bs_persons_ages_cohorts = bs_persons_ages.groupby(['AgeGroup']).size().reset_index(name='counts')
bs_persons_ages_cohorts['AgeGroupLabel'] = bs_persons_ages_cohorts['AgeGroup'].map(age_group_labels)

bs_persons_ages_cohorts.plot.bar(x='AgeGroupLabel', y='counts', ax=ax_bs, color='tab:blue')
ax_bs.set_title("Brunswick")
ax_bs.set_xlabel("Age groups")
ax_bs.set_ylabel("Number of persons")
ax_bs.grid(axis='y', linestyle='--', alpha=0.7)

# Munich
mu_persons_ages = mu[['personID', 'AgeGroup']].drop_duplicates()
mu_persons_ages = mu_persons_ages[mu_persons_ages['AgeGroup'].notna()]
mu_persons_ages_cohorts = mu_persons_ages.groupby(['AgeGroup']).size().reset_index(name='counts')
mu_persons_ages_cohorts['AgeGroupLabel'] = mu_persons_ages_cohorts['AgeGroup'].map(age_group_labels)

mu_persons_ages_cohorts.plot.bar(x='AgeGroupLabel', y='counts', ax=ax_mu, color='tab:orange')
ax_mu.set_title("Munich")
ax_mu.set_xlabel("Age groups")
ax_mu.set_ylabel("Number of persons")
ax_mu.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'age_distribution_comparison.png'), dpi=300)

# 2. Heatmap Leisure Comparison (excluding "not specified")
print("Creating leisure heatmap comparison...")
fig, ax_bs, ax_mu = create_comparison_plot("Leisure Activities Transition Comparison", "heatmap_leisure.png")

# Function to create filtered heatmap excluding "not specified" (key 0)
def create_filtered_heatmap(data, ax, title):
    # Filter out "not specified" category (key 0)
    filtered_data = data[data['ActivityBefore'] != 0]
    filtered_data = filtered_data[filtered_data['ActivityAfter'] != 0]
    
    # Create crosstab
    matrix_leisure = pd.crosstab(
        filtered_data['ActivityBefore'],
        filtered_data['ActivityAfter'],
        normalize='index')
    
    # Get filtered dictionary without "not specified"
    filtered_dict_leisure = {k: v for k, v in dict_leisure.items() if k != 0}
    
    # Create heatmap
    sns.heatmap(matrix_leisure, vmin=0, vmax=1, 
                xticklabels=[filtered_dict_leisure.get(i) for i in matrix_leisure.columns],
                yticklabels=[filtered_dict_leisure.get(i) for i in matrix_leisure.index],
                cmap="Blues", ax=ax)
    ax.set_title(title)
    ax.set_xlabel('Activity after trip')
    ax.set_ylabel('Activity before trip')

# Brunswick
create_filtered_heatmap(bs, ax_bs, "Brunswick")

# Munich
create_filtered_heatmap(mu, ax_mu, "Munich")

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'heatmap_leisure_comparison.png'), dpi=300)

# 3. Number of trips to different location types (with dict_leisure mapping)
print("Creating location types comparison...")
fig, ax_bs, ax_mu = create_comparison_plot("Number of Trips to Different Location Types", "location_types.png")

# Function to map activity codes to names
def map_activity_to_name(series):
    return series.map(lambda x: dict_leisure.get(x, f'Unknown ({x})'))

# Brunswick
bs_location_types = bs.groupby(['ActivityAfter']).size()
bs_location_types.index = map_activity_to_name(bs_location_types.index)
bs_location_types.plot(kind='bar', ax=ax_bs, color='tab:blue')
ax_bs.set_title("Brunswick")
ax_bs.set_xlabel("Location type")
ax_bs.set_ylabel("Number of trips")
ax_bs.grid(axis='y', linestyle='--', alpha=0.7)

# Munich
mu_location_types = mu.groupby(['ActivityAfter']).size()
mu_location_types.index = map_activity_to_name(mu_location_types.index)
mu_location_types.plot(kind='bar', ax=ax_mu, color='tab:orange')
ax_mu.set_title("Munich")
ax_mu.set_xlabel("Location type")
ax_mu.set_ylabel("Number of trips")
ax_mu.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'location_types_comparison.png'), dpi=300)

# 4. Number of trips per person
print("Creating trips per person comparison...")
fig, ax_bs, ax_mu = create_comparison_plot("Number of Trips per Person", "trips_per_person.png")

# Brunswick
bs_persons = bs.groupby(['personID']).size().reset_index(name='counts')
bs_trips_per_person = bs_persons['counts'].value_counts().sort_index()
bs_trips_per_person.plot(kind='bar', ax=ax_bs, color='tab:blue')
ax_bs.set_title("Brunswick")
ax_bs.set_xlabel("Number of trips per day")
ax_bs.set_ylabel("Number of persons")
ax_bs.grid(axis='y', linestyle='--', alpha=0.7)

# Munich
mu_persons = mu.groupby(['personID']).size().reset_index(name='counts')
mu_trips_per_person = mu_persons['counts'].value_counts().sort_index()
mu_trips_per_person.plot(kind='bar', ax=ax_mu, color='tab:orange')
ax_mu.set_title("Munich")
ax_mu.set_xlabel("Number of trips per day")
ax_mu.set_ylabel("Number of persons")
ax_mu.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'trips_per_person_comparison.png'), dpi=300)

# 5. Trip purpose comparison by age group - using pie charts
print("Creating trip purpose comparison by age group...")

# Define age group labels for plot titles
age_group_title_labels = {
    age_group_0_to_4: '0-4 years',
    age_group_5_to_14: '5-14 years',
    age_group_15_to_34: '15-34 years',
    age_group_35_to_59: '35-59 years',
    age_group_60_to_79: '60-79 years',
    age_group_80_plus: '80+ years'
}

# Create a subplot for each age group (2 cities × 6 age groups)
fig, axes = plt.subplots(2, 6, figsize=(24, 10))
fig.suptitle("Trip Purpose Distribution by Age Group", fontsize=22)

# Custom colors
colors = plt.cm.tab10(np.arange(len(dict_leisure)))

# Function to create pie chart with improved percentage display
def create_better_pie(ax, data_dict, title):
    if not data_dict or sum(data_dict.values()) == 0:
        ax.set_title(title)
        ax.axis('off')  # Hide axis if no data
        return
    
    labels = list(data_dict.keys())
    sizes = list(data_dict.values())
    
    # Only show percentages for slices that are large enough (>= 5%)
    total = sum(sizes)
    autopct_func = lambda pct: f'{pct:.1f}%' if pct >= 5 else ''
    
    wedges, texts, autotexts = ax.pie(sizes, labels=None, 
                                      autopct=autopct_func, 
                                      startangle=90, 
                                      colors=colors,
                                      wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    
    # Increase font size for percentage labels and move them inward
    plt.setp(autotexts, size=12, weight='bold')
    for autotext in autotexts:
        # Move text closer to center to avoid overlap
        x, y = autotext.get_position()
        autotext.set_position((1.1*x, 1.1*y))
    
    ax.set_title(title, fontsize=16)
    ax.axis('equal')  # Equal aspect ratio ensures pie is drawn as a circle

# Process each age group
for i, age_group in enumerate([age_group_0_to_4, age_group_5_to_14, age_group_15_to_34, 
                              age_group_35_to_59, age_group_60_to_79, age_group_80_plus]):
    
    # Brunswick - top row
    bs_age_group = bs[bs['AgeGroup'] == age_group]
    bs_trip_purpose = bs_age_group.groupby(['ActivityAfter']).size()
    
    # Map activity codes to names and create dictionary for the pie chart
    bs_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in bs_trip_purpose.items() if idx in dict_leisure}
    
    # Create improved pie chart
    create_better_pie(axes[0, i], bs_trip_purpose_dict, f"Brunswick: {age_group_title_labels[age_group]}")
    
    # Munich - bottom row
    mu_age_group = mu[mu['AgeGroup'] == age_group]
    mu_trip_purpose = mu_age_group.groupby(['ActivityAfter']).size()
    
    mu_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in mu_trip_purpose.items() if idx in dict_leisure}
    
    create_better_pie(axes[1, i], mu_trip_purpose_dict, f"Munich: {age_group_title_labels[age_group]}")

# Create a common legend for all pie charts
all_activities = list(dict_leisure.values())
legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i % len(colors)]) for i in range(len(all_activities))]
fig.legend(legend_handles, all_activities, loc='lower center', ncol=len(all_activities)//2, bbox_to_anchor=(0.5, 0), fontsize=20)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15, top=0.9)
plt.savefig(os.path.join(output_path, 'trip_purpose_by_age_comparison.png'), dpi=300)

# Apply the same improvements to the overall pie charts
print("Creating overall trip purpose comparison...")
fig, (ax_bs, ax_mu) = plt.subplots(1, 2, figsize=(16, 8))
fig.suptitle("Trip Purpose Distribution", fontsize=18)

# Convert numeric activity codes to human-readable labels
# Brunswick
bs_trip_purpose = bs.groupby(['ActivityAfter']).size()
bs_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in bs_trip_purpose.items() if idx in dict_leisure}
create_better_pie(ax_bs, bs_trip_purpose_dict, "Brunswick")

# Munich
mu_trip_purpose = mu.groupby(['ActivityAfter']).size()
mu_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in mu_trip_purpose.items() if idx in dict_leisure}
create_better_pie(ax_mu, mu_trip_purpose_dict, "Munich")

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'trip_purpose_comparison.png'), dpi=300)

# 6. Overall Trip Purpose Distribution with bigger captions
print("Creating overall trip purpose comparison...")
fig, (ax_bs, ax_mu) = plt.subplots(1, 2, figsize=(16, 8))
fig.suptitle("Trip Purpose Distribution", fontsize=20)

# Convert numeric activity codes to human-readable labels
# Brunswick
bs_trip_purpose = bs.groupby(['ActivityAfter']).size()
bs_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in bs_trip_purpose.items() if idx in dict_leisure}
labels_bs = list(bs_trip_purpose_dict.keys())
sizes_bs = list(bs_trip_purpose_dict.values())

# Munich
mu_trip_purpose = mu.groupby(['ActivityAfter']).size()
mu_trip_purpose_dict = {dict_leisure[idx]: count for idx, count in mu_trip_purpose.items() if idx in dict_leisure}
labels_mu = list(mu_trip_purpose_dict.keys())
sizes_mu = list(mu_trip_purpose_dict.values())

# Function for enhanced pie charts with bigger labels and percentages
def create_big_caption_pie(ax, sizes, labels, title):
    # Create the pie chart
    wedges, texts, autotexts = ax.pie(
        sizes, 
        labels=labels, 
        autopct='%1.1f%%', 
        startangle=90, 
        colors=colors,
        textprops={'fontsize': 16},  # Bigger label font size
        wedgeprops={'linewidth': 1, 'edgecolor': 'white'}
    )
    
    # Make percentage labels bigger and bolder
    plt.setp(autotexts, size=16, weight='bold')
    
    # Adjust position of labels for better readability
    for text in texts:
        x, y = text.get_position()
        text.set_position((1.0*x, 1.0*y))  # Move labels further out
    
    ax.set_title(title, fontsize=18)
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle

# Plot enhanced pie charts
create_big_caption_pie(ax_bs, sizes_bs, labels_bs, "Brunswick")
create_big_caption_pie(ax_mu, sizes_mu, labels_mu, "Munich")

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(output_path, 'trip_purpose_comparison.png'), dpi=300)

print(f"All comparison plots saved to: {output_path}")