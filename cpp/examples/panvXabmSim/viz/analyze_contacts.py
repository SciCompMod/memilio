import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_infection_data(data_dir):
    """Load infection timeline data."""
    infection_file = os.path.join(data_dir, 'best_run_detailed_infection.csv')
    if not os.path.exists(infection_file):
        raise FileNotFoundError(f"Infection data not found: {infection_file}")

    print(f"  Loading infection data: {infection_file}")
    df = pd.read_csv(infection_file)
    print(f"  Loaded {len(df)} infection records")
    return df


def load_contact_data(data_dir):
    """Load contact data."""
    contact_file = os.path.join(data_dir, 'best_run_contact_data.csv')
    if not os.path.exists(contact_file):
        raise FileNotFoundError(f"Contact data not found: {contact_file}")

    print(f"  Loading contact data: {contact_file}")
    df = pd.read_csv(contact_file)
    print(f"  Loaded {len(df)} contact records")
    return df


def calculate_rolling_24h_average(data_df, window_size=24):
    """
    Calculate rolling 24-hour averages for contact metrics.

    Args:
        data_df: DataFrame with timestep data
        window_size: Rolling window size in timesteps (default 24 for 24 hours)

    Returns:
        DataFrame with additional rolling average columns
    """
    print(f"  Calculating {window_size}-hour rolling averages...")

    # Make a copy to avoid modifying original data
    rolling_data = data_df.copy()

    # Calculate rolling averages
    rolling_data['total_contacts_24h_avg'] = rolling_data['total_potential_contacts'].rolling(
        window=window_size, min_periods=1, center=True).mean()
    rolling_data['contacts_per_infectious_24h_avg'] = rolling_data['contacts_per_infectious'].rolling(
        window=window_size, min_periods=1, center=True).mean()
    rolling_data['infectious_count_24h_avg'] = rolling_data['infectious_count'].rolling(
        window=window_size, min_periods=1, center=True).mean()

    return rolling_data


def calculate_potential_contacts_per_timestep(infection_df, contact_df, rolling_window=24):
    """
    Calculate potential contacts that infectious persons have per timestep.

    For each timestep:
    1. Find all infectious people
    2. For each infectious person, find all people they could potentially contact
       (people in the same location at the same time)
    3. Sum up all potential contacts

    Args:
        infection_df: DataFrame with infection data
        contact_df: DataFrame with contact data
        rolling_window: Window size for rolling average calculation
    """
    print("  Calculating potential contacts per timestep...")

    # Get all timesteps from both datasets
    all_timesteps = sorted(set(infection_df['Timestep'].unique()) | set(
        contact_df['Timestep'].unique()))

    results = []

    for timestep in all_timesteps:
        # Get infectious people at this timestep
        infectious_people = set(
            infection_df[infection_df['Timestep'] == timestep]['Person_ID'])

        # Get all contacts at this timestep
        timestep_contacts = contact_df[contact_df['Timestep'] == timestep]

        if len(infectious_people) == 0 or len(timestep_contacts) == 0:
            results.append({
                'Timestep': timestep,
                'infectious_count': len(infectious_people),
                'total_potential_contacts': 0,
                'contacts_per_infectious': 0
            })
            continue

        # Calculate potential contacts for infectious people
        total_potential = 0

        # Group contacts by location to find potential contact opportunities
        location_groups = timestep_contacts.groupby(
            'Location_ID')['Person_ID'].apply(set).to_dict()

        # For each infectious person, count potential contacts
        for infectious_person in infectious_people:
            # Find locations where this infectious person is present
            infectious_locations = timestep_contacts[
                timestep_contacts['Person_ID'] == infectious_person
            ]['Location_ID'].unique()

            # Count potential contacts in each location
            for location in infectious_locations:
                if location in location_groups:
                    # All people in this location are potential contacts (excluding the infectious person)
                    potential_contacts = len(location_groups[location]) - 1
                    total_potential += potential_contacts

        results.append({
            'Timestep': timestep,
            'infectious_count': len(infectious_people),
            'total_potential_contacts': total_potential,
            'contacts_per_infectious': total_potential / len(infectious_people) if len(infectious_people) > 0 else 0
        })

    df_results = pd.DataFrame(results)
    print(f"  Calculated contacts for {len(df_results)} timesteps")

    # Add rolling averages with custom window size
    df_results_with_rolling = calculate_rolling_24h_average(
        df_results, window_size=rolling_window)

    return df_results_with_rolling


def create_comparison_visualization(memilio_data, panvadere_data, output_file, title="Contact Analysis Comparison", rolling_window=24):
    """Create comparison visualization between Memilio and Panvadere with rolling averages."""

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 15))
    fig.suptitle(title, fontsize=16, fontweight='bold')

    # Plot 1: Total Potential Contacts Per Timestep (Raw + Rolling Average)
    ax1.plot(memilio_data['Timestep'], memilio_data['total_potential_contacts'],
             'b-', linewidth=1, label='Memilio (Raw)', alpha=0.6)
    ax1.plot(memilio_data['Timestep'], memilio_data['total_contacts_24h_avg'],
             'b-', linewidth=3, label=f'Memilio ({rolling_window}h Avg)', alpha=0.9)
    ax1.plot(panvadere_data['Timestep'], panvadere_data['total_potential_contacts'],
             'r--', linewidth=1, label='Panvadere (Raw)', alpha=0.6)
    ax1.plot(panvadere_data['Timestep'], panvadere_data['total_contacts_24h_avg'],
             'r--', linewidth=3, label=f'Panvadere ({rolling_window}h Avg)', alpha=0.9)

    ax1.set_xlabel('Timestep (hours)')
    ax1.set_ylabel('Total Potential Contacts')
    ax1.set_title('Total Potential Contacts Available to Infectious Persons')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Add statistics for raw data
    memilio_max = memilio_data['total_potential_contacts'].max()
    panvadere_max = panvadere_data['total_potential_contacts'].max()
    memilio_mean = memilio_data['total_potential_contacts'].mean()
    panvadere_mean = panvadere_data['total_potential_contacts'].mean()

    stats_text = f'Raw Data:\nMemilio: Max={memilio_max:.0f}, Mean={memilio_mean:.1f}\nPanvadere: Max={panvadere_max:.0f}, Mean={panvadere_mean:.1f}'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Plot 2: Average Contacts Per Infectious Person (Raw + Rolling Average)
    ax2.plot(memilio_data['Timestep'], memilio_data['contacts_per_infectious'],
             'b-', linewidth=1, label='Memilio (Raw)', alpha=0.6)
    ax2.plot(memilio_data['Timestep'], memilio_data['contacts_per_infectious_24h_avg'],
             'b-', linewidth=3, label='Memilio (24h Avg)', alpha=0.9)
    ax2.plot(panvadere_data['Timestep'], panvadere_data['contacts_per_infectious'],
             'r--', linewidth=1, label='Panvadere (Raw)', alpha=0.6)
    ax2.plot(panvadere_data['Timestep'], panvadere_data['contacts_per_infectious_24h_avg'],
             'r--', linewidth=3, label='Panvadere (24h Avg)', alpha=0.9)

    ax2.set_xlabel('Timestep (hours)')
    ax2.set_ylabel('Average Contacts per Infectious Person')
    ax2.set_title('Average Potential Contacts per Infectious Person')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Add statistics for contacts per infectious
    memilio_max_per = memilio_data['contacts_per_infectious'].max()
    panvadere_max_per = panvadere_data['contacts_per_infectious'].max()
    memilio_mean_per = memilio_data['contacts_per_infectious'].mean()
    panvadere_mean_per = panvadere_data['contacts_per_infectious'].mean()

    per_stats_text = f'Raw Data:\nMemilio: Max={memilio_max_per:.1f}, Mean={memilio_mean_per:.1f}\nPanvadere: Max={panvadere_max_per:.1f}, Mean={panvadere_mean_per:.1f}'
    ax2.text(0.02, 0.98, per_stats_text, transform=ax2.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    # Plot 3: 24-Hour Rolling Average Comparison (Clean comparison)
    ax3.plot(memilio_data['Timestep'], memilio_data['total_contacts_24h_avg'],
             'b-', linewidth=3, label='Memilio (24h Avg)', alpha=0.9)
    ax3.plot(panvadere_data['Timestep'], panvadere_data['total_contacts_24h_avg'],
             'r-', linewidth=3, label='Panvadere (24h Avg)', alpha=0.9)

    ax3.set_xlabel('Timestep (hours)')
    ax3.set_ylabel('24-Hour Rolling Average Total Contacts')
    ax3.set_title('24-Hour Rolling Average: Total Potential Contacts')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Add rolling average statistics
    memilio_rolling_max = memilio_data['total_contacts_24h_avg'].max()
    panvadere_rolling_max = panvadere_data['total_contacts_24h_avg'].max()
    memilio_rolling_mean = memilio_data['total_contacts_24h_avg'].mean()
    panvadere_rolling_mean = panvadere_data['total_contacts_24h_avg'].mean()

    rolling_stats_text = f'24h Rolling Avg:\nMemilio: Max={memilio_rolling_max:.1f}, Mean={memilio_rolling_mean:.1f}\nPanvadere: Max={panvadere_rolling_max:.1f}, Mean={panvadere_rolling_mean:.1f}'
    ax3.text(0.02, 0.98, rolling_stats_text, transform=ax3.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    return {
        'memilio': {
            'max_total_contacts': memilio_max,
            'mean_total_contacts': memilio_mean,
            'max_contacts_per_infectious': memilio_max_per,
            'mean_contacts_per_infectious': memilio_mean_per,
            'total_infectious_events': memilio_data['infectious_count'].sum(),
            'max_total_contacts_24h_avg': memilio_rolling_max,
            'mean_total_contacts_24h_avg': memilio_rolling_mean
        },
        'panvadere': {
            'max_total_contacts': panvadere_max,
            'mean_total_contacts': panvadere_mean,
            'max_contacts_per_infectious': panvadere_max_per,
            'mean_contacts_per_infectious': panvadere_mean_per,
            'total_infectious_events': panvadere_data['infectious_count'].sum(),
            'max_total_contacts_24h_avg': panvadere_rolling_max,
            'mean_total_contacts_24h_avg': panvadere_rolling_mean
        }
    }


def print_comparison_summary(stats):
    """Print a summary comparison of the results."""
    memilio = stats['memilio']
    panvadere = stats['panvadere']

    print("\n" + "="*80)
    print("CONTACT ANALYSIS SUMMARY (with 24-hour Rolling Averages)")
    print("="*80)
    print(f"{'Metric':<40} {'Memilio':<15} {'Panvadere':<15} {'Difference':<15}")
    print("-"*80)

    metrics = [
        ('Max Total Contacts (Raw)', 'max_total_contacts', '{:.0f}'),
        ('Mean Total Contacts (Raw)', 'mean_total_contacts', '{:.1f}'),
        ('Max Total Contacts (24h Avg)',
         'max_total_contacts_24h_avg', '{:.1f}'),
        ('Mean Total Contacts (24h Avg)',
         'mean_total_contacts_24h_avg', '{:.1f}'),
        ('Max Contacts/Infectious (Raw)',
         'max_contacts_per_infectious', '{:.1f}'),
        ('Mean Contacts/Infectious (Raw)',
         'mean_contacts_per_infectious', '{:.1f}'),
        ('Total Infectious Events', 'total_infectious_events', '{:.0f}')
    ]

    for metric_name, key, fmt in metrics:
        m_val = memilio[key]
        p_val = panvadere[key]
        diff = m_val - p_val

        print(
            f"{metric_name:<40} {fmt.format(m_val):<15} {fmt.format(p_val):<15} {diff:+.1f}")

    print("="*80)

    # Additional rolling average insights
    print("\n24-HOUR ROLLING AVERAGE INSIGHTS:")
    print("-"*50)
    memilio_smoothness = memilio['mean_total_contacts_24h_avg'] / \
        memilio['mean_total_contacts'] * 100
    panvadere_smoothness = panvadere['mean_total_contacts_24h_avg'] / \
        panvadere['mean_total_contacts'] * 100
    print(f"Memilio 24h avg represents {memilio_smoothness:.1f}% of raw mean")
    print(
        f"Panvadere 24h avg represents {panvadere_smoothness:.1f}% of raw mean")

    rolling_diff = memilio['mean_total_contacts_24h_avg'] - \
        panvadere['mean_total_contacts_24h_avg']
    print(f"24h rolling average difference: {rolling_diff:+.1f} contacts")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze potential contacts for infectious persons with rolling averages')
    parser.add_argument('--memilio-dir', required=True,
                        help='Directory containing Memilio results')
    parser.add_argument('--panvadere-dir', required=True,
                        help='Directory containing Panvadere results')
    parser.add_argument(
        '--output-file', default='contact_analysis_comparison.png', help='Output PNG file')
    parser.add_argument(
        '--title', default='Contact Analysis Comparison', help='Plot title')
    parser.add_argument('--rolling-window', type=int, default=24,
                        help='Rolling average window size in hours (default: 24)')

    args = parser.parse_args()

    # Validate input directories
    for dir_path, name in [(args.memilio_dir, 'Memilio'), (args.panvadere_dir, 'Panvadere')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
            sys.exit(1)

    try:
        # Load and analyze Memilio data
        print("Loading Memilio data...")
        memilio_infection_df = load_infection_data(args.memilio_dir)
        memilio_contact_df = load_contact_data(args.memilio_dir)
        memilio_analysis = calculate_potential_contacts_per_timestep(
            memilio_infection_df, memilio_contact_df, args.rolling_window)

        # Load and analyze Panvadere data
        print("\nLoading Panvadere data...")
        panvadere_infection_df = load_infection_data(args.panvadere_dir)
        panvadere_contact_df = load_contact_data(args.panvadere_dir)
        panvadere_analysis = calculate_potential_contacts_per_timestep(
            panvadere_infection_df, panvadere_contact_df, args.rolling_window)

        # Create visualization
        print(f"\nCreating comparison visualization: {args.output_file}")
        stats = create_comparison_visualization(
            memilio_analysis, panvadere_analysis, args.output_file, args.title, args.rolling_window)

        # Print summary
        print_comparison_summary(stats)

        print(
            f"\n✓ Analysis complete! Visualization saved as: {args.output_file}")

    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
