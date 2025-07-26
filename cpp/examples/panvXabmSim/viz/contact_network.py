import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import seaborn as sns


def load_contact_data(contact_file):
    """Load contact data from CSV file"""
    return pd.read_csv(contact_file)


def load_infection_data(infection_file):
    """Load infection data from CSV file"""
    return pd.read_csv(infection_file)


def create_contact_network(contact_data, time_window=None):
    """
    Create a network graph from contact data

    Parameters:
    contact_data: DataFrame with columns ['person_1', 'person_2', 'time', 'location']
    time_window: tuple (start_time, end_time) to filter contacts, None for all time
    """
    # Filter by time window if specified
    if time_window:
        start_time, end_time = time_window
        contact_data = contact_data[
            (contact_data['time'] >= start_time) &
            (contact_data['time'] <= end_time)
        ]

    # Create network
    G = nx.Graph()

    # Add edges for each contact
    for _, row in contact_data.iterrows():
        person1, person2 = row['person_1'], row['person_2']
        time, location = row['time'], row['location']

        # Add edge with attributes
        if G.has_edge(person1, person2):
            # If edge exists, increment contact count and add time/location info
            G[person1][person2]['contact_count'] += 1
            G[person1][person2]['times'].append(time)
            G[person1][person2]['locations'].add(location)
        else:
            # Create new edge
            G.add_edge(person1, person2,
                       contact_count=1,
                       times=[time],
                       locations={location})

    return G


def calculate_infection_rates(infection_data, total_runs=100):
    """
    Calculate infection rates for each person across multiple runs

    Parameters:
    infection_data: DataFrame with columns ['person_id', 'run_id', 'infected', 'infection_time']
    total_runs: total number of simulation runs
    """
    # Count infections per person
    infection_counts = infection_data.groupby('person_id')['infected'].sum()

    # Calculate infection rates (percentage)
    infection_rates = (infection_counts / total_runs) * 100

    # Ensure all persons are included (even those never infected)
    all_persons = set(infection_data['person_id'].unique())
    for person in all_persons:
        if person not in infection_rates:
            infection_rates[person] = 0.0

    return infection_rates


def visualize_contact_network(G, infection_rates=None, figsize=(12, 10),
                              node_size=300, edge_width_multiplier=1,
                              layout='spring', save_path=None):
    """
    Visualize the contact network with infection rate coloring

    Parameters:
    G: NetworkX graph
    infection_rates: dict mapping person_id to infection rate (0-100)
    figsize: figure size tuple
    node_size: size of nodes
    edge_width_multiplier: multiplier for edge widths based on contact count
    layout: layout algorithm ('spring', 'circular', 'random', etc.)
    save_path: path to save the figure, None to display only
    """

    plt.figure(figsize=figsize)

    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(G, k=1, iterations=50)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'random':
        pos = nx.random_layout(G)
    else:
        pos = nx.spring_layout(G)

    # Prepare node colors based on infection rates
    if infection_rates is not None:
        node_colors = [infection_rates.get(node, 0) for node in G.nodes()]
        # Create colormap
        cmap = plt.cm.Reds
        vmin, vmax = 0, 100
    else:
        node_colors = 'lightblue'
        cmap = None
        vmin = vmax = None

    # Prepare edge widths based on contact frequency
    edge_widths = [G[u][v]['contact_count'] * edge_width_multiplier
                   for u, v in G.edges()]

    # Draw the network
    nx.draw_networkx_nodes(G, pos,
                           node_color=node_colors,
                           node_size=node_size,
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           alpha=0.8)

    nx.draw_networkx_edges(G, pos,
                           width=edge_widths,
                           alpha=0.6,
                           edge_color='gray')

    nx.draw_networkx_labels(G, pos,
                            font_size=8,
                            font_weight='bold')

    # Add colorbar if using infection rates
    if infection_rates is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=plt.gca(), shrink=0.8)
        cbar.set_label('Infection Rate (%)', rotation=270, labelpad=20)

    plt.title('Agent Contact Network\n(Node color = Infection rate, Edge width = Contact frequency)',
              fontsize=14, pad=20)
    plt.axis('off')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def print_network_stats(G, infection_rates=None):
    """Print basic network statistics"""
    print("=== Network Statistics ===")
    print(f"Number of nodes (agents): {G.number_of_nodes()}")
    print(f"Number of edges (contacts): {G.number_of_edges()}")
    print(f"Network density: {nx.density(G):.3f}")
    print(f"Average clustering coefficient: {nx.average_clustering(G):.3f}")

    # Contact frequency statistics
    contact_counts = [G[u][v]['contact_count'] for u, v in G.edges()]
    print(f"Average contacts per pair: {np.mean(contact_counts):.2f}")
    print(f"Max contacts per pair: {max(contact_counts)}")

    if infection_rates is not None:
        print(f"\n=== Infection Statistics ===")
        rates = list(infection_rates.values())
        print(f"Average infection rate: {np.mean(rates):.2f}%")
        print(f"Max infection rate: {max(rates):.2f}%")
        print(f"Agents never infected: {sum(1 for r in rates if r == 0)}")


def main():
    """Main function to run the analysis"""

    # Load data
    print("Loading data...")
    contact_data = load_contact_data('contact_data.csv')
    infection_data = load_infection_data('infection_data.csv')

    print(f"Loaded {len(contact_data)} contact records")
    print(f"Loaded {len(infection_data)} infection records")

    # Create contact network
    print("\nCreating contact network...")
    G = create_contact_network(contact_data)

    # Calculate infection rates
    print("Calculating infection rates...")
    infection_rates = calculate_infection_rates(infection_data, total_runs=100)

    # Print statistics
    print_network_stats(G, infection_rates)

    # Visualize network
    print("\nCreating visualization...")
    visualize_contact_network(G, infection_rates,
                              figsize=(14, 12),
                              node_size=500,
                              edge_width_multiplier=2,
                              layout='spring',
                              save_path='contact_network.png')

    # Create time-filtered visualization (example: first 24 hours)
    print("\nCreating time-filtered visualization (first 24 hours)...")
    G_filtered = create_contact_network(contact_data, time_window=(0, 24))
    visualize_contact_network(G_filtered, infection_rates,
                              figsize=(12, 10),
                              node_size=400,
                              layout='spring',
                              save_path='contact_network_24h.png')

    print("Analysis complete!")


if __name__ == "__main__":
    main()
