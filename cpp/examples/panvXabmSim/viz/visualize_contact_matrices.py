#!/usr/bin/env python3
"""
Contact Matrix Visualization Script
Visualizes age-stratified contact matrices for different countries and settings
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Age group labels
age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']

# ============================================================================
# GERMANY CONTACT MATRICES
# ============================================================================
germany_home = np.array([
    [0.4413, 0.0504, 1.2383, 0.8033, 0.0494, 0.0017],
    [0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013],
    [0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032],
    [0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014],
    [0.0087, 0.0394, 0.1417, 0.3834, 0.7064, 0.0447],
    [0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544]
])

germany_school = np.array([
    [1.1165, 0.2741, 0.2235, 0.1028, 0.0007, 0.0000],
    [0.1627, 1.9412, 0.2431, 0.1780, 0.0130, 0.0000],
    [0.0148, 0.1646, 1.1266, 0.0923, 0.0074, 0.0000],
    [0.0367, 0.1843, 0.3265, 0.0502, 0.0021, 0.0005],
    [0.0004, 0.0370, 0.0115, 0.0014, 0.0039, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]
])

germany_work = np.array([
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0127, 1.7570, 1.6050, 0.0133, 0.0000],
    [0.0000, 0.0020, 1.0311, 2.3166, 0.0098, 0.0000],
    [0.0000, 0.0002, 0.0194, 0.0325, 0.0003, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]
])

germany_other = np.array([
    [0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428],
    [0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148],
    [0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077],
    [0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178],
    [0.0371, 0.0393, 0.4171, 0.9666, 0.7495, 0.0257],
    [0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180]
])

# ============================================================================
# FRANCE CONTACT MATRICES
# ============================================================================
france_home = np.array([
    [0.6881, 0.6771, 1.2965, 0.9261, 0.0337, 0.0034],
    [0.2257, 1.6804, 0.7570, 1.4088, 0.0235, 0.0022],
    [0.2563, 0.3517, 1.4941, 0.7716, 0.0381, 0.0020],
    [0.2096, 0.6996, 0.8293, 1.2112, 0.0640, 0.0075],
    [0.1710, 0.4608, 0.4983, 0.6181, 0.8598, 0.0665],
    [0.1990, 0.6331, 0.5450, 1.1572, 0.4189, 0.2800]
])

france_school = np.array([
    [2.4233, 0.3734, 0.3248, 0.3342, 0.0032, 0.0000],
    [0.1972, 3.3900, 0.1830, 0.2888, 0.0054, 0.0001],
    [0.0317, 0.4252, 1.3594, 0.1850, 0.0051, 0.0002],
    [0.1477, 0.4910, 0.4562, 0.2218, 0.0061, 0.0001],
    [0.0207, 0.0356, 0.0564, 0.0603, 0.0277, 0.0035],
    [0.0000, 0.0212, 0.0289, 0.0000, 0.0000, 0.0000]
])

france_work = np.array([
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0074, 0.0222, 0.0254, 0.0000, 0.0000],
    [0.0000, 0.0235, 2.7163, 2.4826, 0.0050, 0.0000],
    [0.0000, 0.0249, 2.0297, 3.5135, 0.0052, 0.0000],
    [0.0000, 0.0005, 0.0133, 0.0262, 0.0001, 0.0000],
    [0.0000, 0.0000, 0.0001, 0.0001, 0.0000, 0.0000]
])

france_other = np.array([
    [0.6935, 0.4678, 0.9381, 1.1670, 0.3841, 0.0317],
    [0.2217, 2.2229, 0.8144, 1.0773, 0.2950, 0.0378],
    [0.1051, 0.4410, 2.7757, 1.2325, 0.1793, 0.0230],
    [0.0847, 0.2358, 1.1066, 1.8774, 0.4578, 0.0398],
    [0.0495, 0.1284, 0.7258, 1.5573, 1.0805, 0.1023],
    [0.0446, 0.1150, 0.4335, 1.0023, 0.8198, 0.1387]
])

# ============================================================================
# USA CONTACT MATRICES
# ============================================================================
usa_home = np.array([
    [0.6197, 0.7925, 1.1687, 0.9485, 0.0335, 0.0040],
    [0.2460, 1.7677, 0.8041, 1.3347, 0.0308, 0.0032],
    [0.2832, 0.4564, 1.5297, 0.7018, 0.0449, 0.0028],
    [0.2613, 0.8550, 0.9010, 1.1937, 0.0834, 0.0095],
    [0.2117, 0.5553, 0.5944, 0.7254, 0.6861, 0.0641],
    [0.1787, 0.5972, 0.5194, 0.9997, 0.2481, 0.2714]
])

usa_school = np.array([
    [1.1966, 0.2696, 0.2404, 0.2202, 0.0036, 0.0000],
    [0.1381, 3.9843, 0.2120, 0.2775, 0.0077, 0.0001],
    [0.0240, 0.5101, 1.8745, 0.1804, 0.0066, 0.0003],
    [0.0982, 0.4768, 0.4823, 0.1892, 0.0071, 0.0001],
    [0.0251, 0.0527, 0.0906, 0.0676, 0.0291, 0.0044],
    [0.0000, 0.0211, 0.0289, 0.0000, 0.0000, 0.0000]
])

usa_work = np.array([
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0402, 0.0502, 0.0548, 0.0000, 0.0000],
    [0.0000, 0.0526, 2.5209, 2.3056, 0.0387, 0.0000],
    [0.0000, 0.0664, 1.9258, 3.4223, 0.0442, 0.0000],
    [0.0000, 0.0105, 0.1102, 0.2350, 0.0047, 0.0000],
    [0.0000, 0.0000, 0.0001, 0.0001, 0.0000, 0.0000]
])

usa_other = np.array([
    [0.7819, 0.5386, 1.0918, 1.2430, 0.3503, 0.0309],
    [0.2522, 2.6569, 0.9886, 1.1725, 0.2698, 0.0380],
    [0.1245, 0.5258, 3.4633, 1.3746, 0.1660, 0.0232],
    [0.0913, 0.2606, 1.2291, 1.9120, 0.3962, 0.0368],
    [0.0469, 0.1202, 0.7026, 1.3731, 0.7760, 0.0768],
    [0.0315, 0.0835, 0.3139, 0.6754, 0.4473, 0.0847]
])


def plot_contact_matrix(ax, matrix, title, vmin=None, vmax=None, cmap='YlOrRd'):
    """Plot a single contact matrix as a heatmap."""
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(age_groups)))
    ax.set_yticks(np.arange(len(age_groups)))
    ax.set_xticklabels(age_groups, fontsize=20)
    ax.set_yticklabels(age_groups, fontsize=20)
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add text annotations
    for i in range(len(age_groups)):
        for j in range(len(age_groups)):
            value = matrix[i, j]
            # Choose text color based on background
            text_color = 'white' if value > (vmax * 0.6 if vmax else 1.5) else 'black'
            text = ax.text(j, i, f'{value:.2f}', ha="center", va="center",
                          color=text_color, fontsize=22)
    
    ax.set_title(title, fontsize=24, fontweight='bold', pad=10)
    ax.set_xlabel('Age Group of Contact', fontsize=22)
    ax.set_ylabel('Age Group', fontsize=22)
    
    return im


def create_country_comparison(output_file='contact_matrices_comparison.png'):
    """Create a comprehensive comparison of all contact matrices."""
    
    # Organize data
    countries = {
        'Germany': {
            'Home': germany_home,
            'School': germany_school,
            'Work': germany_work,
            'Other': germany_other
        },
        'France': {
            'Home': france_home,
            'School': france_school,
            'Work': france_work,
            'Other': france_other
        },
        'USA': {
            'Home': usa_home,
            'School': usa_school,
            'Work': usa_work,
            'Other': usa_other
        }
    }
    
    settings = ['Home', 'School', 'Work', 'Other']
    
    # Find global min/max for consistent color scaling within each setting
    setting_ranges = {}
    for setting in settings:
        all_values = []
        for country in countries.values():
            all_values.extend(country[setting].flatten())
        setting_ranges[setting] = (0, max(all_values))
    
    # Create figure
    fig, axes = plt.subplots(4, 3, figsize=(18, 20))
    fig.suptitle('Age-Stratified Contact Matrices by Country and Setting',
                 fontsize=24, fontweight='bold', y=0.995)
    
    # Plot matrices
    for col_idx, (country_name, matrices) in enumerate(countries.items()):
        for row_idx, setting in enumerate(settings):
            ax = axes[row_idx, col_idx]
            matrix = matrices[setting]
            vmin, vmax = setting_ranges[setting]
            
            title = f'{country_name} - {setting}'
            im = plot_contact_matrix(ax, matrix, title, vmin=vmin, vmax=vmax)
            
            # Add colorbar for the rightmost column only
            if col_idx == 2:
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label('Contact rate', rotation=270, labelpad=20, fontsize=18)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved comparison plot to {output_file}")
    plt.close()


def create_setting_comparison(output_file='contact_matrices_by_setting.png'):
    """Create comparison grouped by setting type."""
    
    settings_data = {
        'Home': {'Germany': germany_home, 'France': france_home, 'USA': usa_home},
        'School': {'Germany': germany_school, 'France': france_school, 'USA': usa_school},
        'Work': {'Germany': germany_work, 'France': france_work, 'USA': usa_work},
        'Other': {'Germany': germany_other, 'France': france_other, 'USA': usa_other}
    }
    
    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    fig.suptitle('Age-Stratified Contact Matrices by Setting Type',
                 fontsize=24, fontweight='bold', y=0.995)
    
    for col_idx, (setting, countries) in enumerate(settings_data.items()):
        # Find vmax for this setting
        all_values = []
        for matrix in countries.values():
            all_values.extend(matrix.flatten())
        vmax = max(all_values)
        
        for row_idx, (country, matrix) in enumerate(countries.items()):
            ax = axes[row_idx, col_idx]
            title = f'{country} - {setting}'
            im = plot_contact_matrix(ax, matrix, title, vmin=0, vmax=vmax)
            
            # Add colorbar for bottom row only
            if row_idx == 2:
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label('Contact rate', rotation=270, labelpad=20, fontsize=18)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved setting comparison plot to {output_file}")
    plt.close()


def create_difference_matrices(output_file='contact_matrices_differences.png'):
    """Create heatmaps showing differences between countries."""
    
    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    fig.suptitle('Contact Matrix Differences (France - Germany, USA - Germany)',
                 fontsize=24, fontweight='bold', y=0.995)
    
    settings = ['Home', 'School', 'Work', 'Other']
    comparisons = [
        ('France - Germany', france_home - germany_home, france_school - germany_school,
         france_work - germany_work, france_other - germany_other),
        ('USA - Germany', usa_home - germany_home, usa_school - germany_school,
         usa_work - germany_work, usa_other - germany_other)
    ]
    
    # Calculate global range for diverging colormap
    all_diffs = []
    for _, *matrices in comparisons:
        for matrix in matrices:
            all_diffs.extend(matrix.flatten())
    abs_max = max(abs(min(all_diffs)), abs(max(all_diffs)))
    
    for row_idx, (comparison_name, *matrices) in enumerate(comparisons):
        for col_idx, (setting, diff_matrix) in enumerate(zip(settings, matrices)):
            ax = axes[row_idx, col_idx]
            
            im = ax.imshow(diff_matrix, cmap='RdBu_r', aspect='auto',
                          vmin=-abs_max, vmax=abs_max)
            
            # Set ticks and labels
            ax.set_xticks(np.arange(len(age_groups)))
            ax.set_yticks(np.arange(len(age_groups)))
            ax.set_xticklabels(age_groups, fontsize=18)
            ax.set_yticklabels(age_groups, fontsize=18)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            
            # Add text annotations
            for i in range(len(age_groups)):
                for j in range(len(age_groups)):
                    value = diff_matrix[i, j]
                    text_color = 'white' if abs(value) > abs_max * 0.6 else 'black'
                    text = ax.text(j, i, f'{value:.2f}', ha="center", va="center",
                                  color=text_color, fontsize=18)
            
            ax.set_title(f'{comparison_name} - {setting}', fontsize=18, fontweight='bold')
            ax.set_xlabel('Age of contact', fontsize=18)
            ax.set_ylabel('Age Groupividual', fontsize=18)
            
            # Add colorbar for rightmost column
            if col_idx == 3:
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label('Difference', rotation=270, labelpad=20, fontsize=18)
    
    # Hide the third row (empty)
    for col_idx in range(4):
        axes[2, col_idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved difference plot to {output_file}")
    plt.close()


def create_total_contact_comparison(output_file='total_contacts_comparison.png'):
    """Create bar chart comparing total contact rates across countries and settings."""
    
    countries = ['Germany', 'France', 'USA']
    settings = ['Home', 'School', 'Work', 'Other']
    
    # Calculate total contacts (sum of all matrix elements)
    totals = {
        'Germany': [germany_home.sum(), germany_school.sum(), germany_work.sum(), germany_other.sum()],
        'France': [france_home.sum(), france_school.sum(), france_work.sum(), france_other.sum()],
        'USA': [usa_home.sum(), usa_school.sum(), usa_work.sum(), usa_other.sum()]
    }
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    x = np.arange(len(settings))
    width = 0.25
    
    for i, country in enumerate(countries):
        offset = width * (i - 1)
        bars = ax.bar(x + offset, totals[country], width, label=country, alpha=0.8)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1f}', ha='center', va='bottom', fontsize=14)
    
    ax.set_xlabel('Setting', fontsize=18, fontweight='bold')
    ax.set_ylabel('Total Contact Rate', fontsize=18, fontweight='bold')
    ax.set_title('Total Contact Rates by Country and Setting', fontsize=20, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(settings, fontsize=18)
    ax.legend(fontsize=18)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved total contacts comparison to {output_file}")
    plt.close()


def create_summed_matrices(output_file='contact_matrices_summed_all_locations.png'):
    """Create heatmaps showing sum of all location types for each country."""
    
    # Calculate summed matrices
    germany_total = germany_home + germany_school + germany_work + germany_other
    france_total = france_home + france_school + france_work + france_other
    usa_total = usa_home + usa_school + usa_work + usa_other
    
    # Find global vmax for consistent color scaling
    vmax = max(germany_total.max(), france_total.max(), usa_total.max())

    fig, axes = plt.subplots(3, 1, figsize=(10, 20))
    fig.suptitle('Total Contact Matrices (Sum of All Location Types)',
                 fontsize=28, fontweight='bold', y=0.995)

    countries_data = [
        ('Germany', germany_total),
        ('France', france_total),
        ('USA', usa_total)
    ]

    for idx, (country_name, matrix) in enumerate(countries_data):
        ax = axes[idx]
        im = plot_contact_matrix(ax, matrix, f'{country_name}',
                                vmin=0, vmax=vmax)

        # Add colorbar with larger font
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Contact rate', rotation=270, labelpad=22, fontsize=24)
        cbar.ax.tick_params(labelsize=24)

        # Add sum annotation
        # total_sum = matrix.sum()
        # ax.text(0.02, 0.98, f'Total contacts: {total_sum:.1f}',
        #         transform=ax.transAxes, fontsize=11, fontweight='bold',
        #         verticalalignment='top',
        #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout(h_pad=3.0)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved summed matrices to {output_file}")
    plt.close()


def main():
    """Generate all visualizations."""
    print("Generating contact matrix visualizations...")
    print()
    
    # Create all visualizations
    create_country_comparison('contact_matrices_comparison.png')
    create_setting_comparison('contact_matrices_by_setting.png')
    create_difference_matrices('contact_matrices_differences.png')
    create_summed_matrices('contact_matrices_summed_all_locations.png')
    create_total_contact_comparison('total_contacts_comparison.png')
    
    print()
    print("✓ All visualizations completed!")
    print()
    print("Generated files:")
    print("  - contact_matrices_comparison.png (4x3 grid by country)")
    print("  - contact_matrices_by_setting.png (3x4 grid by setting)")
    print("  - contact_matrices_differences.png (difference heatmaps)")
    print("  - contact_matrices_summed_all_locations.png (total matrices stacked)")
    print("  - total_contacts_comparison.png (bar chart comparison)")


if __name__ == "__main__":
    main()
