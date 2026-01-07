import re
import matplotlib.pyplot as plt
import os

def parse_cpp_header_distributions(file_path):
    """
    Parses a C++ header file to extract age and household size distributions.
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None, None

    def extract_vector_data(vector_name):
        # Regex to find a std::vector<double> and capture its contents
        pattern = re.compile(r'const std::vector<double>\s+' + vector_name + r'\s*=\s*\{(.*?)\};', re.DOTALL)
        match = pattern.search(content)
        if not match:
            print(f"Warning: Vector '{vector_name}' not found in {file_path}")
            return None
        
        # Extract the content of the vector
        values_str = match.group(1)
        
        # Remove C-style comments before extracting numbers
        values_str = re.sub(r'//.*', '', values_str)
        values_str = re.sub(r'/\*.*?\*/', '', values_str, flags=re.DOTALL)

        # Extract all valid floating-point numbers from the cleaned string
        numbers_str = re.findall(r'([0-9\.]+)', values_str)
        
        numbers = []
        for num_s in numbers_str:
            try:
                numbers.append(float(num_s))
            except ValueError:
                print(f"Warning: Could not parse float from '{num_s}' in {file_path}")
                continue
        
        if not numbers:
            print(f"Warning: No data parsed for vector '{vector_name}' in {file_path}")
            return None
            
        return numbers

    age_distribution = extract_vector_data('AGE_DISTRIBUTION')
    household_distribution = extract_vector_data('HOUSEHOLD_SIZE_DISTRIBUTION')

    return age_distribution, household_distribution

def plot_demographics(countries_data, output_path):
    """
    Creates and saves a figure with pie charts for age and household distributions for multiple countries.
    """
    num_countries = len(countries_data)
    if num_countries == 0:
        return

    # Create a horizontal layout: 2 rows, num_countries columns
    fig, axes = plt.subplots(2, num_countries, figsize=(6 * num_countries, 12), squeeze=False)

    fig.suptitle('Demographic Distributions by Country', fontsize=18, fontweight='bold')

    age_labels_template = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    household_labels_template = ['1 Person', '2 People', '3 People', '4 People', '5+ People']

    age_colors = ['#FFA07A', '#20B2AA', '#DA70D6', '#32CD32', "#989606", "#FE8515"]
    household_colors = plt.cm.Set2.colors

    for i, (country_name, data) in enumerate(countries_data.items()):
        age_data, household_data = data
        
        # Top row for age distributions
        ax_age = axes[0, i]
        # Bottom row for household distributions
        ax_hh = axes[1, i]

        # Age Distribution Pie Chart
        if age_data:
            # Ensure labels match data size
            current_age_labels = age_labels_template[:len(age_data)]
            ax_age.pie(age_data, labels=current_age_labels, autopct='%1.1f%%', startangle=90, counterclock=False,
<<<<<<< Updated upstream
                       colors=age_colors, pctdistance=0.85)
            ax_age.set_title(f'{country_name} - Age Distribution', fontweight='bold')
=======
                       colors=age_colors)
            if i == num_countries // 2:
                ax_age.set_title(f'Age Distribution\n{country_name}', fontweight='bold', fontsize=14, pad=10)
            else:
                ax_age.set_title(f'{country_name}', fontweight='bold', pad=10)
>>>>>>> Stashed changes
        else:
            ax_age.text(0.5, 0.5, 'Age data not found', ha='center', va='center')
            if i == num_countries // 2:
                ax_age.set_title(f'Age Distribution\n{country_name}', fontweight='bold', fontsize=14, pad=10)
            else:
                ax_age.set_title(f'{country_name}', fontweight='bold', pad=10)
        ax_age.axis('equal')

        # Household Distribution Pie Chart
        if household_data:
            # Ensure labels match data size
            current_household_labels = household_labels_template[:len(household_data)]
            ax_hh.pie(household_data, labels=current_household_labels, autopct='%1.1f%%', startangle=90, counterclock=False,
                      colors=household_colors)
            if i == num_countries // 2:
                ax_hh.set_title(f'Household Size Distribution\n{country_name}', fontweight='bold', fontsize=14, pad=10)
            else:
                ax_hh.set_title(f'{country_name}', fontweight='bold', pad=10)
        else:
            ax_hh.text(0.5, 0.5, 'Household data not found', ha='center', va='center')
            if i == num_countries // 2:
                ax_hh.set_title(f'Household Size Distribution\n{country_name}', fontweight='bold', fontsize=14, pad=10)
            else:
                ax_hh.set_title(f'{country_name}', fontweight='bold', pad=10)
        ax_hh.axis('equal')

    plt.tight_layout(rect=[0, 0, 1, 0.96], h_pad=4.0)
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    plt.savefig(output_path, dpi=300)
    plt.close()

def main():
    # Assuming the script is run from the 'cpp' directory or similar
    base_path = os.path.join(os.path.dirname(__file__), '..', 'include')
    
    countries = {
        "Germany": os.path.join(base_path, 'city_parameters.h'),
        "France": os.path.join(base_path, 'city_parameters_france.h'),
        "USA": os.path.join(base_path, 'city_parameters_us.h')
    }

    countries_data = {}
    for name, path in countries.items():
        age_dist, hh_dist = parse_cpp_header_distributions(path)
        if age_dist and hh_dist:
            countries_data[name] = (age_dist, hh_dist)

    if not countries_data:
        print("No data could be parsed. Exiting.")
        return

    output_file = os.path.join(os.path.dirname(__file__), '..', 'results', 'demographics_comparison.png')
    
    plot_demographics(countries_data, output_file)

if __name__ == "__main__":
    main()
