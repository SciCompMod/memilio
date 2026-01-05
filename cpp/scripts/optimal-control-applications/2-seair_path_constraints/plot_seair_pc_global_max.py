import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('result_seair_pc_global_max.csv')

# Convert Time to days, and Subjects columns to thousands for plotting
data['Time'] = data['Time']  # Time (in days)
data['Infected'] /= 1000      # Convert to thousands
data['Recovered'] /= 1000     # Convert to thousands
data['Dead'] /= 1000          # Convert to thousands
data['Exposed'] /= 1000       # Convert to thousands
data['Asymptomatic'] /= 1000 # Convert to thousands

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot with linear y-axis (ax1)
ax1.plot(data['Time'], data['Exposed'], label='Exposed', color='b')  # Blue for Exposed
ax1.plot(data['Time'], data['Asymptomatic'], label='Asymptomatic', color='orange')  # Orange for Asymptomatic
ax1.plot(data['Time'], data['Infected'], label='Infected', color='r')  # Red for Infected
ax1.plot(data['Time'], data['Recovered'], label='Recovered', color='g')  # Green for Recovered
ax1.plot(data['Time'], data['Dead'], label='Dead', color='k')  # Black for Dead

ax1.set_title('Epidemic Model over Time (Linear Scale)')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Subjects (Thousands)')
ax1.legend()
ax1.grid(True)

# Plot with logarithmic y-axis (ax2)
ax2.plot(data['Time'], data['Exposed'], label='Exposed', color='b')  # Blue for Exposed
ax2.plot(data['Time'], data['Asymptomatic'], label='Asymptomatic', color='orange')  # Orange for Asymptomatic
ax2.plot(data['Time'], data['Infected'], label='Infected', color='r')  # Red for Infected
ax2.plot(data['Time'], data['Recovered'], label='Recovered', color='g')  # Green for Recovered
ax2.plot(data['Time'], data['Dead'], label='Dead', color='k')  # Black for Dead

ax2.set_title('Epidemic Model over Time (Logarithmic Scale)')
ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Subjects (Thousands)')
ax2.set_yscale('log')
ax2.legend()
ax2.grid(True, which="both", ls="--")

# Adjust layout and show the plots
plt.tight_layout()
plt.savefig('result_seair_pc_global_max.png')  # Save the figure
plt.show()

# Load data from the CSV file
df = pd.read_csv('controls_seair_pc_global_max.csv')

# Plotting the step function with vertical lines at jumps
plt.figure(figsize=(10, 6))

# Plot each control parameter with step function
plt.step(df['ControlInterval'], df['SocialDistancing'], label='SocialDistancing', color='blue', where='post')
plt.step(df['ControlInterval'], df['Quarantined'], label='Quarantined', color='black', where='post')
plt.step(df['ControlInterval'], df['TestingRate'], label='TestingRate', color='red', where='post')

# Add vertical lines at jumps
for col in ['SocialDistancing', 'Quarantined', 'TestingRate']:
    for i in range(1, len(df[col])):
        if df[col].iloc[i] != df[col].iloc[i - 1]:  # Detect a jump
            plt.axvline(x=df['ControlInterval'].iloc[i], color='gray', linestyle='--', linewidth=1)

# Customize the plot
plt.xlabel('Time (days)')
plt.ylabel('Inputs')
plt.title('Control Parameters Over Time')
plt.legend()
plt.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
# Save the figure
plt.savefig('controls_seair_pc_global_max.png')
plt.show()
