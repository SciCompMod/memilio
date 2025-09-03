import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import gamma
import seaborn as sns

# Set style for better plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# Set larger font sizes globally
plt.rcParams.update({
    'font.size': 16,           # Default font size
    'axes.titlesize': 18,      # Subplot titles
    'axes.labelsize': 16,      # Axis labels
    'xtick.labelsize': 14,     # X-axis tick labels
    'ytick.labelsize': 14,     # Y-axis tick labels
    'legend.fontsize': 13,     # Legend text
    'figure.titlesize': 22     # Main figure title
})

# Parameters from the table - Age group specific values
# Age groups: 0-4, 5-14, 15-34, 35-59, 60-79, 80+
age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']

# Viral load parameters (same for all age groups based on image)
vp = 8.1  # Viral load peak
vi = 2    # Viral load incline
vd = -0.17  # Viral load decline
alpha = -7  # Viral shed parameter
beta = 1   # Viral shed parameter
sr = 1.6   # Viral shed factor (Gamma shape parameter for Gamma(1.6, 1/22))
lambda_shed = 22  # Infection rate from viral shed

# Age-specific probability parameters from the table
p_sym_ages = [0.5, 0.55, 0.6, 0.7, 0.83, 0.9]  # Chance to develop symptoms
# Chance to develop severe symptoms
p_sev_ages = [0.02, 0.03, 0.04, 0.07, 0.17, 0.24]
# Chance to develop critical symptoms
p_c_ages = [0.1, 0.11, 0.12, 0.14, 0.33, 0.62]
# Chance to die from critical infection
p_d_ages = [0.12, 0.13, 0.15, 0.26, 0.4, 0.48]


def get_mu_and_sigma(mean, stddev):
    """
    Convert from actual mean and standard deviation to log-normal parameters
    Returns (mu, sigma) for the underlying normal distribution
    """
    mu = np.log(mean * mean / np.sqrt(mean * mean + stddev * stddev))
    sigma = np.sqrt(np.log(1 + stddev * stddev / (mean * mean)))
    return mu, sigma

# Time parameters - actual means and standard deviations from literature
# We need to convert these to log-normal parameters using the conversion function


# Time from Exposed to nonsymptomatic
t_E_mean, t_E_std = 5.1, 1.5  # Example values - replace with actual from literature
t_E_mu, t_E_sigma = get_mu_and_sigma(t_E_mean, t_E_std)
t_E_params = {'mean_log': t_E_mu, 'std_log': t_E_sigma}

# Time to develop symptoms after infection (symptomatic)
t_I_sym_mean, t_I_sym_std = 3.0, 0.9  # Example values
t_I_sym_mu, t_I_sym_sigma = get_mu_and_sigma(t_I_sym_mean, t_I_sym_std)
t_I_sym_params = {'mean_log': t_I_sym_mu, 'std_log': t_I_sym_sigma}

# Time to recover (asymptomatic)
t_I_R_mean, t_I_R_std = 10.0, 2.0  # Example values
t_I_R_mu, t_I_R_sigma = get_mu_and_sigma(t_I_R_mean, t_I_R_std)
t_I_R_params = {'mean_log': t_I_R_mu, 'std_log': t_I_R_sigma}

# Time to develop severe symptoms (severe infection)
t_I_sev_mean, t_I_sev_std = 7.0, 4.9  # Example values
t_I_sev_mu, t_I_sev_sigma = get_mu_and_sigma(t_I_sev_mean, t_I_sev_std)
t_I_sev_params = {'mean_log': t_I_sev_mu, 'std_log': t_I_sev_sigma}

# Time to recover (symptomatic)
t_R_sym_mean, t_R_sym_std = 14.0, 2.0  # Example values
t_R_sym_mu, t_R_sym_sigma = get_mu_and_sigma(t_R_sym_mean, t_R_sym_std)
t_R_sym_params = {'mean_log': t_R_sym_mu, 'std_log': t_R_sym_sigma}

# Time to develop critical symptoms (critical infection)
t_sev_mean, t_sev_std = 4.5, 2.0  # Example values
t_sev_mu, t_sev_sigma = get_mu_and_sigma(t_sev_mean, t_sev_std)
t_sev_params = {'mean_log': t_sev_mu, 'std_log': t_sev_sigma}

# Time to recover (severe infection)
t_R_sev_mean, t_R_sev_std = 21.0, 6.3  # Example values
t_R_sev_mu, t_R_sev_sigma = get_mu_and_sigma(t_R_sev_mean, t_R_sev_std)
t_R_sev_params = {'mean_log': t_R_sev_mu, 'std_log': t_R_sev_sigma}

# Time to die
t_D_mean, t_D_std = 15.0, 4.8  # Example values
t_D_mu, t_D_sigma = get_mu_and_sigma(t_D_mean, t_D_std)
t_D_params = {'mean_log': t_D_mu, 'std_log': t_D_sigma}

# Time to recover (critical infection)
t_R_c_mean, t_R_c_std = 25.0, 6.3  # Example values
t_R_c_mu, t_R_c_sigma = get_mu_and_sigma(t_R_c_mean, t_R_c_std)
t_R_c_params = {'mean_log': t_R_c_mu, 'std_log': t_R_c_sigma}


def viral_load_function(t, vp, vi, vd):
    """
    Simple viral load model: linear rise and decline
    """
    # Peak occurs around day 5, rise with incline vi, decline with rate vd after peak
    t_peak = vp / vi  # Time at which peak occurs
    viral_load = np.where(t <= t_peak,
                          vi * (t),  # Rising phase
                          vp - (vd * (t - t_peak)))  # Declining phase
    return viral_load


def logistic_function(t, s_fp, alpha, beta, vp_t):
    """
    Logistic function: Sp(t) = s_fp / (1 + exp(-(alpha + beta * vp(t))))
    vp_t should be an array of viral load values at time t
    """
    values = 1 / (1 + np.exp(-(alpha + beta * vp_t)))
    return s_fp * values


def gamma_pdf(x, shape, scale):
    """Gamma distribution PDF"""
    return stats.gamma.pdf(x, a=shape, scale=scale)


def lognormal_pdf(x, mean_log, std_log):
    """LogNormal distribution PDF"""
    return stats.lognorm.pdf(x, s=std_log, scale=np.exp(mean_log))


# Create time arrays
t_short = np.linspace(0, 30, 1000)  # For shorter time scales
t_long = np.linspace(0, 70, 1000)   # For longer time scales
t_very_short = np.linspace(0, 10, 1000)  # For very short time scales

# Create a large figure with subplots
fig = plt.figure(figsize=(20, 16))

# 1. Viral shed factor (Gamma distribution)
ax1 = plt.subplot(4, 3, 1)
x_gamma = np.linspace(0, 1, 1000)
gamma_shape = sr  # 1.6
gamma_scale = 1/22  # Based on Gamma(1.6, 1/22)
y_gamma = gamma_pdf(x_gamma, gamma_shape, gamma_scale)
plt.plot(x_gamma, y_gamma, 'b-', linewidth=2,
         label=f'Gamma({gamma_shape}, 1/22)')
plt.title('Viral Shed Factor Distribution', fontweight='bold')
plt.xlabel('x')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 2. Viral load over time and logistic shedding function with multiple shed samples
ax2 = plt.subplot(4, 3, 2)
s_fp = 1.0  # Assuming maximum shedding probability of 1
t_logistic = np.linspace(0, 30, 1000)
vp_t = viral_load_function(t_logistic, vp, vi, abs(vd))

# Generate 100 random viral shed factor samples from Gamma distribution
np.random.seed(42)  # For reproducible results
shed_samples = np.random.gamma(sr, gamma_scale, 100)

# Plot both viral load and multiple shedding probability curves
ax2_twin = ax2.twinx()
ax2.plot(t_logistic, vp_t, 'g-', linewidth=3,
         label='Viral Load V(t)', zorder=10)
max_y_logistic = 0  # To track maximum y value for shedding curves
# Plot 100 shedding curves with different shed factors
for i, shed_factor in enumerate(shed_samples):
    y_logistic_sample = logistic_function(
        t_logistic, shed_factor, alpha, beta, vp_t)
    alpha_val = 1.0 if i < 99 else 1.0  # Make last curve more visible
    linewidth = 1.0 if i < 99 else 1.0
    color = 'lightcoral' if i < 99 else 'lightcoral'
    label = None if i < 99 else 'Shedding Sp(t) samples'
    ax2_twin.plot(t_logistic, y_logistic_sample, color=color,
                  linewidth=linewidth, alpha=alpha_val, label=label)
    max_y_logistic = max(max(y_logistic_sample), max_y_logistic)

# Plot mean shedding curve
mean_shed_factor = np.mean(shed_samples)
y_logistic_mean = logistic_function(
    t_logistic,  mean_shed_factor, alpha, beta, vp_t)
ax2_twin.plot(t_logistic, y_logistic_mean, 'darkred', linewidth=3,
              linestyle='--', label=f'Mean Shedding (sf={mean_shed_factor:.3f})')

ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Viral Load', color='g')
ax2_twin.set_ylabel('Shedding Probability', color='r')
ax2.set_title(
    'Viral Load and Shedding Function\n(100 random shed factor samples)', fontweight='bold')

# Set different y-axis ranges to avoid overlap
ax2.set_ylim(0, max(vp_t) * 1.2)  # Viral load axis
# Shedding probability axis (0-1)
ax2_twin.set_ylim(0,  max_y_logistic * 1.2)

ax2.legend(loc='upper left')
ax2_twin.legend(loc='upper right')
ax2.grid(True, alpha=0.3)

# 3. Time from Exposed to nonsymptomatic
ax3 = plt.subplot(4, 3, 3)
y3 = lognormal_pdf(np.linspace(0, 20, 1000),
                   t_E_params['mean_log'], t_E_params['std_log'])
plt.plot(np.linspace(0, 20, 1000), y3, 'g-', linewidth=2,
         label=f'LogNormal(μ={t_E_params["mean_log"]:.2f}, σ={t_E_params["std_log"]:.2f})')
plt.title(
    f'Time: Exposed → Nonsymptomatic\n(Mean={t_E_mean:.1f}, Std={t_E_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 4. Time to develop symptoms (symptomatic)
ax4 = plt.subplot(4, 3, 4)
y4 = lognormal_pdf(
    t_short, t_I_sym_params['mean_log'], t_I_sym_params['std_log'])
plt.plot(t_short, y4, 'orange', linewidth=2,
         label=f'LogNormal(μ={t_I_sym_params["mean_log"]:.2f}, σ={t_I_sym_params["std_log"]:.2f})')
plt.title(
    f'Time: Nonsymptomatic → Symptomatic\n(Mean={t_I_sym_mean:.1f}, Std={t_I_sym_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 5. Time to recover (asymptomatic)
ax5 = plt.subplot(4, 3, 5)
y5 = lognormal_pdf(t_long, t_I_R_params['mean_log'], t_I_R_params['std_log'])
plt.plot(t_long, y5, 'purple', linewidth=2,
         label=f'LogNormal(μ={t_I_R_params["mean_log"]:.2f}, σ={t_I_R_params["std_log"]:.2f})')
plt.title(
    f'Time: Asymptomatic → Recover\n(Mean={t_I_R_mean:.1f}, Std={t_I_R_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 6. Time to develop severe symptoms
ax6 = plt.subplot(4, 3, 6)
y6 = lognormal_pdf(
    t_long, t_I_sev_params['mean_log'], t_I_sev_params['std_log'])
plt.plot(t_long, y6, 'brown', linewidth=2,
         label=f'LogNormal(μ={t_I_sev_params["mean_log"]:.2f}, σ={t_I_sev_params["std_log"]:.2f})')
plt.title(
    f'Time: Symptomatic → Severe\n(Mean={t_I_sev_mean:.1f}, Std={t_I_sev_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 7. Time to recover (symptomatic)
ax7 = plt.subplot(4, 3, 7)
y7 = lognormal_pdf(
    t_long, t_R_sym_params['mean_log'], t_R_sym_params['std_log'])
plt.plot(t_long, y7, 'pink', linewidth=2,
         label=f'LogNormal(μ={t_R_sym_params["mean_log"]:.2f}, σ={t_R_sym_params["std_log"]:.2f})')
plt.title(
    f'Time: Symptomatic → Recover\n(Mean={t_R_sym_mean:.1f}, Std={t_R_sym_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 8. Time to develop critical symptoms
ax8 = plt.subplot(4, 3, 8)
y8 = lognormal_pdf(t_short, t_sev_params['mean_log'], t_sev_params['std_log'])
plt.plot(t_short, y8, 'cyan', linewidth=2,
         label=f'LogNormal(μ={t_sev_params["mean_log"]:.2f}, σ={t_sev_params["std_log"]:.2f})')
plt.title(
    f'Time: Severe → Critical\n(Mean={t_sev_mean:.1f}, Std={t_sev_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 9. Time to recover (severe infection)
ax9 = plt.subplot(4, 3, 9)
y9 = lognormal_pdf(
    t_long, t_R_sev_params['mean_log'], t_R_sev_params['std_log'])
plt.plot(t_long, y9, 'olive', linewidth=2,
         label=f'LogNormal(μ={t_R_sev_params["mean_log"]:.2f}, σ={t_R_sev_params["std_log"]:.2f})')
plt.title(
    f'Time: Severe → Recover\n(Mean={t_R_sev_mean:.1f}, Std={t_R_sev_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 10. Time to die
ax10 = plt.subplot(4, 3, 10)
y10 = lognormal_pdf(t_long, t_D_params['mean_log'], t_D_params['std_log'])
plt.plot(t_long, y10, 'red', linewidth=2,
         label=f'LogNormal(μ={t_D_params["mean_log"]:.2f}, σ={t_D_params["std_log"]:.2f})')
plt.title(
    f'Time: Critical → Die\n(Mean={t_D_mean:.1f}, Std={t_D_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 11. Time to recover (critical infection)
ax11 = plt.subplot(4, 3, 11)
y11 = lognormal_pdf(t_long, t_R_c_params['mean_log'], t_R_c_params['std_log'])
plt.plot(t_long, y11, 'navy', linewidth=2,
         label=f'LogNormal(μ={t_R_c_params["mean_log"]:.2f}, σ={t_R_c_params["std_log"]:.2f})')
plt.title(
    f'Time: Critical → Recover\n(Mean={t_R_c_mean:.1f}, Std={t_R_c_std:.1f})', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, alpha=0.3)

# 12. Age-specific probability parameters (multiple subplots)
ax12 = plt.subplot(4, 3, 12)

# Create a grouped bar chart for age-specific probabilities
x = np.arange(len(age_groups))
width = 0.2

bars1 = plt.bar(x - 1.5*width, p_sym_ages, width, label='Symptoms', alpha=0.8)
bars2 = plt.bar(x - 0.5*width, p_sev_ages, width, label='Severe', alpha=0.8)
bars3 = plt.bar(x + 0.5*width, p_c_ages, width, label='Critical', alpha=0.8)
bars4 = plt.bar(x + 1.5*width, p_d_ages, width, label='Death', alpha=0.8)

plt.title('Disease Progression Probabilities by Age Group', fontweight='bold')
plt.ylabel('Probability')
plt.xlabel('Age Group')
plt.xticks(x, age_groups, rotation=45)
plt.ylim(0, 1.8)
plt.legend(ncol=2, loc='upper left')
plt.grid(True, alpha=0.3, axis='y')


# Add value labels on bars for better readability


def add_value_labels(bars, values):
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                 f'{value:.2f}', ha='center', va='bottom', fontsize=12, rotation=90)


add_value_labels(bars1, p_sym_ages)
add_value_labels(bars2, p_sev_ages)
add_value_labels(bars3, p_c_ages)
add_value_labels(bars4, p_d_ages)

plt.tight_layout(pad=3.0, h_pad=2.5, w_pad=2.0)
plt.savefig('viral_load_model_parameters.png', dpi=300, bbox_inches='tight')


# Print summary statistics
print("=== MODEL PARAMETERS SUMMARY ===")
print(f"Viral Load Peak (vp): {vp}")
print(f"Viral Load Incline (vi): {vi}")
print(f"Viral Load Decline (vd): {vd}")
print(f"Viral Shed Parameters: α={alpha}, β={beta}")
print(f"Viral Shed Factor: Gamma({sr}, 1/22)")
print(f"Infection Rate from Viral Shed: {lambda_shed}")
print()
print("Age-Specific Disease Progression Probabilities:")
print("Age Group    Symptoms  Severe   Critical  Death")
print("-----------------------------------------------")
for i, age in enumerate(age_groups):
    print(
        f"{age:9s}    {p_sym_ages[i]:.2f}     {p_sev_ages[i]:.2f}     {p_c_ages[i]:.2f}     {p_d_ages[i]:.2f}")
print()
print("Time Distribution Parameters:")
print("Distribution                    Actual Mean  Actual Std   LogNormal μ   LogNormal σ")
print("------------------------------------------------------------------------------")
distributions = [
    ("Exposed → Nonsymptomatic", t_E_mean, t_E_std, t_E_params),
    ("Nonsymptomatic → Symptomatic", t_I_sym_mean, t_I_sym_std, t_I_sym_params),
    ("Asymptomatic → Recover", t_I_R_mean, t_I_R_std, t_I_R_params),
    ("Symptomatic → Severe", t_I_sev_mean, t_I_sev_std, t_I_sev_params),
    ("Symptomatic → Recover", t_R_sym_mean, t_R_sym_std, t_R_sym_params),
    ("Severe → Critical", t_sev_mean, t_sev_std, t_sev_params),
    ("Severe → Recover", t_R_sev_mean, t_R_sev_std, t_R_sev_params),
    ("Critical → Death", t_D_mean, t_D_std, t_D_params),
    ("Critical → Recover", t_R_c_mean, t_R_c_std, t_R_c_params)
]

for name, actual_mean, actual_std, log_params in distributions:
    print(
        f"{name:30s} {actual_mean:8.1f}     {actual_std:8.1f}     {log_params['mean_log']:8.3f}     {log_params['std_log']:8.3f}")
