import pandas as pd
import matplotlib.pyplot as plt

# --- Load data ---
df = pd.read_csv("control_parameters.csv")  # replace with your filename if needed
df = pd.read_csv("control_parameters.csv", skipinitialspace=True)


# --- Create plot ---
plt.figure(figsize=(10, 6))

# Columns to plot, HomeOffice first
cols = [
    "School closure",
    "Face masks & social distancing School",
    "Face masks & social distancing Work",
    "Face masks & social distancing Other",
    "Remote work"
]

lw = 3

plt.step(df["Time"], df["Remote work"], where="post", label="HomeOffice", color="tab:red", linewidth=lw)
plt.step(df["Time"], df["School closure"], where="post", label="SchoolClosure", color="tab:brown", linewidth=lw)
plt.step(df["Time"], df["Face masks & social distancing School"], where="post", label="PhysicalDistancingSchool", color="tab:blue", linewidth=lw)
plt.step(df["Time"], df["Face masks & social distancing Work"], where="post", label="PhysicalDistancingWork", color="tab:orange", linewidth=lw, linestyle='--',
    dashes=(2, 2))
plt.step(df["Time"], df["Face masks & social distancing Other"], where="post", label="PhysicalDistancingOther", color="tab:green", linewidth=lw, linestyle='--',
    dashes=(2, 2))

# --- Styling ---
plt.xlabel("Time (days)")
plt.ylabel("Control intensity")
plt.title("Weekly Control Measures")
# plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(loc="best")
plt.tight_layout()

# --- Show plot ---
plt.show()
plt.savefig("control_measures.png", dpi=300)
