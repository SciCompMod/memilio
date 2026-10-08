import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("active_pais.txt", index_col=0, parse_dates=True)

plt.figure(figsize=(10, 6))
for column in df.columns:
    plt.plot(df.index, df[column], label=column)
plt.legend(["0-4", "5-14", "15-34", "35-59", "60-79", "80+"], title="Age Groups")
plt.title("Active PAIS Over Time by Age Group")
plt.xlabel("Time in Days since 01.01.1970")
plt.ylabel("Number of Individuals")
plt.show()