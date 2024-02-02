import numpy as np
import matplotlib.pyplot as plt

# reads a table written by print_table from the given file
def read_from_terminal(filename):
    with open(filename) as file:
        lines = file.readlines()
        skip_to = 0
        labels = ""
        # find first line of the output table
        for i in range(len(lines)):
            if "Time" in lines[i]:
                skip_to = i + 1
                # read labels
                labels = lines[i].split()[1:]
                break
        # throw error if table was not found
        if labels == "":
            raise EOFError("Could not find results table in " + filename)
        result = []
        for i in range(skip_to, len(lines)):
            result+=[[]]
            for txt in lines[i].split():
                try:
                    result[-1] += [float(txt)]
                except ValueError:
                    # remove entries from failed line (should be empty)
                    result = result[:-1]
                    return np.array(result), labels
            if len(result[-1]) != len(labels) + 1:
                # remove entries from failed line
                result = result[:-1]
                return np.array(result), labels
        return np.array(result), labels

def plot_populations(time, metapopulations, labels, name):
    fig_ctr = 1
    lbl_ctr = 0
    for x in metapopulations:
        ys = [np.zeros_like(time)]
        for y in reversed([x[:, i] for i in range(len(x[0]))]):
            ys = [y + ys[0]] + ys

        plt.figure()
        
        for i in range(len(ys) - 1):
            plt.fill_between(time, ys[i], ys[i+1], label=labels[lbl_ctr])
            plt.plot(time, ys[i], c='black')
            lbl_ctr += 1
        plt.plot(time, ys[-1], c='black')
        
        plt.legend()
        plt.savefig(name+str(fig_ctr - 1)+".png")        
        fig_ctr += 1

if __name__ == "__main__":
    table, labels = read_from_terminal("abm_minimal.txt")
    time = table[:,0]
    # subtables = [table[:, 1:4],table[:, 4:7]]
    subtables = [table[:, 1:]]
    # subtables = [table[:, 1 + i * 6: 2+i*6] for i in range(8)]
    plot_populations(time, subtables, labels, "mpm")
