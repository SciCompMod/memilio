#include <string>
#include <vector>
#include <fstream>
#include <sstream>

/**
 * @brief read table of space separated numerical values from a file
 */
template <class Real>
std::vector<std::vector<Real>> load_test_data_csv(const std::string& filename)
{
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(filename, std::ios::in);

    // Read the Data from the file
    std::vector<std::vector<Real>> data;
    std::vector<Real> row;
    std::string line, word;
    while (getline(fin, line)) {
        row.clear();

        // ignore comments
        if (line[0] == '#') {
            continue;
        }

        //read columns in this row
        std::stringstream s(line);
        double v;
        while (s >> v) {
            row.push_back(v);
        }

        //append non empty rows
        if (row.size() > 0) {
            data.push_back(row);
        }
    }

    return data;
}