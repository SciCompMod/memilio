#include <gtest/gtest.h>
#include <epidemiology/seir.h>

#include <string>
#include <vector>
#include <fstream>
#include <ios>

using real = double;

std::vector<std::vector<real> > loadCSV(const std::string& filename)
{
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(filename, std::ios::in);


    // Read the Data from the file
    // as String Vector
    std::vector<std::vector<real> > data;
    std::vector<real> row;
    std::string line, word, temp;
    int linecount = 0;
    while (fin >> temp) {

        row.clear();

        // read an entire row and
        // store it in a string variable 'line'
        getline(fin, line);
        linecount++;

        // ignore comments
        if (line[0] == '#') {
            continue;
        }

        // used for breaking words
        std::stringstream s(line);


        // read every column data of a row and
        // store it in a string variable, 'word'
        while (getline(s, word, ' ')) {

            // add all the column data
            // of a row to a vector
            row.push_back(atof(word.c_str()));
        }

        if (row.size() == 5)
            data.push_back(row);

    }

    return data;
}

class TestCompareSeirWithJS : public testing::Test
{
protected:
    void SetUp() override
    {
        refData = loadCSV("data/seir-js-compare.csv");
        t0 = 0.;
        tmax = 50.;
        dt = 1.002003997564315796e-01;

        params.E0 = 10000;
        params.I0 = 1000;
        params.N = 1061000;
        params.R0 = 1000;
        params.a = 1. / 5.2;
        params.b = 2.7;
        params.g = 0.5;

        // add two dampings
        damping<real> d1(0, 1.0);
        damping<real> d2(12, 0.4);
        params.dampings = {d1, d2};
    }

public:


    std::vector<std::vector<real> > refData;
    real t0;
    real tmax;
    real dt;
    seirParam<real> params;
};

TEST_F(TestCompareSeirWithJS, integrate)
{
    EXPECT_EQ(500, refData.size());

    std::vector<std::vector<real> > result(0);

    simulate_seir<real>(t0, tmax, dt, params, result);

    EXPECT_EQ(500, result.size());

    for (size_t irow = 0; irow < 500; ++irow) {
        double t =  refData[irow][0];
        for (size_t icol = 0; icol < 4; ++icol) {
            double ref = refData[irow][icol + 1];
            double actual = result[irow][icol];

            double tol = 1e-6 * ref;
            EXPECT_NEAR(ref, actual, tol);
        }
    }
}
