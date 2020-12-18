#include <epidemiology_io/mobility_io.h>

// wrapper function to print out matrix entries by gdb's 'print get_element(M,1,1)'
// (GDB doesn't support calling the overloaded operator())
int get_element(Eigen::MatrixXd const& m, int i, int j)
{
    return m(i, j);
}

int main()
{
    // Place text file needs to be in working directory build/examples/ and
    // start from within examples folder in build directory
    Eigen::MatrixXd twitter_migration_2018 = epi::read_mobility_formatted("2018_lk_matrix.txt");
}
