#include <epidemiology_io/twitter_migration_io.h>

// wrapper function to print out matrix entries by gdb's 'print get_element(M,1,1)'
// (GDB doesn't support calling the overloaded operator())
int get_element(Eigen::MatrixXi const& m, int i, int j)
{
    return m(i, j);
}

int main()
{
    // Place text file needs to be in working directory build/examples/ and
    // start from within examples folder in build directory
    Eigen::MatrixXi twitter_migration_2018 = epi::read_migration("2018_lk_matrix.txt");
}
