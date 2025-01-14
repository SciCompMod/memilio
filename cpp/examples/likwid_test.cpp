#include <cstddef> 

double work(double* a, size_t n) {
    double s = 0;
    for (size_t j=0; j<n; j++){
        s=s+a[j]*a[j];
    }
    return s;
}

int main()
{
    const size_t n=1e9;
    double s = 0;
    double a[n];
     
    for(size_t i=0; i<1e7; i++){
        s += work(a, n);
    }
}