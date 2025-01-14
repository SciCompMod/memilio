#include <vector>

int main()
{
    const size_t n=1e9;
    double s, a[n];
    for(size_t i=0; i<1e5; i++){
        for(size_t j=0; j<n; j++){
            s=s+a[j]*a[j];
        }
    }
}
