#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <epidemiology/euler.h>

// Test for y'(t) = cos(t)
template <typename T>
void integration_test(size_t n, T tmax, T &err)
{

  std::vector<std::vector<T> > y(n,std::vector<T>(1,0));
  std::vector<std::vector<T> > sol(n,std::vector<T>(1,0));

  std::vector<T> f(1,0);

  T dt =  tmax/n;

  sol[0][0] = std::sin(0);
  sol[n-1][0] = std::sin((n-1)*dt);


  for(size_t i=0;i<n-1;i++)
  {
    sol[i+1][0] = std::sin((i+1)*dt);

    f[0] = std::cos(i*dt);

    explicit_euler(y[i], dt, f, y[i+1]); // 

    // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

    err += std::pow(std::abs(y[i+1][0]-sol[i+1][0]),2.0);

  }

}

int main()
{

    const double pi = std::acos(-1);

    size_t n=100000;
    double tmax = 2*pi;
    double err = 0;
 
    integration_test(n,tmax,err);

    err = std::sqrt(err)/n;

    printf("For n=%d the error is %.4e\n",(int)n, err);


}
