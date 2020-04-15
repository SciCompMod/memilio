#ifndef EULER_H
#define EULER_H

#include <vector>

/**
 * Simple explicit euler integration y(t+1) = y(t) + h*f(t,y)
 * for ODE y'(t) = f(t,y)
 * @param[in] yt value of y at t, y(t)
 * @param[in] dt current time step h=dt
 * @param[in] f right hand side of ODE f(t,y)
 * @param[out] ytp1 approximated value y(t+1)
 */
template <typename T>
void explicit_euler(std::vector<T> const &yt, const T dt, std::vector<T> const &f, std::vector<T> &ytp1) {

  for (size_t i=0;i<yt.size();i++) 
  {
    ytp1[i] = yt[i] + dt * f[i];
  }

}

#endif // EULER_H
