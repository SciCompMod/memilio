#include <epidemiology/secir.h>

int main()
{
  double t0 = 0;
  double tmax = 1;
  double dt = 0.1;

  int model = 1;

  if(model == 0)
  {

    seirParam<double> params = seirParam<double>();

    printSeirParams(params);

    std::vector<std::vector<double> > seir(0);

    simulate(t0, tmax, dt, params, seir);
  }else{

    double tinc = 5.2, // R_2^(-1)+R_3^(-1)
           tinfmild = 6, // 4-14  (=R4^(-1))
           tserint = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
           thosp2home = 12, // 7-16 (=R5^(-1))
           thome2hosp = 5, // 2.5-7 (=R6^(-1))
           thosp2icu = 2, // 1-3.5 (=R7^(-1))
           ticu2home = 8, // 5-16 (=R8^(-1))
           tinfasy = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
           ticu2death = 5; // 3.5-7 (=R5^(-1))

    double tinfasy2 = 1.0/(0.5/(tinfmild-tserint)+0.5/tinfmild);
    if(fabs(tinfasy2-tinfasy)>0)
    {
        printf("\n ----> TODO / To consider: In the HZI paper, tinfasy (the asymptomatic infectious time) or R9^(-1)=R_3^(-1)+0.5*R_4^(-1) is directly given by R_3 and R_4 and maybe should not be an 'additional parameter'\n"); 
    }

    double cont_freq = 0.5, // 0.2-0.75
            alpha = 0.09, // 0.01-0.16
            beta = 0.25, // 0.05-0.5
            delta = 0.3, // 0.15-0.77
            rho = 0.2, // 0.1-0.35
            theta = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000,
            nb_exp_t0 = 100,
            nb_inf_t0 = 50,
            nb_car_t0 = 50,
            nb_hosp_t0 = 20,
            nb_icu_t0 = 10,
            nb_rec_t0 = 10,
            nb_dead_t0 = 0;

    seirParam<double> params = seirParam<double>(tinc, tinfmild, tserint, thosp2home, thome2hosp, thosp2icu, ticu2home, tinfasy, ticu2death, cont_freq, alpha, beta, delta, rho, theta,
                  nb_total_t0, nb_exp_t0, nb_car_t0, nb_inf_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0) ;

    printSeirParams(params);

    std::vector<std::vector<double> > seir(0);

    simulate(t0, tmax, dt, params, seir);

  }

  



}
