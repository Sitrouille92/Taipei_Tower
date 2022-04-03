#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#define Nmax 500.0
#define Tmax = 10.0
#define m = 1.0
#define k = 20.5
#define c = 1.0
#define Exn1 = 1.0
#define Evn1 = 0.0
#define RK2xn1 = 1.0
#define RK2vn1 = 0.0
#define RK4xn1 = 1.0
#define RK4vn1 = 0.0

void main(){
  float dt = Tmax/Nmax;
  float delta = ((c/m)*(c/m)) - ((4*k)/m);
  float xana = 0.0;
  float t = 0.0;
  float RK2vn2 = 0.0;
  float RK4vn2 = 0.0;
  float k1 = 0, k2 = 0, k3 = 0, k4 = 0;
  FILE *analytique = fopen("analytique.dat","w");
  FILE *euler = fopen("euler.dat","w");

  if (delta >0){
    for (int i = 1; i <= 100; i++){
      xana = exp(-0.5 * t) * (cos(4.5 * t)+ (1/9) * sin(4.5 * t));
      t += dt;
      printf("Analytique : %f\n", xana);
    }
  }else if (delta == 0){
    for (int i = 1; i <= 100; i++){
      xana = exp(-t) * (2 * t + 1);
      t += dt;
      printf("Analytique : %f\n", xana);
    }
  }else{
    for (int i = 1; i <= 100; i++){
      xana = (-1/18) * exp(-t/2) + (19/18) * exp((-19 * t)/2);
      t += dt;
      printf("Analytique : %f\n", xana);
    }
  }
  for (int i = 1; i <= Nmax; i++){
    Evn1 += dt * (-k/m) * Exn1 - (c/m) * Evn1;
    Exn1 += dt * Evn1;
    printf("Euler : U(%f, %f)\n", Exn1, Evn1);
  }
  for (int i = 1; i <= Nmax; i++){
    RK2vn2 = dt * (-k/m) * (RK2xn1 + (dt/2))- (c/m) * (dt/2) * RK2vn1;
    RK2vn1 += RK2vn2;
    RK2xn1 += dt * RK2vn1;
    printf("RK2 Heun : U(%f, %f)\n", RK2xn1, RK2vn1);
  }
  for (int i = 1; i <= Nmax; i++){
    k1 = (-k/m) * RK4xn1 - (c/m) * RK4vn1;
    RK4vn2 = RK4vn1 + 0.5 * dt * k1;
    k2 = (-k/m) * (RK4xn1 + (0.5 * dt)) - (c/m) * RK4vn2;
    RK4vn2 = RK4vn1 + 0.5 * dt * k2;
    k3 = (-k/m) * (RK4xn1 + (0.5 * dt)) - (c/m) * RK4vn2;
    RK4vn2 = RK4vn1 + dt * k3;
    RK4vn1 += (dt/6) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    RK4xn1 += dt * RK4vn1;
    printf("RK4 : U(%f, %f)\n", RK4xn1, RK4vn1);
  }
}