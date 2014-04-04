#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float d2xdt2(int num_masa, float *pos_x, float *pos_y, float *pos_z);
float d2ydt2(int num_masa, float *pos_x, float *pos_y, float *pos_z);
float d2zdt2(int num_masa, float *pos_x, float *pos_y, float *pos_z);
//faltan primeras derivadas, energia y rungekutta

int main(){

  //Datos (en SI):
  double G=6.67*pow(10,-11);
  double R=1.496*pow(10,13);
  double m=1.989*pow(10,30);
  double tiempo=1.577*pow(10,11);
  double dt=10;
  int num_posiciones=int(tiempo/dt);
  int num_masas=3;

  //Posiciones
  float *xs=malloc(sizeof(float)*num_masas);
  float *ys=malloc(sizeof(float)*num_masas);
  float *zs=malloc(sizeof(float)*num_masas);
  
  //Velocidades
  float *mag_vs=malloc(sizeof(float)*num_masas);

  //Masas
  float m_1=m;
  float m_2=2*m;
  float m_3=3*m;

  //Los ubico en el plano x,y,z - Posiciones Iniciales

  ///Masa 1
  xs[0]=0;
  ys[0]=R;
  zs[0]=0;

  ///Masa 2
  xs[1]=-R*pow(3,1/2)/2;
  ys[1]=-R/2;
  zs[1]=0;

  ///Masa 3
  xs[2]=R*pow(3,1/2)/2;
  ys[2]=-R/2;
  zs[2]=0;

  //Velocidades iniciales (magnitudes)

  mag_vs[0]=sqrt(11*G*m/(3*R));
  mag_vs[1]=sqrt(11*G*m/(3*R));
  mag_vs[2]=sqrt(11*G*m/(3*R));
  
  //Arreglos que guardan toda la trayectoria, la energia y el tiempo

  float *pos_x=malloc(sizeof(float)*num_posiciones);
  float *pos_y=malloc(sizeof(float)*num_posiciones);
  float *pos_z=malloc(sizeof(float)*num_posiciones);

  float *energias=malloc(sizeof(float)*num_posiciones);
  float *tiempos=malloc(sizeof(float)*num_posiciones);

  //Ahora a iterar llenando los arreglos

  int i;

  for (i=0;i<tiempo;i+=dt){

    //la idea es iterar en el tiempo, hallar la fuerza en cada masa usando los arreglos xs,ys y zs que guardan la informacion de ellas y la van cambiando
    //a medida que el tiempo avanza. Quiero guardar el avance de las posiciones en otros arreglos pos_x, pos_y, pos_Z

    //evaluar aceleraciones
    //actualizar posiciones
    //actualizar velocidades

  }



  return 0;
}

//FUNCIONES

//Funcion que calcula la segunda derivada de x_n contra el tiempo. Donde n es el num de la masa (1,2,3)

float d2xdt2(int num_masa, float *pos_x, float *pos_y, float *pos_z){
  if(num_masa==1){
    float d_21=G*m_2*(pos_x[1]-pos_x[0])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_31=G*m_3*(pos_x[2]-pos_x[0])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    return d_21 + d_31;
  }
  if(num_masa==2){
    float d_12=G*m_2*(pos_x[0]-pos_x[1])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_32=G*m_3*(pos_x[2]-pos_x[1])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
  if(num_masa==3){
    float d_13=G*m_3*(pos_x[0]-pos_x[2])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    float d_23=G*m_3*(pos_x[1]-pos_x[2])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
}

//Funcion que calcula la segunda derivada de y_n contra el tiempo. Donde n es el num de la masa (1,2,3)

float d2ydt2(int num_masa, float *pos_x, float *pos_y, float *pos_z){
  if(num_masa==1){
    float d_21=G*m_2*(pos_y[1]-pos_y[0])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_31=G*m_3*(pos_y[2]-pos_y[0])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    return d_21 + d_31;
  }
  if(num_masa==2){
    float d_12=G*m_2*(pos_y[0]-pos_y[1])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_32=G*m_3*(pos_y[2]-pos_y[1])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
  if(num_masa==3){
    float d_13=G*m_3*(pos_y[0]-pos_y[2])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    float d_23=G*m_3*(pos_y[1]-pos_y[2])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
}

//Funcion que calcula la segunda derivada de z_n contra el tiempo. Donde n es el num de la masa (1,2,3)

float d2zdt2(int num_masa, float *pos_x, float *pos_y, float *pos_z){
  if(num_masa==1){
    float d_21=G*m_2*(pos_z[1]-pos_z[0])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_31=G*m_3*(pos_z[2]-pos_z[0])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    return d_21 + d_31;
  }
  if(num_masa==2){
    float d_12=G*m_2*(pos_z[0]-pos_z[1])/(pow((pow(pos_x[1]-pos_x[0],2)+pow(pos_y[1]-pos_y[0],2)+pow(pos_z[1]-pos_z[0],2)),3/2));
    float d_32=G*m_3*(pos_z[2]-pos_z[1])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
  if(num_masa==3){
    float d_13=G*m_3*(pos_z[0]-pos_z[2])/(pow((pow(pos_x[2]-pos_x[0],2)+pow(pos_y[2]-pos_y[0],2)+pow(pos_z[2]-pos_z[0],2)),3/2));
    float d_23=G*m_3*(pos_z[1]-pos_z[2])/(pow((pow(pos_x[2]-pos_x[1],2)+pow(pos_y[2]-pos_y[1],2)+pow(pos_z[2]-pos_z[1],2)),3/2));
    return d_12 + d_32;
  }
}

//Funcion que calcula la primera derivada de x_n contra el tiempo. Donde n es el num de la masa (1,2,3)

float dxdt(){

}

//Funcion que calcula la primera derivada de y_n contra el tiempo. Donde n es el num de la masa (1,2,3)

//Funcion que calcula la primera derivada de z_n contra el tiempo. Donde n es el num de la masa (1,2,3)

//Funcion que calcula la energia del sistema a partir de sus posiciones y velocidades

float energia(){

}

//RungeKutta de cuarto orden en x,y,z

float RungeKuttaCuartoOrden_x(int num_masa, double *x_old, double *y_old, double *z_old, float vxs, float vys, float vzs,float dt, float time){
  float *xtemp = x_old;
  float k_1_prime1 =  vxs[num_masa];
  float k_1_prime2 =  d2xdt2(num_masa,x_old,y_old,z_old);
  
  //first step
  time = time + dt/2;
  x_old[num_masa]= x_old[num_masa] + k_1_prime1*dt/2;
  float x2_1 = vxs + k_1_prime2*dt/2;
  float k_2_prime1 = x2_1;
  float k_2_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //second step
  x_old[num_masa] = xtemp[num_masa] + k_2_prime1*dt/2;
  float x2_2 = vxs + k_2_prime2*dt/2;
  float k_3_prime1 = x2_2;
  float k_3_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //third step
  x_old[num_masa] = xtemp[num_masa] + k_3_prime1*dt;
  float x3_2 = vxs + k_3_prime2*dt;
  float k_4_prime1 = x3_2;
  float k_4_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //fourth step
  float k1 = (1/6)*(k_1_prime1+2*k_2_prime1+2*k_3_prime1+k_4_prime1);
  float k2 = (1/6)*(k_1_prime2+2*k_2_prime2+2*k_3_prime2+k_4_prime2);
  
  float x1 = x_old[num_masa] + k1*dt;
  float x2 = vxs + k2*dt;
  
  float r[3] = [time, x1, x2];
  return r
}
float RungeKuttaCuartoOrden_x(int num_masa, double *x_old, double *y_old, double *z_old, float vxs, float vys, float vzs,float dt, float time){
  float *xtemp = x_old;
  float k_1_prime1 =  vxs[num_masa];
  float k_1_prime2 =  d2xdt2(num_masa,x_old,y_old,z_old);
  
  //first step
  time = time + dt/2;
  x_old[num_masa]= x_old[num_masa] + k_1_prime1*dt/2;
  float x2_1 = vxs + k_1_prime2*dt/2;
  float k_2_prime1 = x2_1;
  float k_2_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //second step
  x_old[num_masa] = xtemp[num_masa] + k_2_prime1*dt/2;
  float x2_2 = vxs + k_2_prime2*dt/2;
  float k_3_prime1 = x2_2;
  float k_3_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //third step
  x_old[num_masa] = xtemp[num_masa] + k_3_prime1*dt;
  float x3_2 = vxs + k_3_prime2*dt;
  float k_4_prime1 = x3_2;
  float k_4_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //fourth step
  float k1 = (1/6)*(k_1_prime1+2*k_2_prime1+2*k_3_prime1+k_4_prime1);
  float k2 = (1/6)*(k_1_prime2+2*k_2_prime2+2*k_3_prime2+k_4_prime2);
  
  float x1 = x_old[num_masa] + k1*dt;
  float x2 = vxs + k2*dt;
  
  float r[3] = [time, x1, x2];
  return r
}
float RungeKuttaCuartoOrden_x(int num_masa, double *x_old, double *y_old, double *z_old, float vxs, float vys, float vzs,float dt, float time){
  float *xtemp = x_old;
  float k_1_prime1 =  vxs[num_masa];
  float k_1_prime2 =  d2xdt2(num_masa,x_old,y_old,z_old);
  
  //first step
  time = time + dt/2;
  x_old[num_masa]= x_old[num_masa] + k_1_prime1*dt/2;
  float x2_1 = vxs + k_1_prime2*dt/2;
  float k_2_prime1 = x2_1;
  float k_2_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //second step
  x_old[num_masa] = xtemp[num_masa] + k_2_prime1*dt/2;
  float x2_2 = vxs + k_2_prime2*dt/2;
  float k_3_prime1 = x2_2;
  float k_3_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //third step
  x_old[num_masa] = xtemp[num_masa] + k_3_prime1*dt;
  float x3_2 = vxs + k_3_prime2*dt;
  float k_4_prime1 = x3_2;
  float k_4_prime2 = d2xdt2(num_masa,x_old,y_old,z_old);

  //fourth step
  float k1 = (1/6)*(k_1_prime1+2*k_2_prime1+2*k_3_prime1+k_4_prime1);
  float k2 = (1/6)*(k_1_prime2+2*k_2_prime2+2*k_3_prime2+k_4_prime2);
  
  float x1 = x_old[num_masa] + k1*dt;
  float x2 = vxs + k2*dt;
  
  float r[3] = [time, x1, x2];
  return r
}
float RungeKuttaCuartoOrden_y(int num_masa, double *x_old, double *y_old, double *z_old, float vxs, float vys, float vzs,float dt, float time){
  float *ytemp = y_old;
  float k_1_prime1 =  vys[num_masa];
  float k_1_prime2 =  d2ydt2(num_masa,x_old,y_old,z_old);
  
  //first step
  time = time + dt/2;
  y_old[num_masa]= y_old[num_masa] + k_1_prime1*dt/2;
  float y2_1 = vys + k_1_prime2*dt/2;
  float k_2_prime1 = y2_1;
  float k_2_prime2 = d2ydt2(num_masa,x_old,y_old,z_old);

  //second step
  y_old[num_masa] = ytemp[num_masa] + k_2_prime1*dt/2;
  float y2_2 = vys + k_2_prime2*dt/2;
  float k_3_prime1 = y2_2;
  float k_3_prime2 = d2ydt2(num_masa,x_old,y_old,z_old);

  //third step
  y_old[num_masa] = ytemp[num_masa] + k_3_prime1*dt;
  float y3_2 = vys + k_3_prime2*dt;
  float k_4_prime1 = y3_2;
  float k_4_prime2 = d2ydt2(num_masa,x_old,y_old,z_old);

  //fourth step
  float k1 = (1/6)*(k_1_prime1+2*k_2_prime1+2*k_3_prime1+k_4_prime1);
  float k2 = (1/6)*(k_1_prime2+2*k_2_prime2+2*k_3_prime2+k_4_prime2);
  
  float y1 = y_old[num_masa] + k1*dt;
  float y2 = vys + k2*dt;
  
  float r[3] = [time, x1, x2];
  return r
}
float RungeKuttaCuartoOrden_z(int num_masa, double *x_old, double *y_old, double *z_old, float vxs, float vys, float vzs,float dt, float time){
  float *ztemp = z_old;
  float k_1_prime1 =  vzs[num_masa];
  float k_1_prime2 =  d2zdt2(num_masa,x_old,y_old,z_old);
  
  //first step
  time = time + dt/2;
  z_old[num_masa]= z_old[num_masa] + k_1_prime1*dt/2;
  float z2_1 = vzs + k_1_prime2*dt/2;
  float k_2_prime1 = z2_1;
  float k_2_prime2 = d2zdt2(num_masa,x_old,y_old,z_old);

  //second step
  z_old[num_masa] = ztemp[num_masa] + k_2_prime1*dt/2;
  float z2_2 = vzs + k_2_prime2*dt/2;
  float k_3_prime1 = z2_2;
  float k_3_prime2 = d2zdt2(num_masa,x_old,y_old,z_old);

  //third step
  z_old[num_masa] = ztemp[num_masa] + k_3_prime1*dt;
  float z3_2 = vzs + k_3_prime2*dt;
  float k_4_prime1 = z3_2;
  float k_4_prime2 = d2zdt2(num_masa,x_old,y_old,z_old);

  //fourth step
  float k1 = (1/6)*(k_1_prime1+2*k_2_prime1+2*k_3_prime1+k_4_prime1);
  float k2 = (1/6)*(k_1_prime2+2*k_2_prime2+2*k_3_prime2+k_4_prime2);
  
  float z1 = z_old[num_masa] + k1*dt;
  float z2 = vzs + k2*dt;
  
  float r[3] = [time, x1, x2];
  return r
}
