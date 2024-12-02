#include <Wire.h>



// Intervalle de temps entre les mises à jour en secondes
float dt = 0.008; 
const long interval = 8;
unsigned long previousMillis = 0;
unsigned long currentMillis = 0;

float x_hat[5] = {0,0,0,0,0}; // État estimé [Position x,Position y, Vitesse x,Vitesse y,theta]
float P_hat[5][5] = {{1,0,0,0,0}, {0, 1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}}; // Matrice de covariance de l'état

float PN = 0.005; // Bruit du processus(IMU)
float MN = 0.005; // Bruit de mesure(ODO)

//IMU
const int MPU6050_ADDR = 0x68; // MPU-6050 I2C address
//odometry
float xp=0,yp=0,thetap=0,r,DM ;
volatile long right_wheel_pulse_count = 0, left_wheel_pulse_count = 0;
#define R_ENCount 480
#define L_ENCount 480
#define R_CHA 2  //jaune right
#define R_CHB 7  //blanc right + bleu 5v + vert gnd

#define L_CHA 3   // jaune left
#define L_CHB 4   // blanc left + bleu 5v + vert gnd 

#define R_MOTORB 9  //in2 violet
#define R_MOTORA 10  //in1 orange

#define L_MOTORA 6   //in3 gris
#define L_MOTORB 5  //in4 blanc
boolean Direction_right, Direction_left = true;



void setup() {
  Serial.begin(9600);
  //IMU
  Wire.begin();
  MPU6050_init();
  //ODO
  attachInterrupt(digitalPinToInterrupt(R_CHA), right_wheel_pulse, RISING);
  attachInterrupt(digitalPinToInterrupt(L_CHA), left_wheel_pulse, RISING);
  pinMode(R_CHA, INPUT);
  pinMode(R_CHB, INPUT);
  pinMode(L_CHA, INPUT);
  pinMode(L_CHB, INPUT);
  pinMode(R_MOTORA, OUTPUT);
  pinMode(R_MOTORB, OUTPUT);
  pinMode(L_MOTORA, OUTPUT);
  pinMode(L_MOTORB, OUTPUT);
}

void loop() {
  currentMillis = millis();
    if (currentMillis - previousMillis > interval)
  {
    previousMillis = currentMillis;
  // Mesurer l'accélération linéaire et la vitesse angulaire avec l'IMU

  float ax,ay,gyro_z= updateIMU();
  // Mesurer la vitesse linéaire avec l'odométrie (encoder)
  float xp,yp,Thetap = updateODO();

  // Prédiction du modèle dynamique
  predictState(ax, ay,gyro_z);

  // Correction avec les mesures de l'IMU
  correctWithODO(xp,yp,Thetap);

  }

  // Mettre à jour l'affichage ou envoyer les données à d'autres composants
  updateDisplay();

motor_write_r(110);
motor_write_l(110);

}






void predictState(float ax, float ay,float gyro_z)
{
  // Modèle dynamique pour la prédiction
  float A[5][5] = {{1,0, dt,0,0}, {0,1,0,dt,0},{0,0, 1,0,0},{0,0,0,1,0},{0,0,0,0,0}}; // Matrice d'état
  float B[5][3] = {{0.5*dt*dt,0,0}, {0,0.5*dt*dt,0},{dt,0, 0},{0,dt,0},{0,0,1}}; // Vecteur d'entrée
  float u[3]= {ax,ay,gyro_z}; // Commande (vitesse angulaire)

        for (int i = 0; i < 5; ++i) 
        {
            float sum = 0;
            for (int j = 0; j < 5; ++j) 
            {
        sum += A[i][j] * x_hat[j];
            }
            x_hat[i] = sum + B[i][0] * u[0] + B[i][1] * u[1] + B[i][2] * u[2];
        }

// Mise à jour de la matrice de covariance
float AT[5][5];  // Transposée de la matrice A
transpose(A, AT);  // Vous devez implémenter la fonction de transposition

float AP_hat[5][5];  // Produit A * P_hat
matrixProduct(A, P_hat, AP_hat);  // Vous devez implémenter la fonction de produit matriciel

float APT[5][5];  // Produit (A * P_hat) * A^T
matrixProduct(AP_hat, AT, APT);

      for (int i = 0; i < 5; ++i) 
      {
          for (int j = 0; j < 5; ++j) 
          {
          P_hat[i][j] = APT[i][j] + Q[i][j];  // Ajoutez la matrice de bruit du processus
          }
      }
}





void correctWithODO(float xp,float yp,float thetap)
{
  // Mesures de l'IMU

  float z[3] = {xp,yp,thetap}; // Mesures de l'accéléromètre pour l'angle

  // Matrices pour la correction
  float H[3][5] = {{1,0,0,0,0},{0,1,0,0,0},{0,0,0,0,1}};; // Matrice de mesure
  float R[3][3] = {{MN, 0,0}, {0,MN,0},{0, 0,MN}}; // Matrice de covariance de mesure
  float Q[5][5] = {{PN,0,0,0,0}, {0,PN,0,0,0},{0,0,PN,0,0},{0,0,0,PN,0},{0,0,0,0,PN}};
  float K[5][3];
  transpose(H, HT); // Transposée de la matrice H

  // Calcul de la résiduelle et de la matrice de gain de Kalman
  float S[3][3], SI[3][3], Ktemp[5][3], Hx[3], y[3];
  matrixProduct(H, P_hat, Ktemp); // Ktemp = H * P_hat
  matrixProduct(Ktemp, HT, S);    // S = H * P_hat * H^T + R

  inverse(S, SI);    // SI = S^-1
  matrixProduct(P_hat, HT, Ktemp); // Ktemp = P_hat * H^T
  matrixProduct(Ktemp, SI, K);    // K = P_hat * H^T * (H * P_hat * H^T + R)^-1

  // Correction de l'état estimé
  matrixProduct(H, x_hat, Hx); // Hx = H * x_hat
  subtract(z, Hx, y);          // y = z - H * x_hat
  matrixProduct(K, y, Ktemp);  // Ktemp = K * y
  add(x_hat, Ktemp, x_hat);    // x_hat = x_hat + K * (z - H * x_hat)

  // Correction de la matrice de covariance
  matrixProduct(K, H, Ktemp);      // Ktemp = K * H
  subtract(I, Ktemp, Ktemp);       // Ktemp = (I - K * H)
  matrixProduct(Ktemp, P_hat, P_hat); // P_hat = (I - K * H) * P_hat


}






float updateODO()
{   
right_wheel_pulse_count = 0;
left_wheel_pulse_count = 0;
float RPS = (float)(right_wheel_pulse_count / (R_ENCount * 0.08));
float LPS =  (float)(left_wheel_pulse_count  / (L_ENCount * 0.08));
float v_left =LPS*2*3.14*r;
float v_right =RPS*2*3.14*r;
float VL =(v_right+v_left)/2;
float VA =(v_right-v_left)/DM;
xp=xp+VL*dt*cos(thetap);
yp=yp+VL*dt*sin(thetap);
thetap=thetap+VA*dt;
return xp,yp,thetap;
}







void left_wheel_pulse()
{
  int val = digitalRead(L_CHB);

  if (val == LOW)
  {
    Direction_left = false; // Reverse
  }
  else
  {
    Direction_left = true; // Forward
  }

  if (Direction_left)
  {
    left_wheel_pulse_count++;
  }
  else
  {
    left_wheel_pulse_count--;
  }
}






void right_wheel_pulse()
{
  int val = digitalRead(R_CHB);

  if (val == LOW)
  {
    Direction_right = true; // Reverse
  }
  else
  {
    Direction_right = false; // Forward
  }

  if (Direction_right)
  {
    right_wheel_pulse_count++;
  }
  else
  {
    right_wheel_pulse_count--;
  }
}






void motor_write_r(float pwm_mtr)
{
  if (pwm_mtr >= 0)
  {
    analogWrite(R_MOTORB, pwm_mtr);
    analogWrite(R_MOTORA, 0);
  }
  else
  {
    analogWrite(R_MOTORA, abs(pwm_mtr));
    analogWrite(R_MOTORB, 0);
  }
}

void motor_write_l(float pwm_mtr)
{
  if (pwm_mtr >= 0)
  {
    analogWrite(L_MOTORB, pwm_mtr);
    analogWrite(L_MOTORA, 0);
  }
  else
  {
    analogWrite(L_MOTORA, abs(pwm_mtr));
    analogWrite(L_MOTORB, 0);
  }
}






float updateIMU()
{
int16_t ax,ay,az;
  
  // Read accelerometer data
  readAccelerometer(&accelerometer_x, &accelerometer_y, &accelerometer_z);
  
  // Convert raw data to acceleration in m/s^2
  float ax = (float)ax / 16384.0 * 9.81;
  float ay = (float)ay / 16384.0 * 9.81;
  float az = (float)az / 16384.0 * 9.81;
  return ax,ay,gyro_z
}





void MPU6050_init() 
{
  Wire.beginTransmission(MPU6050_ADDR);
  Wire.write(0x6B); // PWR_MGMT_1 register
  Wire.write(0);    // set to zero (wakes up the MPU-6050)
  Wire.endTransmission(true);
}






void readAccelerometer(int16_t* x, int16_t* y, int16_t* z) 
{
  Wire.beginTransmission(MPU6050_ADDR);
  Wire.write(0x3B); // starting with register 0x3B (ACCEL_XOUT_H)
  Wire.endTransmission(false);
  Wire.requestFrom(MPU6050_ADDR, 6, true); // request a total of 6 registers

  // Read accelerometer data
  *x = Wire.read() << 8 | Wire.read();
  *y = Wire.read() << 8 | Wire.read();
 *z = Wire.read() << 8 | Wire.read();
}
