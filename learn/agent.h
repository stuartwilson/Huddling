#include "morph/display.h"
#include "morph/tools.h"
#include "morph/HdfData.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>

using namespace morph;


class sensor{
public:
    double x, y, angle;
    bool onset, offset, val, vPrev, tmpValue;
    sensor(void) {
        x=0.;
        y=0.;
        val=false;
        onset=false;
        offset=false;
        angle=0.;
        tmpValue=false;
        vPrev = false;
    }
    void setAngle(double angle){
        this->angle = angle;
    }
    void setTmpValue(bool v){
        tmpValue = v;
    }
    void setValue(void){
        vPrev = val;
        val = tmpValue;
        onset = (val & !vPrev);
        offset = ((vPrev) && (!val));
    }
};

class brain{
public:
    unsigned int N;
    double alpha, gamma, turn;
    sensor sensorL, sensorR;
    vector<double> w, L, R, L2, R2;
    brain(double angle, unsigned int N,double alpha, double gamma,double weightsigma){
        this->N = N;
        this->alpha = alpha;
        this->gamma = gamma;

        w.resize(N);
        for(int i=0;i<N;i++){
            //w.push_back((morph::Tools::randDouble()-0.5)*weightsigma);
            w[i]= (morph::Tools::randDouble()-0.5)*weightsigma;
        }
        sensorL.setAngle(-angle);
        sensorR.setAngle(+angle);
        turn = 0.;
        L.resize(N,0.);
        R.resize(N,0.);
        L2.resize(N,0.);
        R2.resize(N,0.);
    }
    void update(double reward, double dt){
        L.insert(L.begin(),sensorL.onset*1.);
        L.pop_back();
        R.erase(R.begin());
        R.push_back(sensorR.onset*1.);
        for(int i=0;i<N;i++){
            L2[i] += (L[i]-gamma*L2[i]);// DELTA not *dt
            R2[i] += (R[i]-gamma*R2[i]);// DELTA not *dt
        }
        for (int i=0;i<N;i++){
            w[i] += (alpha*reward*(L2[i]-R2[i]))*dt;
        }
        turn = 0.;
        for (int i=0;i<N;i++){
            turn += w[i]*(L2[i]+R2[i]);
        }
    }
};



class agent
{
public:

    int identity;
    bool free;
    sensor contact;
    double x, y, rProx, rDist, rProxSq, rDistSq, theta, rewardInput, rewardOutput, noise, speed;
    double vx, vy, Tb, Tp, dTheta;
    vector<brain> brains;
    int n;
    vector<double> DK, tau, xk, yk, phik, LR;
    vector<float> x_store, y_store;
    vector<double> col, IDcolor;

    double overdn, piOver2, nnorm, sigma;
    double Tc, A, TL, TR, G, k1, k2;
    agent(double x,double y,double speed, double rProx,double rDist, int identity, double rewardOutput, double noise, double k1, double k2, double G, int free, int n){

        this-> n = n;
        this-> x = x;
        this-> y = y;
        this-> rProx = rProx;
        this-> rDist = rDist;
        this-> identity = identity;
        this-> rewardOutput = rewardOutput;
        this->noise = noise;
        this->speed = speed;
        this->k1 = k1;
        this->k2 = k2;
        this->G = G;
        this->free = (bool)free;
        rProxSq = rProx*rProx;
        rDistSq = rDist*rDist;
        theta = morph::Tools::randDouble()*M_PI*2.;
        dTheta = 0.;

        // thermodynamics-related things
        Tp = 37.0;
        Tb = Tp;
        DK.resize(n);
        tau.resize(n);
        xk.resize(n);
        yk.resize(n);
        phik.resize(n);
        LR.resize(n);

        overdn = 1./(double)n;
        piOver2 = M_PI / 2.;
        nnorm = overdn*0.5;
        sigma = -1./100.;
        for(int k=0;k<n;k++){
            double Phi = (double)k*2.*M_PI/(double)n;//-M_PI;
            phik[k] = Phi;
            xk[k] = rProx*cos(Phi);
            yk[k] = rProx*sin(Phi);
        }


        col.resize(3,0.);
        IDcolor.resize(3,0.9);

        store();
    }

    void add_brain(double angle, int N, double alpha, double gamma, double weightsigma){
        brains.push_back(brain(angle,N,alpha,gamma,weightsigma));
    }

    void move(double dt){

        if(speed>0.){
            double xp = x;
            double yp = y;
            theta += dTheta*dt;
            theta = fmod(theta,2*M_PI);
            x += speed*cos(theta)*dt;
            y += speed*sin(theta)*dt;
        }

        for(int i=0;i<brains.size();i++){
            brains[i].sensorL.x = x+cos(theta+brains[i].sensorL.angle);
            brains[i].sensorL.y = y+sin(theta+brains[i].sensorL.angle);
            brains[i].sensorR.x = x+cos(theta+brains[i].sensorR.angle);
            brains[i].sensorR.y = y+sin(theta+brains[i].sensorR.angle);
        }
    }



    void resetSensors(double Ta){
        rewardInput = 0.;
        contact.setTmpValue(false);
        for(int i=0;i<brains.size();i++){
            brains[i].sensorL.setTmpValue(0);
            brains[i].sensorR.setTmpValue(0);
        }

        for(int k=0;k<n;k++){
            DK[k] = 1e9;
            tau[k] = Ta;
        }
    }

    void learn(double dt){

        contact.setValue();
        for(int i=0;i<brains.size();i++){
            brains[i].sensorL.setValue();
            brains[i].sensorR.setValue();
            brains[i].update(rewardInput*contact.onset, dt);
        }
    }


    void calculateExposure(void){

        Tc=0.;
        int sumContacts = 0;
        for(int k=0;k<n;k++){
            if(DK[k] < 1e9){
                Tc += tau[k];
                sumContacts++;
            }
        }
        if(sumContacts){
            Tc /= (double)sumContacts;
            A = 1.-((double)sumContacts*overdn);
        } else {
            Tc = 0.;
            A = 1.;
        }
    }

/*
    void reorientThermometers(void){

        TL=0.;
        TR=0.;
        for(int k=0;k<n;k++){
            LR[k]=(int)(M_PI-fabs(M_PI-fabs(fmod(theta+piOver2,2.*M_PI)-phik[k]))<piOver2);
            if(LR[k]){
                TL += tau[k];
            } else {
                TR += tau[k];
            }
        }
        TL *= nnorm;
        TR *= nnorm;

    }
    */
void reorientThermometers(void){
         TL=0.;
        TR=0.;
        for(int k=0;k<n;k++){
            LR[k]=round((M_PI-fabs(M_PI-fabs( fmod(theta+piOver2,2.*M_PI)-phik[k]) )<piOver2));

            //LR[k]=round((M_PI-fabs(M_PI-fabs( fmod(theta-phik[k],2.*M_PI))<piOver2));

            //LR[k] = (phik[k]+theta+M_PI*0.5)>(2.*M_PI);

            //LR[k]=round((M_PI-fabs(M_PI-fabs( fmod(theta+piOver2,2.*M_PI-phik[k])) )<piOver2));

            if(LR[k]){
                TL += tau[k];
            } else {
                TR += tau[k];
            }
        }
        TL *= nnorm;
        TR *= nnorm;

    }
//}


    void updateBodyTemp(double Ta, double dt){

        Tb += (k2*(1.-A)*(Tc-Tb)-k1*A*(Tb-Ta)+G)*dt;
    }

    void resetDtheta(void){
        dTheta = 0.;
    }

    void computeDthetaTemp(double Vr){


        double sR = 1./(1.+exp(sigma*(Tp-Tb)*TR));
        double sL = 1./(1.+exp(sigma*(Tp-Tb)*TL));

        //dTheta += atan(Vr*(sL-sR)/(sL+sR));
        dTheta += atan(Vr*(sL-sR)/(sL+sR));

    }

    void computeDthetaNetwork(double Vr){
        double sumTurn = 0.;
        for(int i=0;i<brains.size();i++){
            sumTurn += brains[i].turn;
        }
        dTheta += Vr*(2.0*atan(sumTurn)/M_PI + (morph::Tools::randDouble()-0.5)*noise);
   }

    void resetVelocity(void){

        vx = 0.;
        vy = 0.;
    }

    void applyVelocity(double dt){
        if (free){
            x += vx*dt;
            y += vy*dt;
        }
    }

    void applyBoundary(double ra, double dt){

        double rho = sqrt(x*x+y*y);
        if((rho+rProx)>=ra){
            x += (ra-rho-rProx)*x/rho*dt;
            y += (ra-rho-rProx)*y/rho*dt;
            theta += M_PI;
            theta = fmod(theta,2.*M_PI);
        }
    }

    ~agent (void) {

    }

    void plot(morph::Gdisplay& plt){

       double colscale = fmax(0.,fmin(1.,(Tb-0.)/(50.-0.)));
        col = morph::Tools::getJetColor(colscale);

        plt.drawSphere(x,y,0.,rProx,col,24);
        //plt.drawCylinder(x,y,-0.1,x,y,0.15,rProx,rProx,24,col);
        plt.drawCylinder(x,y,-0.1,x,y,0.1,rDist,rDist,24,IDcolor);

    }

    void store(){
        x_store.push_back((float)x);
        y_store.push_back((float)y);
    }

};
