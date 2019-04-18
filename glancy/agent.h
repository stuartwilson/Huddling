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

class agent
{
public:

    double x, y, r, v, theta, Tb, Tc, A, TL, TR, overdn, nnorm, piOver2;
    double k1, k2, G;
    vector<float> x_store, y_store, b_store;
    vector<double> col;
    int n; // number of thermometers
    vector<double> tau;
    vector<double> DK;
    vector<double> phik, xk, yk;
    vector<int> LR;

    double V, Vr, sigma, Tp, vx, vy;


    agent (double x, double y, double k1, double k2, double G, int n) {

        this->x = x;
        this->y = y;
        this->k1 = k1;
        this->k2 = k2;
        this->G = G;
        this->n = n;
        overdn = 1./(double)n;
        piOver2 = M_PI / 2.;
        nnorm = overdn*0.5;

        theta = morph::Tools::randDouble()*M_PI*2.;
        V = 0.3;
        Vr = 200.;
        sigma = -1./100.;
        Tp = 37.;
        LR.resize(n);
        DK.resize(n);
        tau.resize(n);
        phik.resize(n);
        xk.resize(n);
        yk.resize(n);
        r = 1.;
        for(int k=0;k<n;k++){
            double Phi = (double)k*2.*M_PI/(double)n;//-M_PI; //NOTE THIS IS DIFFERENT FROM SUPPLEMENT
            phik[k] = Phi;
            xk[k] = r*cos(Phi);
            yk[k] = r*sin(Phi);
        }

        Tb = Tp;
        col.resize(3,0.);

        store();
    }


    ~agent (void) {

    }

    /*
    void step(double v){

        this->v = v;
        theta = morph::Tools::randDouble()*2.*M_PI;
        x += v*cos(theta);
        y += v*sin(theta);
        col = morph::Tools::getJetColor(Tb);
    }
    */

    void resetThermometers(double Ta){
        for(int k=0;k<n;k++){
            DK[k] = 1e9;
            tau[k] = Ta;
        }
    }

    void contacts(void){

        Tc=0.;
        int contact = 0;
        for(int k=0;k<n;k++){
            if(DK[k] < 1e9){
                Tc += tau[k];
                contact++;
            }
        }
        if(contact){
            Tc /= (double)contact;
            A = 1.-((double)contact*overdn);
        } else {
            Tc = 0.;
            A = 1.;
        }
    }

    void reorient(){

        TL=0.;
        TR=0.;
        for(int k=0;k<n;k++){
            LR[k]=round((M_PI-fabs(M_PI-fabs(fmod(theta+piOver2,2.*M_PI)-phik[k]))<piOver2));
            if(LR[k]){
                TL += tau[k];
            } else {
                TR += tau[k];
            }
        }
        TL *= nnorm;
        TR *= nnorm;

    }

    void updateBodyTemp(double Ta, double dt){

        Tb += (k2*(1.-A)*(Tc-Tb)-k1*A*(Tb-Ta)+G)*dt;
    }


    void move(double ra, double dt){

        double sR = 1./(1.+exp(sigma*(Tp-Tb)*TR));
        double sL = 1./(1.+exp(sigma*(Tp-Tb)*TL));

        theta += atan(Vr*(sL-sR)/(sL+sR))*dt;

        theta = fmod(theta,2.*M_PI);

        x += cos(theta)*V*dt;
        y += sin(theta)*V*dt;

        double rho = sqrt(x*x+y*y);
        if((rho+r)>=ra){
            x += (ra-rho-r)*x/rho*dt;
            y += (ra-rho-r)*y/rho*dt;
        }
    }

    void resetVelocity(void){

        vx = 0.;
        vy = 0.;
    }

    void applyVelocity(double dt){
        x += vx*dt;
        y += vy*dt;
    }

    void plot(morph::Gdisplay& plt){

        double colscale = fmax(0.,fmin(1.,(Tb-0.)/(50.-0.)));
        col = morph::Tools::getJetColor(colscale);

        plt.drawSphere(x,y,0.,r,col,24);

    }

    void store(void){
        x_store.push_back((float)x);
        y_store.push_back((float)y);
        b_store.push_back((float)Tb);
    }

};
