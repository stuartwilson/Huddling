#include <vector>
#include "morph/HdfData.h"
#include "morph/config.h"
#include "morph/tools.h"

#include "morph/RecurrentNetwork.h"
#include "morph/RecurrentNetworkTools.h" // only needed for max


class RecurrentNetworkBespoke : public RecurrentNetwork {
public:
    RecurrentNetworkBespoke(void){ }
    void step(std::vector<int> inputID, std::vector<double> inputs, std::vector<int> outputID, std::vector<double> targets, bool learn){
        reset();
        for(int i=0; i<inputID.size(); i++){
            Input[inputID[i]] = inputs[i];
        }
        convergeForward();

        std::fill(J.begin(),J.end(),0.);
        for(int i=0;i<outputID.size();i++){
            J[outputID[i]] = targets[i]-X[outputID[i]];
        }
        if(learn){
            convergeBackward();
            weightUpdate();
        }
    }
};



 // MAYBE NOT NECESSARY NOW DELETED 'BRAIN'
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

class agent
{
public:

    RecurrentNetworkBespoke P;
    std::vector<int> inputID,outputID;
    std::vector<double> inputs, targets;

    int identity;
    bool free;
    sensor contact;
    double x, y, rProx, rDist, rProxSq, rDistSq, theta, rewardInput, rewardOutput, noise, speed, error;
    double vx, vy, Tb, Tp, dTheta;
    int n;
    std::vector<double> DK, tau, xk, yk, phik, LR;
    std::vector<float> x_store, y_store, e_store, sL_store, TL_store;
    std::vector<double> col, IDcolor;

    double overdn, piOver2, nnorm, sigma;
    double Tc, A, TL, TR, G, k1, k2, sL, sR;
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
        nnorm = overdn*2.0; // NOTE THIS WAS *0.5 AND THEN CORRECTED TO *2.0
        sigma = -1./100.;   // NOTE THIS WAS -1/100.
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

    void setupNetwork(std::string logpath){

        morph::Config conf;
        { std::stringstream ss; ss << logpath <<"/network.json"; conf.init (ss.str()); }

        float netdt = conf.getFloat("dt",1.0);
        float tauW = conf.getFloat("tauW",32.0);
        float tauX = conf.getFloat("tauX",1.0);
        float tauY = conf.getFloat("tauY",1.0);
        float divergenceThreshold = conf.getFloat("divergenceThreshold",0.000001);
        int maxConvergenceSteps = conf.getInt("maxConvergenceSteps",400);

        // input/output
        const Json::Value inp = conf.getArray ("inputID");
        for (unsigned int i = 0; i < inp.size(); i++) {
            inputID.push_back(inp[i].asInt());
        }
        const Json::Value out = conf.getArray ("outputID");
        for (unsigned int i = 0; i < out.size(); i++) {
            outputID.push_back(out[i].asInt());
        }
        inputs.resize(inputID.size());
        targets.resize(outputID.size());

        // network connectivity
        std::vector<int> pre, post;
        const Json::Value preIn = conf.getArray ("pre");
        for (unsigned int i = 0; i < preIn.size(); i++) {
            pre.push_back(preIn[i].asInt());
        }
        const Json::Value postIn = conf.getArray ("post");
        for (unsigned int i = 0; i < postIn.size(); i++) {
            post.push_back(postIn[i].asInt());
        }

        int netN = tools::getMax(pre);
        if(tools::getMax(post)>netN){
            netN=tools::getMax(post);
        }
        netN++;

        P.init (netN,netdt,tauW,tauX,tauY,divergenceThreshold,maxConvergenceSteps);
        for(int i=0;i<pre.size();i++){ P.connect(pre[i],post[i]); }
        P.addBias();
        P.setNet();

        P.randomizeWeights(-1.0,+1.0);
    }

    void move(double dt){

        if(speed>0.){
            //double xp = x;
            //double yp = y;
            theta += dTheta*dt;
            theta = fmod(theta,2*M_PI);
            x += speed*cos(theta)*dt;
            y += speed*sin(theta)*dt;
        }

    }



    void resetSensors(double Ta){
        rewardInput = 0.;
        contact.setTmpValue(false);

        for(int k=0;k<n;k++){
            DK[k] = 1e9;
            tau[k] = Ta;
        }
    }

    void learn(double dt){

        contact.setValue();
        std::vector<double> map = getDthetaTemp(200.0);
        ///
        inputs[0] = map[0];//TL;
        inputs[1] = map[1];//TR;
        inputs[2] = map[2];//TB;
        //targets[0] = (map[3]+M_PI)/(2.0*M_PI);
        targets[0] = (map[3] + M_PI/2.0 ) /M_PI;// NEEDS TO BE MAPPED INTO 0 TO 1 +M_PI)/(2.0*M_PI);
        P.step(inputID,inputs,outputID,targets,true);
        error = P.getError();
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
        TL *= nnorm;    // average temperature on left sensors
        TR *= nnorm;    // average temperature on right sensors


        sR = 1./(1.+exp(sigma*(Tp-Tb)*TR));
        sL = 1./(1.+exp(sigma*(Tp-Tb)*TL));

    }
//}


    void updateBodyTemp(double Ta, double dt){
        Tb += (k2*(1.-A)*(Tc-Tb)-k1*A*(Tb-Ta)+G)*dt;
    }

    void resetDtheta(void){
        dTheta = 0.;
    }

    void computeDthetaTemp(double Vr){

        //double sR = 1./(1.+exp(sigma*(Tp-Tb)*TR));
        //double sL = 1./(1.+exp(sigma*(Tp-Tb)*TL));

        //dTheta += atan(Vr*(sL-sR)/(sL+sR));
        dTheta += atan(Vr*(sL-sR)/(sL+sR));

    }

    std::vector<double> getDthetaTemp(double Vr){ // HACK: NEEDS TO BE SAME AS COMPUTE TO PROVIDE NET WITH A COPY OF THE COMMANDS TO LEARN

        //double sR = 1./(1.+exp(sigma*(Tp-Tb)*TR));
        //double sL = 1./(1.+exp(sigma*(Tp-Tb)*TL));

        //dTheta += atan(Vr*(sL-sR)/(sL+sR));

        std::vector<double> map (3);
        map[0] = TR;
        map[1] = TL;
        map[2] = Tb;
        map[3] = atan(Vr*(sL-sR)/(sL+sR));

        return map;

    }

    void applyNetControl(void){

        std::vector<double> map = getDthetaTemp(200.0);

        ///
        inputs[0] = map[0];//TL;
        inputs[1] = map[1];//TR;
        inputs[2] = map[2];//Tb;
        targets[0] = map[3];
        P.step(inputID,inputs,outputID,targets,false); // do forward dynamics to get output

        //dTheta += P.X[outputID[0]]*2.0*M_PI-M_PI;
        dTheta += P.X[outputID[0]]*M_PI-(M_PI/2.0);

    }


    void computeDthetaNetwork(double Vr){

        /*
        double sumTurn = 0.;
        for(int i=0;i<brains.size();i++){
            sumTurn += brains[i].turn;
        }
        dTheta += Vr*(2.0*atan(sumTurn)/M_PI + (morph::Tools::randDouble()-0.5)*noise);
        */
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

    /*
    void plot(morph::Gdisplay& plt, double zOffset){

       double colscale = fmax(0.,fmin(1.,(Tb-0.)/(50.-0.)));
        col = morph::Tools::getJetColor(colscale);

        plt.drawSphere(x,y,zOffset,rProx,col,24);
        //plt.drawCylinder(x,y,-0.1,x,y,0.15,rProx,rProx,24,col);
        plt.drawCylinder(x,y,zOffset-0.1,x,y,zOffset+0.1,rDist,rDist,24,IDcolor);

    }
    */

    void store(){
        x_store.push_back((float)x);
        y_store.push_back((float)y);
        e_store.push_back((float)error);
        sL_store.push_back((float)sL);
        TL_store.push_back((float)TL);
    }

};
