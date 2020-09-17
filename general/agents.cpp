/*!
 * Glancy et al., (2015) model with communication
 */

#ifdef __OSX__
#include "OpenGL/gl3.h"
#endif

#include "agent.h"
#include "morph/Visual.h"
#include "morph/ColourMap.h"
#include "morph/QuiverVisual.h"
#include "morph/Vector.h"


using morph::ColourMapType;
using morph::Visual;
using morph::QuiverVisual;
using morph::Vector;


typedef morph::VisualDataModel<FLT>* VdmPtr;

std::vector<int> randPerm(int N){

    std::vector<int> x(N,0);
    for(int i=0;i<N;i++){
        x[i]=i;
    }
    std::vector<int> y;
    while(x.size()){
        int i = floor(morph::Tools::randDouble()*x.size());
        y.push_back(x[i]);
        x.erase(x.begin()+i);
    }
    return y;
}

double getDistSquared(double x1, double y1, double x2, double y2){ return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);}


void communicate(agent* A,agent* B){

    double distSq = getDistSquared(A->x,A->y,B->x,B->y);

    if (distSq < (A->rProx+B->rProx)*(A->rProx+B->rProx)){
        A->contact.setTmpValue(true);
        A->rewardInput = B->rewardOutput;

        // TEMPERATURE EXCHANGE ACROSS THERMOMETERS
        double dkx,dky,dk2;
        double r2 = 1.;
        for(int k=0;k<A->n;k++){
            dkx = B->x-(A->x+A->xk[k]);
            dky = B->y-(A->y+A->yk[k]);
            dk2 = dkx*dkx+dky*dky;
            if(dk2 < r2 && dk2 < A->DK[k]){
                A->DK[k]=dk2;
                A->tau[k]=B->Tb;
            }
        }

    }

}

void springs(agent* i, agent* j){

    if (j->free){
        double r = 1.0; // assumes agent radius = 1.
        double r2x4 = 4.;   // NOTE: assumes radii 1
        double dx = j->x-i->x;
        double dy = j->y-i->y;
        double d = dx*dx+dy*dy;
        if((d<=r2x4)){
            double f = fmin(r-sqrt(d)*0.5,r)/sqrt(d);
            j->vx += f*dx;
            j->vy += f*dy;
        }
    }
}

int main (int argc, char **argv)
{
    if (argc < 4) {
        std::cerr << "\nUsage: ./test simname configfile logdirectory \n\n";
        return -1;
    }

     int steps = std::stoi(argv[2]);
    srand(std::stoi(argv[3]));

    std::string logpath = argv[1];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;

    morph::Config conf;
    { std::stringstream ss; ss << logpath <<"/config.json"; conf.init (ss.str()); }
    double dt = conf.getFloat("dt",0.05);



    /*
     * Get simulation-wide parameters from JSON
     */

    const unsigned int storageRate = conf.getUInt ("storageRate", 1);
    const unsigned int refreshRate = conf.getUInt ("refreshRate", 0);
    int control = conf.getUInt ("control", 0);
    float Ta = conf.getFloat ("TA", 10.);
    const float rA = conf.getFloat ("Rarena", 10.);
    const float param = conf.getFloat ("param", 0.);

    const Json::Value A = conf.getArray ("agents");

    // Instantiate the model object
    std::vector<agent> agents;
    std::vector<float> alphas, gammas;
    std::vector<int> units;

    // iterate over agent types
    for(int i=0;i<A.size();i++){
        Json::Value a = A[i];
        int n = a.get ("N", 1).asUInt();
        double xin = a.get ("x", -1e9).asFloat();               // x location
        double yin = a.get ("y", -1e9).asFloat();               // y location
        double speed = a.get ("speed", 1.).asFloat();           // speed
        double rProx = a.get ("rProx", 1.).asFloat();           // size of agent
        double rDist = a.get ("rDist", 0.).asFloat();           // radius over which it can be sensed
        int  ident = a.get ("identity", 0).asUInt();  // identity INT
        double rward = a.get ("reward", 0.).asFloat();          // reward value
        double noise = a.get ("noise", 0.).asFloat();           // noise
        int  free  = a.get ("free", 0.).asUInt();     // free? INT
        double k1 = a.get ("k1", 0.).asFloat();
        double k2 = a.get ("k2", 0.).asFloat();
        double G = a.get ("G", 0.).asFloat();

        // create n agents of type i
        for(int j=0;j<n;j++){
            double x=xin;
            double y=yin;
            if(xin==-1e9){
                x = (morph::Tools::randDouble()-0.5)*2.;
            }
            if(yin==-1e9){
                y = (morph::Tools::randDouble()-0.5)*2.;
            }

            agents.push_back(agent(x,y,speed,rProx,rDist,ident,rward,noise,k1,k2,G,free,1000));
            alphas.push_back(A[i].get ("alpha", 0.0).asFloat());
            gammas.push_back(A[i].get ("gamma", 0.7 ).asFloat());
            units.push_back(A[i].get ("units", 1).asUInt());

        }

    }

    const int N = agents.size();

    for(int i=0;i<agents.size();i++){
        agents[i].setupNetwork(logpath);
    }

    /*
    // add a 'brain' for each unique agent identity (assumes id's are supplied contiguously)
    int agentTypes = -1;
    for(int i=0;i<N;i++){
        if((agents[i].identity)>agentTypes){
            agentTypes = agents[i].identity;
        }
    }
    for(int i=0;i<N;i++){
            for(int j=0;j<=agentTypes;j++){
                agents[i].add_brain(M_PI*0.5,units[i],alphas[i],gammas[i],0.01);
            }
    }
    for(int i=0;i<N;i++){
        agents[i].IDcolor = morph::Tools::getJetColor((double)agents[i].identity/(double)(agentTypes+1));
    }
    */

    int stepCount = 0;

    // SETUP VISUAL
    Visual v(1024, 768, "Visualization");
    unsigned int visId;
    VdmPtr ptr;
    if(refreshRate){
        v.zNear = 0.001;
        v.showCoordArrows = true;
        v.backgroundWhite();
        morph::Vector<FLT, 3> offset = {0.0,0.0,0.0};
        std::vector<morph::Vector<FLT,3>> zerovecs(agents.size());
        visId = v.addVisualModel (new morph::QuiverVisual<FLT> (v.shaderprog, &zerovecs, offset, &zerovecs, ColourMapType::Fixed, 0.1f));
        ptr = (VdmPtr)v.getVisualModel (visId);
    }

    // Call the init function, which can allocate variables and run
    // any pre-stepping computations.
    try {

    } catch (const std::exception& e) {
        std::cerr << "Exception initialising agent object: " << e.what() << std::endl;
    }


    // Start the loop
    bool finished = false;
    while (!finished) {
        // Step the model
        try {


            if(stepCount>steps*0.75){
                control = 2; // SWITCH TO NETWORK IN CHARGE!
                Ta = 10.0;
            }

            // LEARNING

            std::vector<int> I = randPerm(N);
            for(int i=0;i<N;i++){
                agents[I[i]].resetSensors(Ta);
                for(int j=0;j<N;j++){
                    if(i!=j){
                        communicate(&agents[I[i]],&agents[I[j]]);
                    }
                }
                if(!(control==2)){ // DON'T LEARN IF NET IN CHARGE
                    agents[I[i]].learn(dt);
                }
            }


            // THERMODYNAMICS
            for(int i=0;i<N;i++){
                agents[i].calculateExposure();
            }
            for(int i=0;i<N;i++){
                agents[i].reorientThermometers();
            }
            for(int i=0;i<N;i++){
                agents[i].updateBodyTemp(Ta,dt);
            }
            for(int i=0;i<N;i++){
                agents[i].resetDtheta();
            }
            // MOVE
            switch(control){
                case(0):{
                    for(int i=0;i<N;i++){
                        agents[i].computeDthetaTemp(200.); // 200.
                    }
                    break;
                }
                case(1):{
                    for(int i=0;i<N;i++){
                        agents[i].computeDthetaNetwork(param);
                    }
                    break;
                }
                case(2):{
                    for(int i=0;i<N;i++){
                        agents[i].applyNetControl();
                    }
                    break;
                }
            }

            for(int i=0;i<N;i++){
                agents[i].move(dt);
            }

            // COLLISIONS
            for(int i=0;i<N;i++){
                agents[i].resetVelocity();
            }
            //vector<int> I = randPerm(10);
            for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                    if(i!=j){
                        springs(&agents[I[i]],&agents[I[j]]);
                    }
                }
            }
            for(int i=0;i<N;i++){
                agents[i].applyVelocity(dt);
            }

            // BOUNDARY
            for(int i=0;i<N;i++){
                agents[i].applyBoundary(rA,dt);
            }

            stepCount ++;

        } catch (const std::exception& e) {
            std::cerr << "Caught exception calling agent.step(): " << e.what() << std::endl;
            finished = true;
        }

        if(refreshRate){
            if (!(stepCount%refreshRate)){
                std::vector<morph::Vector<FLT, 3>> c;
                std::vector<morph::Vector<FLT, 3>> q;
                for(int i=0;i<agents.size();i++){
                    c.push_back({(FLT)agents[i].x,0.0f,(FLT)agents[i].y});
                    q.push_back({0.0f,0.5f,0.0f});
                }

                ptr->updateData (&c, &q);
                glfwWaitEventsTimeout (0.018);
                v.render();
            }
        }

        // store the data
        if(storageRate){
            if (!(stepCount%storageRate)) {
                try {
                    for(int i=0;i<agents.size();i++){
                        agents[i].store();
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Caught exception calling agent.store(): " << e.what() << std::endl;
                    finished = true;
                }
            }
        }

        // Halt after how every many iterations suits your model:
        if (stepCount > steps) {
            finished = true;
        }
    }


    /*
    // saving here
    std::stringstream fname;
    fname << logpath << "/agents.h5";
    morph::HdfData data(fname.str());
    for(int i=0;i<agents.size();i++){

        std::stringstream path;
        path << "/a" << i <<"_";

        std::stringstream px;
        px<<path.str()<<"x";
        data.add_contained_vals (px.str().c_str(), agents[i].x_store);

        std::stringstream py;
        py<<path.str()<<"y";
        data.add_contained_vals (py.str().c_str(), agents[i].y_store);

        for(int j=0;j<agents[i].brains.size();j++){
            std::stringstream pw;
            pw<<path.str()<<"w_"<<j;
            data.add_contained_vals (pw.str().c_str(), agents[i].brains[j].w);
        }
    }
    */

    // saving here
    std::stringstream fname;
    fname << logpath << "/agents.h5";
    morph::HdfData data(fname.str());

    std::stringstream px;
    px<<"x";

    std::vector<double> X, Y, E, SL, TL;
    for(int i=0;i<agents.size();i++){
        for(int j=0;j<agents[i].x_store.size();j++){
            X.push_back(agents[i].x_store[j]);
        }
        for(int j=0;j<agents[i].y_store.size();j++){
            Y.push_back(agents[i].y_store[j]);
        }

        for(int j=0;j<agents[i].e_store.size();j++){
            E.push_back(agents[i].e_store[j]);
        }

        for(int j=0;j<agents[i].sL_store.size();j++){
            SL.push_back(agents[i].sL_store[j]);
        }

        for(int j=0;j<agents[i].TL_store.size();j++){
            TL.push_back(agents[i].TL_store[j]);
        }

        /*

        //std::stringstream path;
        //path << "/a" << i <<"_";

        //std::stringstream px;
        px<<path.str()<<"x";
        data.add_contained_vals (px.str().c_str(), agents[i].x_store);

        std::stringstream py;
        py<<path.str()<<"y";
        data.add_contained_vals (py.str().c_str(), agents[i].y_store);
        */
        /*
        for(int j=0;j<agents[i].brains.size();j++){
            std::stringstream pw;
            pw<<path.str()<<"w_"<<j;
            data.add_contained_vals (pw.str().c_str(), agents[i].brains[j].w);
        }
        */
    }

    data.add_contained_vals ("x", X);
    data.add_contained_vals ("y", Y);
    data.add_contained_vals ("e", E);
    data.add_contained_vals ("sl", SL);
    data.add_contained_vals ("tl", TL);

    return 0;
};

