/*!
 * A template huddle
 */
#include "agent.h"
#include <iostream>
#include <vector>
#include <string>
#include <json/json.h>
#include <json/value.h>
#include <fstream>

using namespace std;
using namespace morph;


vector<int> randPerm(int N){

    vector<int> x(N,0);
    for(int i=0;i<N;i++){
        x[i]=i;
    }
    vector<int> y;
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
    if(distSq < (A->rProx+B->rDist)*(A->rProx+B->rDist)){
        double dABL = getDistSquared(A->brains[B->identity].sensorL.x,A->brains[B->identity].sensorL.y,B->x,B->y);
        if(dABL<(B->rDist*B->rDist)){
            A->brains[B->identity].sensorL.setTmpValue(true);
        }
        double dABR = getDistSquared(A->brains[B->identity].sensorR.x,A->brains[B->identity].sensorR.y,B->x,B->y);
        if(dABR<(B->rDist*B->rDist)){
            A->brains[B->identity].sensorR.setTmpValue(true);
        }
    }
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
        cerr << "\nUsage: ./test simname configfile logdirectory \n\n";
        return -1;
    }

    string paramsfile (argv[2]);
    //ifstream configFile("config.json", std::ifstream::binary);

    string logpath = argv[3];
    morph::Tools::createDir (logpath);


    /*
     * Set up JSON code for reading the parameters
     */

    // Test for existence of the JSON file.
    ifstream jsonfile_test;
    int srtn = system ("pwd");
    if (srtn) {
        cerr << "system call returned " << srtn << endl;
    }

    jsonfile_test.open (paramsfile, ios::in);
    if (jsonfile_test.is_open()) {
        // Good, file exists.
        jsonfile_test.close();
    } else {
        cerr << "json config file " << paramsfile << " not found." << endl;
        return 1;
    }

    // Parse the JSON
    ifstream jsonfile (paramsfile, ifstream::binary);
    Json::Value root;
    string errs;
    Json::CharReaderBuilder rbuilder;
    rbuilder["collectComments"] = false;
    bool parsingSuccessful = Json::parseFromStream (rbuilder, jsonfile, &root, &errs);
    if (!parsingSuccessful) {
        // report to the user the failure and their locations in the document.
        cerr << "Failed to parse JSON: " << errs;
        return 1;
    }

    /*
     * Get simulation-wide parameters from JSON
     */
    srand(root.get ("seed", 0).asUInt());

    const unsigned int storageRate = root.get ("storageRate", 1).asUInt();

    const unsigned int refreshRate = root.get ("refreshRate", 1).asUInt();

    const unsigned int control = root.get ("control", 0).asUInt();

    const unsigned int steps = root.get ("steps", 1000).asUInt();
    //const float speed = root.get ("speed", 0.01).asFloat();
    const float Ta = root.get ("TA", 10.).asFloat();
    const float rA = root.get ("Rarena", 10.).asFloat();

    const Json::Value A = root["agents"];
    unsigned int Ntypes = static_cast<unsigned int>(A.size());

    const float param = root.get ("param", 0.).asFloat();
    //unsigned int units = root.get ("units", 60).asUInt();

    // Instantiate the model object
    vector<agent> agents;
    vector<float> alphas, gammas;
    vector<int> units;
    for(int i=0;i<Ntypes;i++){
        const unsigned int n = A[i].get ("N", 1).asUInt();
        const float xin = A[i].get ("x", -1e9).asFloat();                 // x location
        const float yin = A[i].get ("y", -1e9).asFloat();                 // y location
        const float speed = A[i].get ("speed", 1.).asFloat();           // speed
        const float rProx = A[i].get ("rProx", 1.).asFloat();           // size of agent
        const float rDist = A[i].get ("rDist", 0.).asFloat();           // radius over which it can be sensed
        const unsigned int  ident = A[i].get ("identity", 0).asUInt();  // identity INT
        const float rward = A[i].get ("reward", 0.).asFloat();          // reward value
        const float noise = A[i].get ("noise", 0.).asFloat();           // noise
        const unsigned int  free  = A[i].get ("free", 0.).asUInt();     // free? INT
        const float k1 = A[i].get ("k1", 0.).asFloat();
        const float k2 = A[i].get ("k2", 0.).asFloat();
        const float G = A[i].get ("G", 0.).asFloat();

        for(int j=0;j<n;j++){
            double x=xin;
            double y=yin;
            if(xin==-1e9){
                x = (morph::Tools::randDouble()-0.5)*2.;
            }
            if(yin==-1e9){
                y = (morph::Tools::randDouble()-0.5)*2.;
            }

            agents.push_back(agent((double)x,
                                    (double)y,
                               (double)speed,
                               (double)rProx,
                               (double)rDist,
                               (int)ident,
                               (double)rward,
                               (double)noise,
                               (double)k1,
                               (double)k2,
                               (double)G,
                               (int)free,
                               1000)
                               );


                alphas.push_back(A[i].get ("alpha", 0.0).asFloat());
                gammas.push_back(A[i].get ("gamma", 0.7 ).asFloat());
                units.push_back(A[i].get ("units", 1).asUInt());

            }



    }

    const int N = agents.size();


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





    int stepCount = 0;


    // Create some displays
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);

    double rhoInit = root.get ("camera", 25.).asFloat();;
    string worldName(argv[1]);
    string winTitle = worldName + ": window name";

    if(refreshRate){
        displays.push_back (morph::Gdisplay (600, 600, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
        displays.back().resetDisplay (fix, eye, rot);
        displays.back().redrawDisplay();
    }

    double dt = 0.05;                       // integration time constant

    // Call the init function, which can allocate variables and run
    // any pre-stepping computations.
    try {

    } catch (const exception& e) {
        cerr << "Exception initialising agent object: " << e.what() << endl;
    }


    // Start the loop
    bool finished = false;
    while (!finished) {
        // Step the model
        try {

            // LEARNING
            vector<int> I = randPerm(N);
            for(int i=0;i<N;i++){
                agents[I[i]].resetSensors(Ta);
                for(int j=0;j<N;j++){
                    if(i!=j){
                        communicate(&agents[I[i]],&agents[I[j]]);
                    }
                }
                agents[I[i]].learn(dt);
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
                        agents[i].computeDthetaTemp(200.); // 200.
                        agents[i].computeDthetaNetwork(param);
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

        } catch (const exception& e) {
            cerr << "Caught exception calling agent.step(): " << e.what() << endl;
            finished = true;
        }

        if(refreshRate){
            if (!(stepCount%refreshRate)){
                displays[0].resetDisplay (fix, eye, rot);
                displays[0].drawCylinder(0,0,-0.1,0.,0.,0.,rA,rA,24,vector<double>(3,0.9));
                try {
                    for(int i=0;i<agents.size();i++){
                        agents[i].plot(displays[0]);
                    }
                } catch (const exception& e) {
                    cerr << "Caught exception calling agent.plot(): " << e.what() << endl;
                    finished = true;
                }
                displays[0].redrawDisplay();
            }
        }
        // store the data
        if(storageRate){
            if (!(stepCount%storageRate)) {
                try {
                    for(int i=0;i<agents.size();i++){
                        agents[i].store();
                    }
                } catch (const exception& e) {
                    cerr << "Caught exception calling agent.store(): " << e.what() << endl;
                    finished = true;
                }
            }
        }

        // Halt after how every many iterations suits your model:
        if (stepCount > steps) {
            finished = true;
        }
    }

    // saving here
    stringstream fname;
    fname << logpath << "/agents.h5";
    HdfData data(fname.str());
    for(int i=0;i<agents.size();i++){

        stringstream path;
        path << "/a" << i <<"_";

        stringstream px;
        px<<path.str()<<"x";
        data.add_contained_vals (px.str().c_str(), agents[i].x_store);

        stringstream py;
        py<<path.str()<<"y";
        data.add_contained_vals (py.str().c_str(), agents[i].y_store);

        for(int j=0;j<agents[i].brains.size();j++){
            stringstream pw;
            pw<<path.str()<<"w_"<<j;
            data.add_contained_vals (pw.str().c_str(), agents[i].brains[j].w);
        }
    }

    return 0;
};

