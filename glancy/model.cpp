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



void overlap(agent* i, agent* j){
    double r2 = 1.;     // NOTE: assumes radii 1
    double r2x4 = 4.;   // NOTE: assumes radii 1
    double dx, dy, dkx, dky, dk2;
    dx = j->x-i->x;
    dy = j->y-i->y;
    if(dx*dx+dy*dy<=r2x4){
        for(int k=0;k<i->n;k++){
            dkx = j->x-(i->x+i->xk[k]);
            dky = j->y-(i->y+i->yk[k]);
            dk2 = dkx*dkx+dky*dky;
            if(dk2 < r2 && dk2 < i->DK[k]){
                i->DK[k]=dk2;
                i->tau[k]=j->Tb;
            }
        }
    }
}

void springs(agent* i, agent* j){

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


int main (int argc, char **argv)
{
    if (argc < 4) {
        cerr << "\nUsage: ./test simname configfile logdirectory \n\n";
        return -1;
    }

    string paramsfile (argv[2]);
    //ifstream configFile("config.json", std::ifstream::binary);

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

    const unsigned int steps = root.get ("steps", 1000).asUInt();
    const float speed = root.get ("speed", 0.01).asFloat();
    const float Ta = root.get ("TA", 10.).asFloat();
    const float rA = root.get ("Rarena", 10.).asFloat();

    const Json::Value A = root["agents"];
    unsigned int Ntypes = static_cast<unsigned int>(A.size());


    // Instantiate the model object
    vector<agent> agents;
    for(int i=0;i<Ntypes;i++){
        const unsigned int n = A[i].get ("N", 1).asUInt();
        const float xin = A[i].get ("x", -1e9).asFloat();
        const float yin = A[i].get ("y", -1e9).asFloat();
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
            agents.push_back(agent(x,y,(double)k1,(double)k2,(double)G,1000));
        }
    }
    const int N = agents.size();



    string logpath = argv[3];
    morph::Tools::createDir (logpath);


    // Set RNG seed
    int rseed = 1;
    srand(rseed);

    int stepCount = 0;

    // Create some displays
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);

    double rhoInit = 19.;
    string worldName(argv[1]);
    string winTitle = worldName + ": window name";
    displays.push_back (morph::Gdisplay (300, 300, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

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
                for(int i=0;i<N;i++){
                    for(int i=0;i<N;i++){
                    agents[i].resetThermometers(Ta);
                    for(int j=0;j<N;j++){
                        if(i!=j){
                            overlap(&agents[i],&agents[j]);
                        }
                    }
                }

                for(int i=0;i<N;i++){
                    agents[i].contacts();
                }

                for(int i=0;i<N;i++){
                    agents[i].reorient();
                }

                for(int i=0;i<N;i++){
                    agents[i].updateBodyTemp(Ta,dt);
                }

                for(int i=0;i<N;i++){
                    agents[i].move(rA,dt);
                }

                for(int i=0;i<N;i++){
                    agents[i].resetVelocity();
                }

                for(int i=0;i<N;i++){
                    for(int j=0;j<N;j++){
                        if(i!=j){
                            springs(&agents[i],&agents[j]);
                        }
                    }
                }

                for(int i=0;i<N;i++){
                    agents[i].applyVelocity(dt);
                }
            }
            stepCount++;

        } catch (const exception& e) {
            cerr << "Caught exception calling agent.step(): " << e.what() << endl;
            finished = true;
        }

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

        // store the data
        if (stepCount % 10 == 0) {
            try {
                for(int i=0;i<agents.size();i++){
                    agents[i].store();
                }
            } catch (const exception& e) {
                cerr << "Caught exception calling agent.store(): " << e.what() << endl;
                finished = true;
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

        stringstream pb;
        pb<<path.str()<<"b";
        data.add_contained_vals (pb.str().c_str(), agents[i].b_store);

    }

    return 0;
};

