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

    const Json::Value A = root["agents"];
    unsigned int N = static_cast<unsigned int>(A.size());


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

    double rhoInit = 15.;
    string worldName(argv[1]);
    string winTitle = worldName + ": window name";
    displays.push_back (morph::Gdisplay (300, 300, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();


    // Instantiate the model object
    vector<agent> agents;
    for(int i=0;i<N;i++){
        const float x = A[i].get ("x", 0.).asFloat();
        const float y = A[i].get ("y", 0.).asFloat();
        agents.push_back(agent((double)x,(double)y));
    }

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
            for(int i=0;i<agents.size();i++){
                agents[i].step(speed);
            }
            stepCount ++;

        } catch (const exception& e) {
            cerr << "Caught exception calling agent.step(): " << e.what() << endl;
            finished = true;
        }

        displays[0].resetDisplay (fix, eye, rot);
        displays[0].drawCylinder(0,0,-0.1,0.,0.,0.,10.,10.,24,vector<double>(3,0.9));
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
    }

    return 0;
};

