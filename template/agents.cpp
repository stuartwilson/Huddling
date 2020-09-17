/*!
 * A template huddle
 */

#ifdef __OSX__
#include "OpenGL/gl3.h"
#endif

#include "agent.h"
#include "morph/Visual.h"
#include "morph/ColourMap.h"
#include "morph/QuiverVisual.h"
#include "Vector.h"

using morph::ColourMapType;
using morph::Visual;
using morph::QuiverVisual;
using morph::Vector;


typedef morph::VisualDataModel<FLT>* VdmPtr;


int main (int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "\nUsage: ./build/agents logdirectory seed \n\n";
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
    double dt = conf.getFloat("dt",1.0);


    double speed = conf.getFloat ("speed", 0.01);

    std::vector<agent> agents;

    const Json::Value A = conf.getArray ("agents");
    //int N = static_cast<unsigned int>(A.size());
    for (unsigned int i = 0; i < A.size(); i++) {
        Json::Value a = A[i];
        double x = a["x"].asFloat();
        double y = a["y"].asFloat();
        agents.push_back(agent(x,y));
    }

    // SETUP VISUAL
    Visual v(1024, 768, "Visualization");
    v.zNear = 0.001;
    v.showCoordArrows = true;
    v.backgroundWhite();
    morph::Vector<FLT, 3> offset = {0.0,0.0,0.0};
    std::vector<morph::Vector<FLT,3>> zerovecs(agents.size());

    unsigned int visId = v.addVisualModel (new morph::QuiverVisual<FLT> (v.shaderprog, &zerovecs, offset, &zerovecs, ColourMapType::Fixed, 0.1f));

    VdmPtr ptr = (VdmPtr)v.getVisualModel (visId);

    // Start the loop
    unsigned int stepCount = 0;
    bool finished = false;
    while (!finished) {
        try {
            for(int i=0;i<agents.size();i++){
                agents[i].step(speed);
            }

            std::vector<morph::Vector<FLT, 3>> c;
            std::vector<morph::Vector<FLT, 3>> q;
            for(int i=0;i<agents.size();i++){
                c.push_back({(FLT)agents[i].x,0.0f,(FLT)agents[i].y});
                q.push_back({0.0f,0.5f,0.0f});
            }

            ptr->updateData (&c, &q);
            glfwWaitEventsTimeout (0.018);
            v.render();

            stepCount ++;

        } catch (const std::exception& e) {
            std::cerr << "Caught exception calling agent.step(): " << e.what() << std::endl;
            finished = true;
        }

        // store the data
        if (stepCount % 10 == 0) {
            try {
                for(int i=0;i<agents.size();i++){
                    agents[i].store();
                }
            } catch (const std::exception& e) {
                std::cerr << "Caught exception calling agent.store(): " << e.what() << std::endl;
                finished = true;
            }
        }

        // Halt after how every many iterations suits your model:
        if (stepCount > steps) {
            finished = true;
        }
    }

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
    }

    return 0;

};
