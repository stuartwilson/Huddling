#include <vector>
#include "morph/HdfData.h"
#include "morph/config.h"
#include "morph/tools.h"

class agent
{
public:

    double x, y, r, v, theta, TB;
    std::vector<float> x_store, y_store;
    std::vector<double> col;

    agent (double x, double y) {

        this->x = x;
        this->y = y;
        r = 1.;
        TB = 0.5;

        store();
    }


    ~agent (void) {

    }

    void step(double v){

        this->v = v;
        theta = morph::Tools::randDouble()*2.*M_PI;
        x += v*cos(theta);
        y += v*sin(theta);
        col = morph::Tools::getJetColor(TB);
    }

    void store(){
        x_store.push_back((float)x);
        y_store.push_back((float)y);
    }

};
