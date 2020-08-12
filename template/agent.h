/*
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


*/

#include <vector>
#include "morph/HdfData.h"
#include "morph/config.h"
#include "morph/tools.h"

using namespace std;
//using namespace morph;

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

    //void plot(morph::Gdisplay& plt){

    //    plt.drawSphere(x,y,0.,r,col,24);

    //}



    void store(){
        x_store.push_back((float)x);
        y_store.push_back((float)y);
    }

};
