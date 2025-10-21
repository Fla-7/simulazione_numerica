
#ifndef __City__
#define __City__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"

using namespace std;
using namespace arma;

class City {

private:
//Random rnd;
vector<double> coords;

public: 
City(bool circl, Random& rnd); //costruttore; se crcl==true costruisce una città sul cerchio, altrimenti dentro al quadrato
City() { }
vector <double> get_coords(); //restituisce le coordinate della città 
double get_dist2(City& other);
bool operator==(const City& other) const {
        return coords[0] == other.coords[0] && coords[1] == other.coords[1];
    }

City& operator=(const City& other) {
        if (this != &other) {  // Evita auto-assegnazione
            this->coords = other.coords;
        }
        return *this;
    }


};



#endif 

