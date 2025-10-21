
#ifndef __Population__
#define __Population__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"
#include "chromosome.h"

using namespace std;
using namespace arma;

class Population {

private:
//Random rnd; //mi serve in quasi tutte le funzioni
int N=200; //numero di individui nella popolazione (ovvero percorsi possibili)
int Ng=300; //numero di generazioni
double p_cross=0.7; //70% probabilità di fare crossover
vector<Chromosome> pop; 

//arma::field<Chromosome> pop; non so bene come gestirla

public: 
Population (bool crcl, Random& rnd); //se crcl==true crea popolazione sul cerchio, altrimenti nel quadrato
Population (){ }
void swap_chromosomes(int i, int j); //mi serve per poi ordinare i chromosomes
void order(); //ordina i percorsi in ordine crescente
int get_size();
Chromosome& get_element(int i);
void set_element(Chromosome chrom, int i);
void push_back(Chromosome chrom);
bool check(); //invoca il check sui cromosomi e verifica che la prima città di ogni percorso sia la stessa
Chromosome selection(Random& rnd); //seleziona un elemento della popolazione in modo semi random
arma::field<Chromosome> crossover(Chromosome mom, Chromosome dad, Random& rnd);  //passo due percosi che vengono selezionati prima
void average2(int i); // stampa media di L^2 ed errore sulla prima metà di popolazione
Population new_gen(Random& rnd);
// a)funzione che copre tutti gli indici e li accoppia con uno random?
// b)oppure funzione che accoppia random con successivo? 
//per sostituire la generazione conviene mdificare quella esistente o crearla nuova? 
//se a), sono sicura che tutti gli individui vengono sostituiti
//se b), posso ripetere finchè N volte ma non ho la certezza che vengano tutti rimpiazzati
//oppure vedi quaderno

Population evolution(Random&rnd);
};

#endif 

