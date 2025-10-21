
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include "population.h"


using namespace std;
using namespace arma;

Population :: Population(bool crcl, Random& rnd){
  //creo un chromosome 
  Chromosome chrom(crcl, rnd );
  
  //modifico il cromosoma creato per creare N percorsi diversi
  for(int i=0; i<N; i++){
    pop.push_back(chrom.swap(rnd));
  }
}
Chromosome& Population :: get_element(int i){
  return pop.at(i);
}

void Population :: set_element(Chromosome chrom, int i){
  pop.at(i)=chrom;
}

void Population :: push_back(Chromosome chrom){
  pop.push_back(chrom);
}

int Population :: get_size(){
  return N;
}

bool Population::check() {
  
    // Controllo che la prima città sia la stessa in tutti i cromosomi
    City reference_city = pop.at(0).get_city(0);
  
    for (int i = 0; i < N; i++) {
        Chromosome& chrom = pop.at(i);

        // Check che ogni cromosoma sia valido (nessuna città ripetuta)
        if (!chrom.check()) {
            cerr << "Il cromosoma " << i << " contiene città duplicate." << endl;
            return false;
        }

        // Check che la prima città sia uguale alla reference
        if (!(chrom.get_city(0) == reference_city)) {
            cerr << "Il cromosoma " << i << " inizia da una città diversa." << endl;
            return false;
        }
    }

    return true; // Se nessun errore, tutto ok
}


void Population::swap_chromosomes(int i, int j) {
  if (i >= 0 && j >= 0 && i < pop.size() && j < pop.size()) {
      Chromosome c=pop[i];
      pop[i]=pop[j];
      pop[j]=c;
  } else {
      cerr << "Indici non validi nello swap." << endl;
  }
}

void Population::order(){
  for(int i=0; i<N; i++){
    for(int j=i+1; j<N; j++){
      if(pop[i].get_length2()>pop[j].get_length2()) swap_chromosomes(i,j);
    }
  }
  cout<<pop[0].get_length2()<<endl;
}

Chromosome Population :: selection(Random& rnd){
  double r=rnd.Rannyu();
  double p=3.; 
  //cout<<pow(r,p)<<endl;
  int j=(N*pow(r,p));
  //cout<<j<<endl;
  if(j>=N) cout<<"slection out of range, j="<<j<<endl;
  return pop.at(j);
}

arma::field<Chromosome> Population::crossover(Chromosome mom, Chromosome dad, Random& rnd){
  if(rnd.Rannyu()<p_cross){
  arma::field<Chromosome> offspring(2); 
  int dim=mom.get_size();
  cout<<dim<<endl;
  Chromosome child_one(dim), child_two(dim);
  //cout << child_one.get_p_mutation() <<endl;
  City prova=mom.get_city(6);
  int attempts = 0;
  const int max_attempts = 100;
do{
  int cut=int(rnd.Rannyu(1., mom.get_size()-1)); //restituisce l'indice della prima posizione da tagliare
  //cout<<"indice di taglio per cross "<<cut<<endl;
  vector<City> block_mom(dim), block_dad(dim);
  for(int i=0; i<mom.get_size()-cut; i++){
    //cout<<"prendendo  "<<i<<endl;
    block_mom.at(i)=mom.get_city(cut+i);
    block_dad.at(i)=dad.get_city(cut+i);
    //block_mom.push_back(mom.get_city(cut+i)); //vettore con parte finale di mom
    //block_dad.push_back(dad.get_city(cut+i)); //vettore di parte finale di dad
  }
  //cout<<"creati vettori per crossover"<<endl;
  /*for(int j=cut; j<mom.get_size(); j++){
      mom.set_city(j, block_dad.at(j-cut));
      dad.set_city(j, block_mom.at(j-cut));
  }*/
  for(int i=0; i<mom.get_size(); i++){
    if(i<cut){
      child_one.set_city(i, mom.get_city(i));
      child_two.set_city(i, dad.get_city(i));
    }else{
      child_one.set_city(i, block_dad.at(i-cut));
      child_two.set_city(i, block_mom.at(i-cut));
    }
  }

  attempts++;
    if (attempts >= max_attempts) {
        std::cerr << " Raggiunto numero massimo di tentativi nel crossover. Uso i genitori." << std::endl;
        return {mom, dad};
    }
} while (child_one.check()==false || child_two.check()==false);
//cout<<"figli riempiti"<<endl;
  if(child_one.check() && child_two.check()) {
    offspring.at(0)=child_one;
    offspring.at(1)=child_two;
  } else {
    offspring.at(0)=mom;
    offspring.at(1)=dad;
    //cerr << "Uso i genitori" <<endl;
  }
 //cout<<"avvenuto crossover"<<endl;
  return offspring;
  } else return {mom, dad};

}

Population Population :: new_gen(Random& rnd){
  Population new_gen;
  //cout<<"creata popolazione con costruttore di default"<<endl; //funziona
  if(N%2!=0) {
    cout<<"numero dispari di individui!!"<<endl;
  }
  for(int i=0; i<N/2; i++){ //potrei avere problemi se N non è un numero pari
    Chromosome mom=selection(rnd);
    Chromosome dad=selection(rnd);
    //cout<<"genitori selezionati"<<endl;
    arma :: field<Chromosome> offspring=crossover(mom, dad, rnd);
    //cout<<"avvenuto crossover"<<endl;
    
    new_gen.push_back(offspring.at(0));
    new_gen.push_back(offspring.at(1));
    
  }
  return new_gen;
}

void Population :: average2(int i){ //rivedere sia media che errore
  ofstream out("ave_per_gen.dat", ios:: app);
  arma:: field<double> dati(2);
  double ave2=0;
  for(int j=0; j<int(get_size()/2); j++){
    dati.at(0)+=get_element(j).get_length2();
    ave2+=pow(dati.at(0),2);
  }
  dati.at(0)=dati.at(0)/int(get_size()/2);
  ave2=ave2/int(get_size()/2);
  dati.at(1)=sqrt((ave2-pow(dati.at(0),2))/int(get_size()/2)); //palesemente sbagliato è troppo alto

  out<<i+1<<"   "<<dati.at(0)<<"    "<<dati.at(1)<<endl;
  out.close();

  return ;

}

Population Population :: evolution(Random& rnd){
  ofstream out("lunghezza_per_gen.dat");
Population generation=new_gen(rnd); //creo la nuova generazione
  //cout << "Size: " << generation.get_size() << endl;
  for(int i=0; i<Ng; i++){ //numero di generazioni che devo creare
    cout<<"generazione "<<i+1<<endl;
    //muto la generazione
    for(int j=0; j<get_size(); j++){
      //cout << " - Mutating element " << j << endl; 
      
      generation.get_element(j).mutation(rnd); 
      
    }
    //cout<<"mutazioni avvenute per la generazione "<<i+1<<endl;
    //la ordino
    generation.order();
    out<<i+1<<"   " << generation.get_element(0).get_length2()<<endl;
    //cout<<"nuova generazione ordinata"<<endl;
    average2(i);
    for(int k=0; k<get_size();k++){
      pop.at(k)=generation.get_element(k);
    }
    if(i!=Ng-1){
    //ne creo una nuova
    generation=new_gen(rnd);
    //cout<<"nuova generazione creata"<<endl;
    //cout << "Size: " << generation.get_size() << endl;
    }
    
  }
  out.close();
  return generation;
 //CAPIRE STAMPA GRAFICI
}