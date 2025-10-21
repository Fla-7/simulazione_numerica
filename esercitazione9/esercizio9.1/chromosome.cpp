
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include "chromosome.h"
#include "city.h"

using namespace std;
using namespace arma;

Chromosome :: Chromosome(bool crcl, Random& rnd) {
  for(int i=0; i<n_cities; i++){
    chrom.push_back(City(crcl, rnd));
    //cout<<chrom.at(i).get_coords().at(0)<<endl;;
  }
  return;
}

bool Chromosome :: check(){
  for(int i=0; i<n_cities; i++){
    for (int j = i + 1; j < n_cities; j++){
      if (i!=j && chrom.at(i)==chrom.at(j)) {
        //cerr<<"città ripetuta in posizione "<<i<<" e "<<j<<endl;
        return false;}
    }
  }

  return true;
}

double Chromosome :: get_p_mutation(){
  return p_mutation;
}
int Chromosome :: get_size(){
  return n_cities;
}

City Chromosome :: get_city(int i){
  return chrom.at(i);
}

void Chromosome :: set_city(int i, City city){
  chrom.at(i)=city;
}

double Chromosome :: get_length2(){
  double L=0;
  for(int i=0; i<n_cities-1; i++){
    L+=chrom[i].get_dist2(chrom[i+1]);
  }
  L+=chrom.back().get_dist2(chrom.front());
  return L;
}


Chromosome Chromosome::swap(Random& rnd) {
  Chromosome new_chrom = *this;
    for(int i=0; i<20; i++){
      int k = floor(rnd.Rannyu(1., n_cities - 2));
      City c = new_chrom.chrom.at(k);
      new_chrom.chrom.at(k) = new_chrom.chrom.at(k + 1);
      new_chrom.chrom.at(k + 1) = c;
    }

    return new_chrom;
}

void Chromosome :: pair_permutation(Random& rnd){
  //cout<<"sta avvenendo pair_permutation"<<endl;
  int i = int(rnd.Rannyu(1., n_cities - 3)); //evito che possa scambaire la prima città partendo da uno
  City c=chrom[i];
  chrom[i]=chrom[i+1];
  chrom[i+1]=c;

  //cout<<"avvenuta pair_permutation"<<endl;
}

void Chromosome::shift(Random& rnd) { 
    //cout<<"sta avvenendo shift"<<endl;
    int m = int(rnd.Rannyu(1, n_cities-1)); // numero di città da spostare (ne sposto almeno una ma non tutte perchè la prima rimane fissa)
    int start = int(rnd.Rannyu(1, n_cities - m-1)); //indice punto di partenza del blocco (non include città fissa e lascio il posto per avere le m città da shiftare)
    int shift = int(rnd.Rannyu(1, n_cities - m - start)); // di quante posizioni spostare il blocco
    //cout<<m<<"    "<<start<<"   "<<shift<<endl;
    int dest = start + shift;

    // Copia le città da shiftare
    vector<City> block(chrom.begin() + start, chrom.begin() + start + m);

    // Rimuovi il blocco
    chrom.erase(chrom.begin() + start, chrom.begin() + start + m);

    // Inserisci il blocco alla nuova posizione
    chrom.insert(chrom.begin() + dest, block.begin(), block.end());
  
    if(!check())cout<<"error shift"<<endl;
    //if(check()==true) cout<<"avvenuto shift"<<endl;
    //else cout<<"la funzione shift non funziona"<<endl;
    return;

}



void Chromosome :: permutation(Random& rnd){ 
 //cout<<"sta avvenendo permutazione"<<endl;
  int m = int(rnd.Rannyu(1, n_cities/2)); // numero di città da permutare (la prima rimane fissa e non possono essere più della metà)
  int startA = 1; // punto di partenza del blocco A è fisso a uno 
  int startB= int(rnd.Rannyu(m+1, n_cities-m)); //punto di partenza del blocco B
  cout<<"numero città "<<m<<"   inizio blocco a "<<startA<<"  inizio blocco B"<<startB<<endl;
 
  // Copia le città da permutare
  vector<City> blockA(chrom.begin() + startA, chrom.begin() + startA + m);
  vector<City> blockB(chrom.begin() + startB, chrom.begin() + startB + m);

  // Rimuovi il blocco A e B
  //chrom.erase(chrom.begin() + startA, chrom.begin() + startA + m);
  //chrom.erase(chrom.begin()+startB,chrom.begin() + startB + m );

  //inserisco i blocchi A e B permutandoli
  //chrom.insert(chrom.begin()+startA, blockB.begin(), blockB.end());
  //chrom.insert(chrom.begin()+startB, blockA.begin(), blockA.end());
  for(int i=0; i<m; i++){
    chrom.at(i+startA)=blockB.at(i);
    chrom.at(i+startB)=blockA.at(i);
  }

  if(!check()) cerr<<"m_permutation fallita"<<endl;
    return;
}

void Chromosome :: inversion (Random& rnd){
  //cout<<"sta avvenendo inversione"<<endl;
  int m=int(rnd.Rannyu(1., n_cities)); //numero di città da invertire 
  int start=int(rnd.Rannyu(1., n_cities-m));
  //cout<<m<<"    "<<start<<endl;
  reverse(chrom.begin()+start, chrom.begin()+start+m);
  //cout<<"avvenuta inversione"<<endl;
  if(!check()) cerr<<"inversion fallita"<<endl;

}

void Chromosome :: mutation (Random& rnd){
  double p=rnd.Rannyu();
  //cout <<p<<endl;
  //cout<<p_mutation<<endl;
  //if(p>p_mutation) cout<<"nessuna mutazione"<<endl;
  if(p<p_mutation){
    double t=rnd.Rannyu();
    //cout<<t<<endl;
    if(t<1./4.) {
      pair_permutation(rnd);
      return;
    }else if(t<1./2.){
      shift(rnd);
      return;
    } else if(t<3./4.) {
      permutation(rnd);
      return;
    } else {
      inversion(rnd);
      return;
    }
  }
  return;
}




