
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
//RIVEDERE QUESTO CONTROLLO 
/*bool Population :: check(){
  bool ctrl_chrom=false;
  bool ctrl_1=false; 
  for(int i=0; i<N; i++){
   //devo controllare che la prima città sia uguale per tutti, forse meglio fare una funzione nei cromosomi 
   Chromosome c1=pop.at(i);
   Chromosome c2=pop.at(i+1); 
   if(c1.get_city(0)==c2.get_city(0)){
    ctrl_chrom=pop.at(i).check(); //check sui cromosomi (verifica che non ci sono città ripetute)
   }else {
    cerr<<"i percorsi "<<i<<" e "<<i+1<<"iniziano da città diverse"<<endl;
    ctrl_1=false;
   }
  }

  if(ctrl_chrom && ctrl_1) return true;
  else return false;
  
}*/

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
  //cout<<pop[0].get_length2()<<endl;
}

Chromosome Population :: selection(Random& rnd){
  double r=rnd.Rannyu();
  double p=3.; //capire come scegliere P
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
  //cout<<dim<<endl;
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
    //cout<<"MIAO"<<endl;
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
        //cerr << "⚠️ Raggiunto numero massimo di tentativi nel crossover. Uso i genitori." << std::endl;
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
    cerr << "Uso i genitori" <<endl;
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

void Population :: migration(int rank, int size){
  int n_city=get_element(0).get_size(); 
  int dest=(rank+1)%size; //destinatario
  int source=(rank-1+size)%size; //mittente
  //cout<<n_city<<endl;
  

  //serializzo i dati da inviare 
  int size_tot=N_migr*2*n_city; //numero di elememti da migrare * numero di città nel cromosoma * 2 
  vector<double> send_buffer;
  
  for(int i=0; i<N_migr; i++){
    for(int j=0; j<get_element(i).get_size(); j++){
auto coords = get_element(i).get_city(j).get_coords();
  // Stampa di debug aggressiva (verrà stampata prima del crash)
  /*if(rank==3){
        cerr << "Rank " << rank << " | Ind " << i << ", Citta " << j << " | Size: " << coords.size() << endl;
        
        if (coords.size() < 2) {
            // ESEGUI L'ABORT IMMEDIATAMENTE QUI, SENZA ASPETTARE .at(1)
            cerr << "!!! Rank " << rank << " DETECTED CORRUPT DATA!" << endl;
            // opzionale: throw std::runtime_error("Corrupt data found");
        }
      }*/
      //send_buffer.insert(send_buffer.end(), get_element(i).get_city(j).get_coords().begin(), get_element(i).get_city(j).get_coords().end());
      send_buffer.push_back(get_element(i).get_city(j).get_coords().at(0));
      send_buffer.push_back(get_element(i).get_city(j).get_coords().at(1));
      
    }
  }
  //if (rank==0)cerr<<rank<<" "<<"send_ buffer riempito"<<endl;
 /*for(int k=0; k<size_tot; k+=2){
  for(int i=0; i<N_migr; i++){
    for(int j=0; j<n_city; j++){
      //send_buffer.insert(send_buffer.end(), get_element(i).get_city(j).get_coords().begin(), get_element(i).get_city(j).get_coords().end());
      send_buffer.at(k)=get_element(i).get_city(j).get_coords().at(0);
      send_buffer.at(k+1)=get_element(i).get_city(j).get_coords().at(1);
    }
  }
 }*/
/*int b_index=0;
 for(int i = 0; i < N_migr; i++){
    
    // 3. CICLO SULLE CITTÀ (j)
    for(int j = 0; j < n_city; j++){
      
      // SCRIVI LA COORDINATA X E AVANZA L'INDICE
      send_buffer.at(b_index) = get_element(i).get_city(j).get_coords().at(0);
      b_index++; 
      
      // SCRIVI LA COORDINATA Y E AVANZA L'INDICE
      send_buffer.at(b_index) = get_element(i).get_city(j).get_coords().at(1);
      b_index++;
    }
}
cout<<"indice "<<b_index<<endl;*/

  //cout<<send_buffer.size()<<endl;
  //cout<<size_tot<<endl;

vector<double> recv_buffer(size_tot);

MPI_Request send_request, recv_request;
MPI_Status status;

// Inizia la ricezione: ricevi dal 'source' nel 'recv_buffer'
MPI_Irecv(recv_buffer.data(), size_tot, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &recv_request);

// Inizia l'invio: invia il 'send_buffer' al 'dest'
MPI_Isend(send_buffer.data(), size_tot, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &send_request);

MPI_Wait(&send_request, &status);
MPI_Wait(&recv_request, &status);

//if(rank==0) cerr<<rank<<" avvenuti passaggi MPI"<<endl;

//deserializzo gli individui
pop.erase(pop.begin() + (pop.size() - N_migr), pop.end());//libero gli ultimi 50 posti
//cout<<rank<<" liberati ultimi posti"<<endl;
int buffer_index=0;
for(int k=0; k<N_migr; k++){
  Chromosome new_chrom(n_city); // ricreo i cromosomi dai double ricevuti nel buffer
  for(int c=0; c<n_city; c++){
    double x=recv_buffer.at(buffer_index);
    buffer_index++;
    //cout<<buffer_index<<endl;
    double y=recv_buffer.at(buffer_index);
    buffer_index++; //per passare alla città successiva
    //cerr<<rank<<" prese coordinate per la città "<<c<<"dell'individuo "<<k<<endl;
    City new_city(x,y);
    //cerr<<rank<<"città creata"<<endl;
    new_chrom.set_city(c, new_city);
    //cerr<<rank<<"città settata"<<endl;
    
  }
  //cerr<<"ERRORE"<<endl;
  pop.push_back(new_chrom);
}

//cerr<<rank<<"avvenuta de serializzazione"<<endl;

//Attesa dell'invio (garantisce che il buffer di invio sia riutilizzabile)
    MPI_Wait(&send_request, &status);
    
    // Stampa di conferma (opzionale)
    //cerr << "Rank " << rank << ": Migrazione completata. Inviato a " << dest 
              //<< ", Ricevuto da " << source << "." << endl;
    //fflush(stderr);

  //mando i primi 50 individui negli ultimi del rank successivo 
  //dal rank presente riempio un vector in cui copio i primi 50 cromosomi
  //lo spedisco al rank successivo e lo copio negli ultimi 50 posti
}

Population Population :: evolution(Random& rnd, int rank, int size){
  
  ofstream out("lunghezza_per_gen_"+to_string(rank)+".dat");
Population generation=new_gen(rnd); //creo la nuova generazione
  //cout << "Size: " << generation.get_size() << endl;
  for(int i=0; i<Ng; i++){ //numero di generazioni che devo creare
    if(rank==0) cout<<"generazione "<<i+1<<endl;
    //muto la generazione
    for(int j=0; j<get_size(); j++){
      //cout << " - Mutating element " << j << endl; 
      
      generation.get_element(j).mutation(rnd,rank); 
      
    }
    //cout<<"mutazioni avvenute per la generazione "<<i+1<<endl;
    //la ordino
    generation.order();
    
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
  if((i + 1) % Ng_migr == 0){
  //ho creato la nuova generazione e l'ho ordinata, ora avviene la migrazione 
      generation.migration(rank, size); //capire se devo scrivere generation.migration
      //cout<<"avvenuta migrazione"<<endl;
    //riordinare poi 
    generation.order();
    
  }
  out<<i+1<<"   " << generation.get_element(0).get_length2()<<endl; //per ogni rank stampo l'andamendo delle lunghezze
    
  }
  
  out.close();
  return generation;

}