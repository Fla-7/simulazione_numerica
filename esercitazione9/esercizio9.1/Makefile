CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o city.o chromosome.o population.o
	$(CC) population.o chromosome.o city.o random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
chromosome.o : chromosome.cpp chromosome.h
	$(CC) -c chromosome.cpp -o chromosome.o $(CFLAGS)
population.o : population.cpp population.h 
	$(CC) -c population.cpp -o population.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out best_path.dat ave_per_gen.dat lunghezza2_per_gen.dat
