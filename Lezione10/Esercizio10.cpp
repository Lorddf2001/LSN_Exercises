#include "mpi.h"
#include "map.h"
#include "route.h"
#include "population.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Random shared by ranks
    Random rnd_shared;
    int seed[4];
    ifstream Primes, Seed;
    int p1, p2, ptemp;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd_shared.SetRandom(seed,p1,p2);
    Seed.close();

    //Random for each rank
    Random rnd;
    Primes.open("Primes");
    for(int i = 0; i < rank; i++){      //change primes for each rank
        Primes >> ptemp >> ptemp;       
    }
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    seed[3] += rank;                    //change seed for each rank
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    //Input 
    ifstream Input;
    Input.open("input.in");
    int migrate, n_population, nstep, n_migration;
    Input >> migrate;
    Input >> n_population;
    Input >> nstep;
    Input >> n_migration;
    Input.close();
    
    //Simulation
    Population pop(n_population, rnd);             //population
    Route sending_order(rnd_shared, size);         //order list for sending best individuals
    Route best;                                    //Best route
    vector<int> best_route;                        //Decompose best route
    int best_size;
    int current_rank, send_to, received_from;

    if (migrate != 0 )  //without migration
    {     

        ofstream Length;
        Length.open("output_lengthmean_wom_" + to_string(rank) + "_.dat");  

        for(int i  = 0; i < nstep; i++){
        pop.nextgen();
        pop.mutation();
        pop.sorting();
        Length << pop.mean_length() << endl;
        }
        Length.close();

    }

    else      //with migration

    {
        
        ofstream Length;
        Length.open("output_lengthmean_wm_" + to_string(rank) + "_.dat");

        for(int i  = 0; i < nstep; i++){

            pop.nextgen();
            pop.mutation();
            pop.sorting();
            Length << pop.mean_length() << endl;

            //Exchange best individuals
            if(i % n_migration == 0 && i != 0){

                //First send-receive: The last one sends to the first one and receves from the last but one.
                current_rank = sending_order.getroute()[size - 1];             
                send_to = sending_order.getroute()[0];                 // Rank to which the current rank sends the route
                received_from = sending_order.getroute()[size - 2];    // Rank from which current rank receives
                if(rank == current_rank){
                    best = pop.getbest();
                    best_route = best.getroute();
                    best_size = best.getncity();
                    // Receive the best route data from "received from"
                    MPI_Recv(best_route.data(), best_size, MPI_INT, received_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // Send the best route data to "send_to"
                    MPI_Send(best_route.data(), best_size, MPI_INT, send_to, 0, MPI_COMM_WORLD);

                    // Update the best individual with the received route
                    best.updateroute(best_route);
                    pop.exchangebest(best);
                }

                //Second send-receive: The first one sends to the second one and receves from the last one.
                current_rank = sending_order.getroute()[0];  
                send_to = sending_order.getroute()[1];  
                received_from = sending_order.getroute()[size-1];
                if(rank == current_rank){
                    best = pop.getbest();
                    best_route = best.getroute();
                    best_size = best.getncity();
                    MPI_Send(best_route.data(), best_size, MPI_INT, send_to, 0, MPI_COMM_WORLD);
                    MPI_Recv(best_route.data(), best_size, MPI_INT, received_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    best.updateroute(best_route);
                    pop.exchangebest(best);
                }

                //Sending for loop: the j_th one sends to the (j+1)-th and receives from the (j-1)-th
                for(int j = 1; j < size - 1; j++){
                    current_rank = sending_order.getroute()[j];  
                    send_to = sending_order.getroute()[j+1];  
                    received_from = sending_order.getroute()[j-1];
                
                    if(rank == current_rank){
                        best = pop.getbest();
                        best_route = best.getroute();
                        best_size = best.getncity();
                        MPI_Send(best_route.data(), best_size, MPI_INT, send_to, 0, MPI_COMM_WORLD);
                        MPI_Recv(best_route.data(), best_size, MPI_INT, received_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        best.updateroute(best_route);
                        pop.exchangebest(best);
        
                    }

                }
        
                sending_order.reshuffle(rnd_shared);    //reshuffle sending order
                pop.sorting();
            }
        }

        Length.close();

    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    //Population of bests

    best = pop.getbest();
    best_route = best.getroute();
    best_size = best.getncity();

    // Declare send and receive buffers for the gathering operation
    vector<int> send_buffer(best_size);
    vector<int> receive_buffer(best_size * size);

    // Copy the information of the best individual into the send buffer
    copy(best_route.begin(), best_route.end(), send_buffer.begin());

    // Use MPI_Gather() to collect the information of the best individuals from all nodes
    MPI_Gather(send_buffer.data(), best_size, MPI_INT, receive_buffer.data(), best_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {    // Root process
        // Population to store the best individuals from each node
        Population bestpop(size, rnd_shared);
        for (int i = 0; i < size; i++) {
            best.updateroute(vector<int>(receive_buffer.begin() + (i * best_size), receive_buffer.begin() + ((i + 1) * best_size)));
            bestpop.insert(best, i);
        }

        bestpop.sorting();
        bestpop.print_bestroute();
    }
    
    MPI_Finalize();

    return 0;
}



