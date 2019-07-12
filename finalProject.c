#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <allegro5/allegro.h>

void update (int, int*, int, int, int*);

//row/col 9 = 9 procs
//row/col 6 = 4 procs
//row/col 12 = 16 procs
//row/col 33 = 121 procs

#define ROW 12
#define COL 12
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
int main(int argv, char** argc)
{

    MPI_Init(&argv, &argc);
    int dimension = sqrt (COL * ROW / 9); 
    int localRow = 3;
    int localCol = 3;
    int localSum = 0;    
    int nbrsSize = 3;
    int rank, dims[2] = {dimension, dimension}, nbrs[nbrsSize+1], periods  [2] = {0, 0}, reorder = 0;
    int* inbuf = (int*)malloc(sizeof(int)* nbrsSize*4);
    int* outbuf = (int*)malloc(sizeof(int)* nbrsSize*4);
    int* sendToMaster = (int*)malloc(sizeof(int)* (localRow*localCol+1));
    int* world = (int*)malloc(sizeof(int)* ROW * COL);
    int centerProcess = ROW % 2 == 0 ? dimension * dimension / 2 + ROW / 6 : dimension * dimension / 2 ;    
    int* localMatrix = (int*) malloc (sizeof(int)* localRow * localCol);
    MPI_Comm greed;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &greed);
    MPI_Comm_rank(greed, &rank);
    MPI_Cart_shift(greed, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(greed, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
    int centerCell = ROW % 2 == 1 ? 4 : 0;
    int decompositionRow = ROW/3;
    int decompositionCol = ROW/3;
    MPI_Request reqs[4];
    MPI_Status status;

    //needed for the corrispondence between inbuf and localMatrices    
    int nbrsRightShift = 7;
    int nbrsLeftShift = 6;
    int shiftCoefficient = 0; 
    
    //define master and initialize matrices (center or another external process)
    if (rank == centerProcess)
    {
        for (int i = 0; i< ROW * COL; i++)
            world [i] = 0;

//DISPLAY ALLEGRO
/*
        ALLEGRO_DISPLAY *display = NULL;

   if(!al_init()) {
      fprintf(stderr, "failed to initialize allegro!\n");
      return -1;
   }

   display = al_create_display(640, 480);
   if(!display) {
      fprintf(stderr, "failed to create display!\n");
      return -1;
   }

   al_clear_to_color(al_map_rgb(0,0,0));
   
   al_flip_display();
    
    
    */
    }    
    for (int i = 0; i <localRow * localCol; i++)
    {
        inbuf [i] = 0;
        outbuf [i] = 0;
        localMatrix [i] = 0;    
        sendToMaster[i] = 0;    
    }

    sendToMaster[localRow * localCol] = 0;
        
    //define the matrix split
    MPI_Datatype matrixTransposed;
    MPI_Datatype vector;

    //localRow block of length localCol with COL stride
    MPI_Type_vector(localRow, localCol, COL, MPI_INT, &vector);
    
    //starting from the previous one create a new type that start from 0 and goes to 0 + sizeof(int)     
    MPI_Type_create_resized( vector, 0, sizeof(int), &matrixTransposed);

    MPI_Type_commit(&matrixTransposed);
    //MPI_Type_commit(&vector);

    int disps[decompositionRow*decompositionCol];
    int counts[decompositionRow*decompositionCol];
    
    for (int i=0; i<decompositionRow; i++) {
        for (int j=0; j<decompositionCol; j++) {
            disps[i*decompositionCol+j] = i*COL*localRow+j*localCol;
            counts [i*decompositionCol+j] = 1;
        }
    }

    //send the submatrices to the processes
    MPI_Scatterv(world, counts, disps, matrixTransposed, localMatrix, localRow * localCol, MPI_INT, centerProcess, MPI_COMM_WORLD);

    //stampa di prova
    printf("rank %d Local Matrix:\n", rank);
            for (int i=0; i<localRow; i++) {
                for (int j=0; j<localCol; j++) {
                    printf("%d ",(int)localMatrix[i*localCol+j]);
                }
                printf("\n");
            }
            printf("\n");    

    //needed to synchronized the prints
    MPI_Barrier(greed);
    
    for(int u = 0; u<150; u  ++){
      
    //MASTER
    if (rank == centerProcess){

         //value of centerCell ++; localMatrix[4]++;
         localMatrix[centerCell]++;
        //update recursive function
        for (int i = 0; i<localRow*localCol; i++) 
            if (localMatrix[i] > 3)
                update(i, localMatrix, 1, 1, outbuf);        

        //send my state to nbrs
        for (int i = 0; i<nbrsSize*4; i+=3)
            MPI_Isend(&outbuf[i], 3, MPI_INT, nbrs[((i+3) /3)-1],1, greed, &reqs[((i+3) /3)-1]);
        
        for (int i = 0; i<nbrsSize*4;i++)
            outbuf[i] = 0;

        // //receive nbrs states 
        for (int i = 0; i<12;i+=3)
            MPI_Recv(&inbuf[i],3, MPI_INT, nbrs[((i+3) /3)-1], 1, greed, &status);

    
        for (int i = 0; i<nbrsSize+1; i++)
            MPI_Wait (&reqs[i], &status);

        //prepare next iteration and send 
        //corrispondenze:
        //UP localMatrix [i] <-> inbuf [i] ; 0-3  ++1        
        for (int i = 0; i<3; i++) 
            localMatrix[i]+=inbuf[i];
        
        //DOWN localMatrix [i] <-> inbuf [i-3] ; 6-9  ++1
        for (int i =6; i< 9; i++)
            localMatrix[i]+=inbuf[i-3];

    
        //LEFT localMatrix [i] <-> inbuf [i+nbrsLeftShift-shiftCoefficient] ; 0-9 ++3; x=6, y=0 ++2
        for (int i = 0; i< 9; i=i+3)
        {
            localMatrix[i]+= inbuf[i+nbrsLeftShift - shiftCoefficient];        
            shiftCoefficient+=2 ;
        }
        
        shiftCoefficient = 0;
             
        //RIGHT localMatrix [i] <-> inbuf [i+nbrsRightShift-shiftCoefficient] ; 2-9 ++3; x=7, y=0 ++2
        for (int i = 2; i<9; i= i+3)
        {
            localMatrix[i]+= inbuf[i+nbrsRightShift - shiftCoefficient];
            shiftCoefficient+=2;        
        }

        shiftCoefficient = 0;
        
        int printOk = 1;
        for (int i = 0; i<ROW*COL; i++)
            if (world[i] > 3)
                printOk = 0;
        
             
        //if no messagge is flying print the world WORK IN PROGRESS  /TODO CHECK IF NEEDED, MAYBE NOT 
        if (printOk == 1){        
        printf ("giro %d\n", u);
        for (int i = 0; i< ROW; i++){
            for (int j = 0; j<COL; j++){
                    printf ("%d  ", world[i*COL+j]);} printf ("\n");
        }
        printf ("\n");
       }

            
    }

    //local operations 
    if (rank != centerProcess){
    
       // localMatrix[4]++;
        //update recursive function
          for (int i = 0; i<9; i++) 
            if (localMatrix[i] > 3)
                update(i, localMatrix, 1, 1, outbuf);        

        //send my state to nbrs
        for (int i = 0; i<nbrsSize*4; i+=3)
            MPI_Isend(&outbuf[i], 3, MPI_INT, nbrs[((i+3) /3)-1],1, greed, &reqs[((i+3) /3)-1]);
        
        for (int i = 0; i<nbrsSize*4;i++)
            outbuf[i] = 0;
        //receive
        for (int i = 0; i<nbrsSize*4;i+=3)
            MPI_Recv(&inbuf[i],3, MPI_INT, nbrs[((i+3) /3)-1], 1, greed, &status);
        
        for (int i = 0; i<4; i++)
            MPI_Wait (&reqs[i], &status);
        //corrispondenze:
        //UP localMatrix [i] <-> inbuf [i] ; 0-3  ++1        
        for (int i = 0; i<3; i++) 
            localMatrix[i]+=inbuf[i];
        
        //DOWN localMatrix [i] <-> inbuf [i-3] ; 6-9  ++1
        for (int i =6; i< 9; i++)
            localMatrix[i]+=inbuf[i-3];

    
        //LEFT localMatrix [i] <-> inbuf [i+nbrsLeftShift-shiftCoefficient] ; 0-9 ++3; x=6, y=0 ++2
        for (int i = 0; i< 9; i=i+3)
        {
            localMatrix[i]+= inbuf[i+nbrsLeftShift - shiftCoefficient];        
            shiftCoefficient+=2 ;
        }
        
        shiftCoefficient = 0;
             
        //RIGHT localMatrix [i] <-> inbuf [i+nbrsRightShift-shiftCoefficient] ; 2-9 ++3; x=7, y=0 ++2
        for (int i = 2; i<9; i= i+3)
        {
            localMatrix[i]+= inbuf[i+nbrsRightShift - shiftCoefficient];
            shiftCoefficient+=2;        
        }

        shiftCoefficient = 0;

        //prepareNextSend

    
        }

       //sendToMaster the actual state + sum where sum is the next send (maybe gather) //TODO IL SUM
       MPI_Gatherv(localMatrix, localRow * localCol, MPI_INT, world, counts, disps, matrixTransposed, centerProcess, MPI_COMM_WORLD);   

    }

    MPI_Finalize();       

}

void update (int cell, int* localMatrix, int updateNum, int center, int* outbuf)
{

    
    if (localMatrix[cell] <3)
        localMatrix[cell] += updateNum;

    else
    {
        localMatrix[cell] = 0;
        if (cell == 0){
            outbuf[0]+=1;
            outbuf[6]+=1;
            update(1, localMatrix, 1, 0, outbuf);
            update(3, localMatrix, 1, 0, outbuf);
        }
        
        if (cell == 1){
            outbuf[1]+=1;
            update(0, localMatrix, 1, 0, outbuf);
            update(2, localMatrix, 1, 0, outbuf);
            update(4, localMatrix, 1, 0, outbuf);        
            }   

        if (cell == 2){
            outbuf[2]+=1;
            outbuf[9]+=1;
            update(1, localMatrix, 1, 0, outbuf);
            update(5, localMatrix, 1, 0, outbuf);  
          }   
    
        if (cell == 3){
            outbuf[7]+=1;
            update(0, localMatrix, 1, 0, outbuf);
            update(4, localMatrix, 1, 0, outbuf);
            update(6, localMatrix, 1, 0, outbuf);        
        }
        
        if (cell == 4){
            update(1, localMatrix, 1, 0, outbuf);
            update(3, localMatrix, 1, 0, outbuf);
            update(5, localMatrix, 1, 0, outbuf);
            update(7, localMatrix, 1, 0, outbuf);        
        }
    
        if (cell == 5){
            outbuf[10]+=1;
            update(2, localMatrix, 1, 0, outbuf);
            update(4, localMatrix, 1, 0, outbuf);
            update(8, localMatrix, 1, 0, outbuf);
        }
    
        if (cell == 6){
            outbuf[3]+=1;
            outbuf[8]+=1;
            update(3, localMatrix, 1, 0, outbuf);
            update(7, localMatrix, 1, 0, outbuf);
        }
 
        if (cell == 7){
            outbuf[4]+=1;
            update(4, localMatrix, 1, 0, outbuf);
            update(6, localMatrix, 1, 0, outbuf);
            update(8, localMatrix, 1, 0, outbuf);
        }

        if (cell == 8){
        outbuf[5]+=1;
        outbuf[11]+=1;
        update(5, localMatrix, 1, 0, outbuf);
        update(7, localMatrix, 1, 0, outbuf);
         }              
    }
}


