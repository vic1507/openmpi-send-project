#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>

void update (int, int*, int, int, int*, int, int, int);
//void update (int, int*, int, int, int*);

//row/col 9 = 9 procs
//row/col 6 = 4 procs
//row/col 12 = 16 procs
//row/col 27 = 81 procs
//row/col 24 = 64 procs
//row/col 54 = 324 procs //TOO MUCH! NEED BIGGER DECOMPOSITION!

#define ROW 54
#define COL 54
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int numProcs;
    MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    int localRow = 27;//sqrt(COL * ROW / numProcs); //TODO DEVE FARLO IN AUTOMATICO
    int localCol = 27; //sqrt(COL * ROW / numProcs); //TODO DEVE FARLO IN AUTOMATICO
    int dimension = sqrt (COL * ROW / (localRow * localCol)); 
    int localSum = 0;    
    //TODO IL DIMS DEVE FARLO IN AUTOMATICO
    int rank, dims[2] = {2,2/*dimension, dimension*/}, nbrs[4], periods  [2] = {0, 0}, reorder = 0;
    int* inbuf = (int*)malloc(sizeof(int)* (localRow * 2 + localCol *2));
    int* outbuf = (int*)malloc(sizeof(int)* (localRow * 2 + localCol *2));
    int* sendToMaster = (int*)malloc(sizeof(int)* (localRow*localCol));
    int* world = (int*)malloc(sizeof(int)* ROW * COL);
    
    //TODO DEFINIRE IL PROCESSO CENTRALE
    int centerProcess = 3;//numProcs % 2 == 0 ? numProcs / 2 +1 : numProcs / 2;  
    
    int* localMatrix = (int*) malloc (sizeof(int)* localRow * localCol);
    MPI_Comm greed;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &greed);
    MPI_Comm_rank(greed, &rank);
    MPI_Cart_shift(greed, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(greed, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
        
    //TODO DEFINIRE LA CELLA CENTRALE    
    int centerCell = 0;//numProcs % 2 == 0 ? 0 : localRow * localCol / 2;

    int decompositionRow = ROW/localRow;
    int decompositionCol = COL/localCol;
    MPI_Request reqs[4];
    int sendSize = localRow;
    MPI_Status status;

    //needed for the corrispondence between inbuf and localMatrices    
    int nbrsRightShift = localRow + localCol + 1;
    int nbrsLeftShift = localCol * 2;
    int shiftCoefficient = 0; 
    
    //define master and initialize matrices (center or another external process)
    if (rank == centerProcess)
    {
        for (int i = 0; i< ROW * COL; i++)
            world [i] = 0;

//DISPLAY ALLEGRO

        ALLEGRO_DISPLAY *display = NULL;

   if(!al_init()) {
      fprintf(stderr, "failed to initialize allegro!\n");
      return -1;
   }

   display = al_create_display(640, 640);
   if(!display) {
      fprintf(stderr, "failed to create display!\n");
      return -1;
   }

       al_clear_to_color(al_map_rgb(0,0,0));
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

    
    //stampa di prova//
   // printf("rank %d Local Matrix:\n", rank);
   //         for (int i=0; i<localRow; i++) {
   //             for (int j=0; j<localCol; j++) {
//                  al_draw_filled_circle (50.0 , 50.0 , 10.0, al_map_rgb(255.0, 255.0, 255.0));
    //                printf("%d ",(int)localMatrix[i*localCol+j]);
    //            }
      //          printf("\n");
        //    }
        //   printf("\n");    
    //needed to synchronized the decomposition and the send
    MPI_Barrier(greed);
    
    while (1){
      
    //MASTER
    if (rank == centerProcess){

         //value of centerCell ++; localMatrix[4]++;
         localMatrix[centerCell]++;
        //update recursive function
        for (int i = 0; i<localRow*localCol; i++) 
            if (localMatrix[i] > 3)
                update(i, localMatrix, 1, 1, outbuf, localCol, localRow, 0);        

         int waitArray[4] = {0,0,0,0};

        //send my state to nbrs
        if (localRow != localCol)
        {
            for (int i = 0; i< localCol * 2; i+= localCol){
                if (nbrs[((i+localCol) /localCol)-1] >= 0 && nbrs[((i+localCol) /localCol)-1]<numProcs){
                waitArray[((i+localCol) /localCol)-1] = 1;
                 MPI_Isend(&outbuf[i], localCol, MPI_INT, nbrs[((i+localCol) /localCol)-1],1, greed, &reqs[((i+localCol) /localCol)-1]);
                }
            }
           for (int i = localCol * 2; i < localRow * 2 + localCol * 2; i+=localRow){
            if (nbrs[(i - localCol * 2 + localRow) / localRow + 1] >= 0 && nbrs[(i - localCol * 2 + localRow) / localRow + 1]<numProcs){
                waitArray[(i - localCol * 2 + localRow) / localRow + 1] = 1;
                 MPI_Isend(&outbuf[i], localRow, MPI_INT, nbrs[(i - localCol * 2 + localRow) / localRow + 1],1, greed, &reqs[(i - localCol * 2 + localRow) / localRow + 1]); }
                }          
            }
        else
            for (int i = 0; i<localRow * 2 + localCol *2; i+=localRow){   
                 MPI_Isend(&outbuf[i], localRow, MPI_INT, nbrs[((i+localRow) /localRow)-1],1, greed, &reqs[((i+localRow) /localRow)-1]);
                 waitArray[((i+localRow) /localRow)-1] = 1;                
              }


        for (int i = 0; i<localCol * 2 + localRow * 2;i++)
            outbuf[i] = 0;

    
        // //receive nbrs states
        if(localRow != localCol){
            for (int i = 0; i< localCol *2 ; i+= localCol){
            if (nbrs[((i+localCol) /localCol)-1] >= 0 && nbrs[((i+localCol) /localCol)-1]<numProcs)
             MPI_Recv(&inbuf[i], localCol, MPI_INT, nbrs[((i+localCol) /localCol)-1], 1, greed, &status);
            }
            for (int i = localCol*2; i<localRow * 2 +  localCol*2; i+=localRow){
                if (nbrs[(i - localCol * 2 + localRow) / localRow + 1] >= 0 && nbrs[(i - localCol * 2 + localRow) / localRow + 1]<numProcs)
                     MPI_Recv(&inbuf[i], localRow, MPI_INT, nbrs[(i - localCol * 2 + localRow) / localRow + 1], 1, greed, &status);      }
        }
        else 
            for (int i = 0; i<localRow * 2 + localCol *2;i+=localRow)
                MPI_Recv(&inbuf[i], localRow, MPI_INT, nbrs[((i+localRow) /localRow)-1], 1, greed, &status);

    
    for (int i = 0; i<4; i++){if (waitArray[i] ==1)
            MPI_Wait (&reqs[i], &status);}

        //prepare next iteration and send 
        //corrispondenze:
        
        //UP localMatrix [i] <-> inbuf [i] ; 0-3  ++1        
        for (int i = 0; i<localCol; i++) 
            localMatrix[i]+=inbuf[i];
        
        //DOWN localMatrix [i] <-> inbuf [i-3] ; 6-9  ++1
        for (int i = localCol*(localRow-1); i< localCol*localRow; i++){
        if(i- localCol* (localRow - 2) <0)
            printf("sono la stella\n");
                  localMatrix[i]+=inbuf[i- localCol* (localRow - 2)];
}
    
        //LEFT localMatrix [i] <-> inbuf [i+nbrsLeftShift-shiftCoefficient] ; 0-9 ++3; x=6, y=0 ++2
        for (int i = 0; i<= localRow * localCol - localCol; i+=localCol)
        {
            localMatrix[i]+= inbuf[i+nbrsLeftShift - shiftCoefficient];        
            shiftCoefficient+= localCol -1 ;
            if (i+nbrsLeftShift - shiftCoefficient <0)
                printf("mi chiamo virgola\n");
        }
        
        shiftCoefficient = 0;
             
        //RIGHT localMatrix [i] <-> inbuf [i+nbrsRightShift-shiftCoefficient] ; 2-9 ++3; x=7, y=0 ++2
        for (int i = localCol-1; i<localRow * localCol; i+=localCol)
        {

            if (i+nbrsRightShift - shiftCoefficient <0)
                printf("sono un gattino\n");
            localMatrix[i]+= inbuf[i+nbrsRightShift - shiftCoefficient];
            shiftCoefficient+=localCol - 1;        
        }        

          shiftCoefficient = 0;
        
        int printOk = 1;
        for (int i = 0; i<ROW*COL; i++)
            if (world[i] > 3)
                printOk = 0;
            
    
         if (printOk == 1){  
        float c1, c2, c3 = 0.0;
        for (int i = 0; i<ROW; i++){
                for (int j = 0; j<COL; j++){
            if (world[i*COL+j] == 0){ c1 = 255.0; c2 = 229.0; c3 = 204.0;  c1 = 0.0; c2 = 0.0; c3 = 0.0;}
            else if (world[i*COL+j] == 1){c1 = 204.0; c2 = 102.0; c3 = 0.0;}
            else if (world[i*COL+j] == 2){c1 = 153.0; c2 = 76.0; c3 = 0.0;}
            else if (world[i*COL+j] == 3){c1 = 102.0; c2= 51.0; c3 = 0.0;}
            al_draw_filled_circle (j * 10.0 + 20  , i * 10.0  + 10.0 , 5.0, al_map_rgb(c1, c2, c3));
                   
            }
        }
        al_flip_display();/* 
        for (int i = 0; i< ROW; i++){
            for (int j = 0; j<COL; j++){ 
                    printf ("%d  ", world[i*COL+j]);
            } printf ("\n");
          }
        printf ("\n");*/
       }

}

    //local operations 
    if (rank != centerProcess){
       // localMatrix[4]++;
        //update recursive function
          for (int i = 0; i<localRow * localCol; i++) 
            if (localMatrix[i] > 3)
                update(i, localMatrix, 1, 1, outbuf, localCol, localRow, 0);        

          int waitArray[4] = {0,0,0,0};


      //send my state to nbrs
        if (localRow != localCol)
        {
            for (int i = 0; i< localCol * 2; i+=localCol){
                if (nbrs[((i+localCol) /localCol)-1] >= 0 && nbrs[((i+localCol) /localCol)-1] < numProcs){
                 MPI_Isend(&outbuf[i], localCol, MPI_INT, nbrs[((i+localCol) /localCol)-1],1, greed, &reqs[((i+localCol) /localCol)-1]);
                waitArray[((i+localCol) /localCol)-1] = 1;
                }
            }
            for (int i = localCol * 2; i < localRow * 2 + localCol * 2; i+=localRow){
            if (nbrs[(i - localCol * 2 + localRow) / localRow + 1] >= 0 && nbrs[(i - localCol * 2 + localRow) / localRow + 1]<numProcs){
                waitArray[((i+localCol) /localCol)-1] = 1;
                 MPI_Isend(&outbuf[i], localRow, MPI_INT, nbrs[(i - localCol * 2 + localRow) / localRow + 1],1, greed, &reqs[(i - localCol * 2 + localRow) / localRow + 1]); }}          
            }
        else
            for (int i = 0; i<localRow * 2 + localCol *2; i+=localRow){
                MPI_Isend(&outbuf[i], localRow, MPI_INT, nbrs[((i+localRow) /localRow)-1],1, greed, &reqs[((i+localRow) /localRow)-1]);
                   waitArray[((i+localRow) /localRow)-1] = 1;        
        }
        
        for (int i = 0; i<localRow*2 + localCol *2 ;i++)
            outbuf[i] = 0;
        
        //receive
        if(localRow != localCol){
            for (int i = 0; i< localCol *2 ; i+= localCol){
                if (nbrs[((i+localCol) /localCol)-1] >= 0 && nbrs[((i+localCol) /localCol)-1]<numProcs)
                    MPI_Recv(&inbuf[i], localCol, MPI_INT, nbrs[((i+localCol) /localCol)-1], 1, greed, &status);
                }
            for (int i = localCol*2; i<localRow * 2 +  localCol*2; i+=localRow){
                if (nbrs[(i - localCol * 2 + localRow) / localRow + 1] >= 0 && nbrs[(i - localCol * 2 + localRow) / localRow + 1]<numProcs){
                    MPI_Recv(&inbuf[i], localRow, MPI_INT, nbrs[(i - localCol * 2 + localRow) / localRow + 1], 1, greed, &status);}
           }        
        }
        else 
            for (int i = 0; i<localRow * 2 + localCol *2;i+=localRow)
                MPI_Recv(&inbuf[i], localRow, MPI_INT, nbrs[((i+localRow) /localRow)-1], 1, greed, &status);
        
        for (int i = 0; i<4; i++){
          if (waitArray[i] ==1)
                MPI_Wait (&reqs[i], &status);
        }

        //corrispondenze:
        //UP localMatrix [i] <-> inbuf [i] ; 0-3  ++1        
        for (int i = 0; i<localCol; i++) 
            localMatrix[i]+=inbuf[i];
        
        //DOWN localMatrix [i] <-> inbuf [i-3] ; 6-9  ++1
        for (int i = localCol*(localRow-1); i< localCol*localRow; i++)
            localMatrix[i]+=inbuf[i- localCol* (localRow - 2)];

    
        //LEFT localMatrix [i] <-> inbuf [i+nbrsLeftShift-shiftCoefficient] ; 0-9 ++3; x=6, y=0 ++2
        for (int i = 0; i<= localRow * localCol - localCol; i+=localCol)
        {
            localMatrix[i]+= inbuf[i+nbrsLeftShift - shiftCoefficient];        
            shiftCoefficient+= localCol -1 ;
        }
        
        shiftCoefficient = 0;
             
        //RIGHT localMatrix [i] <-> inbuf [i+nbrsRightShift-shiftCoefficient] ; 2-9 ++3; x=7, y=0 ++2
        for (int i = localCol-1; i<localRow * localCol; i+=localCol)
        {

            localMatrix[i]+= inbuf[i+nbrsRightShift - shiftCoefficient];
            shiftCoefficient+=localCol - 1;        
        }        

          shiftCoefficient = 0;
    
        }

       //sendToMaster the actual state + sum where sum is the next send (maybe gather) 
       MPI_Gatherv(localMatrix, localRow * localCol, MPI_INT, world, counts, disps, matrixTransposed, centerProcess, MPI_COMM_WORLD);   

    }

    MPI_Finalize();       

}


void update (int cell, int* localMatrix, int updateNum, int center, int* outbuf, int sizeCol, int sizeRow, int lap)
{
    
    

    if (localMatrix[cell] <3)
        localMatrix[cell] += updateNum;

    else
    {

        int sendLeft = 0;
        int sendRight = 0;
        int sendUp = 0;
        int sendDown = 0;
        localMatrix[cell] = 0;
        
        int start = sizeCol-1;
        int stride = 1;
        while (start <= sizeRow*sizeCol-1){
           if (cell == start){
                    sendLeft = 1;
                break;
            }
            start+=sizeCol;     
            stride++;        
        }   
        
        //DOWN IN LOCAL, UP IN GENERAL   IT'S OK!
                //pirma riga                   //prima colonna tranne prima e ultima cella                                ultima colonna tranne prima e ultima cella
        if (cell >= 0 && cell < sizeCol || (cell % sizeCol == 0 && cell < sizeCol * sizeRow - sizeCol && cell > 0) || (sendLeft == 1 && cell > sizeCol-1 && cell<sizeRow*sizeCol-1))
        {
            sendDown = 1;   
            if (cell >= 0 && cell < sizeCol )
                outbuf[cell]+=1;         
            update (cell+sizeCol, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
        }        

        //UP IN LOCAL, DOWN IN GENERAL  IT'S OK!
        if (cell < sizeRow * sizeCol && cell>= sizeRow*sizeCol-sizeCol || (cell % sizeCol == 0 && cell < sizeCol * sizeRow - sizeCol && cell > 0) || (sendLeft == 1 && cell > sizeCol-1 && cell<sizeRow*sizeCol-1) )
        {

            sendUp = 1;
            if (cell < sizeRow * sizeCol && cell>= sizeRow*sizeCol-sizeCol)
                outbuf[cell-  (sizeCol * (sizeRow -2)) ]+=1; 
            
            update (cell-sizeCol, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
        }

        //RIGHT IN LOCAL, LEFT IN GENERAL IT'S OK!    
        if (cell % sizeCol == 0 || (cell > 0 && cell< sizeCol-1) || (cell<sizeRow*sizeCol-1 && cell> sizeRow*sizeCol-sizeCol))
        {
            if (cell%sizeCol == 0){
                outbuf[sizeCol*2 + (cell/sizeCol)]+=1;
            }
            update (cell+1, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
            
            sendRight = 1;
            //send right
        }


    //SINISTRA LOCALE, DESTRA GENERALE IT'S OK!!!
        if (sendLeft == 1 || (cell <sizeRow * sizeCol -1 && cell > sizeRow * sizeCol - sizeCol) || (cell>0 && cell<sizeCol-1)){
            //sendLeft
            if (sendLeft == 1)
                outbuf[sizeCol * 2 + sizeRow + cell - stride * (sizeCol-1)]+=1;
            update (cell-1, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
        }

        if (sendLeft + sendRight + sendUp + sendDown == 0){

              //send everywhere
            update (cell-1, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
            update (cell+1, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
            update (cell+sizeCol, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
            update (cell-sizeCol, localMatrix, 1, 0, outbuf, sizeCol, sizeRow, lap);
        }
    }
}


