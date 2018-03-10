#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include"omp.h"

/*
 Global parameters
 */
#define N 194
#define MAX 0x7fffffff
#define MASTER_TAG 0
#define WORKER_TAG 1

/*
 Parameters for GA
 */
#define P_CROSS 0.6
#define P_MUT 0.1
#define MAX_GEN 2000

/*
 Parameters for AC
 */
#define NUM_ANT 2000
#define RC_ITER 1000
#define IN_PHE 1
#define tao_max 20
#define tao_min 0.0001
#define mig_rate 0.1

/*
 Structure to represent cities
 */
struct coordinate{
    int city;
    float x;
    float y;
}coords[N];
double graph[N][N];
/*
 global loop variables
 */
int i,j,k;
int aver_popu = NUM_ANT/4;
int aver_ant = NUM_ANT/4;
/*
 variables used in AC
 */
double phe[N][N];
double add[N][N];
double yita[N][N];
int vis[NUM_ANT][N];
int map[NUM_ANT][N];
double solution[NUM_ANT];
int bestway[N];
double bestsolution=MAX;
int NcMax;
double alpha;
double beta;
double rou;
double Q;
clock_t start;
/*
 global methods
 */
void Inputcoords(FILE *fp);
void CreateGraph();
double Distance(int *p);
void Result();

/*
 GA methods
 */
void init_gen(int chrom[NUM_ANT][N]);
void Choice(int chrom[NUM_ANT][N]); //choose based on fitness
void Cross(int chrom[NUM_ANT][N]); //random cross inside chrome
void Mutation(int chrom[NUM_ANT][N]); //random mutation inside chrome
void Reverse(int chrom[NUM_ANT][N]); //reverse to ensure

/*
 AC methods
 */
void init_ac();

/*
 main function
 */
int main(){
    /*
     Read map and generate matrix
     */
    FILE *fp;
    fp = fopen("map_qatar.txt","r+");
    Inputcoords(fp);
    CreateGraph();
    fclose(fp);
    
    
    
    /* MPI related variables*/
    int subpros,rank,provided;
    
    
    int chrom[NUM_ANT][N];
    int NcChunk = 200;
    int GenChunk = 2000;
    
    
    MPI_Status status;
    
    /* initialize MPI procedure */
    MPI_Init_thread(NULL,NULL,MPI_THREAD_FUNNELED,&provided);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    subpros = 3;
    
    int subchrom[aver_popu][N];
    double subbestsolution = MAX;
    int subbestway[N];
    
    if (rank == 0) {
        
        
        /*
         initial variables
         */
        memset(solution,0.0,sizeof(solution));
        memset(bestway, 0, sizeof(bestway));
        
	start = MPI_Wtime();
        /*
         GA part
         */
        
        memset(chrom, 0, sizeof(chrom));
        init_gen(chrom);
        
        for(i=0;i<aver_popu;i++){
            for(j=0;j<N;j++){
                subchrom[i][j] = chrom[i][j];
            }
        }
	int j;
#pragma omp parallel for
        for(j=0;j<NUM_ANT;j++){
            solution[j] = Distance(chrom[j]);
        }
        #pragma omp critical


        for(j=0;j<NUM_ANT;j++){
            if(solution[j] < bestsolution){
                bestsolution = solution[j];
                for(i=0;i<N;i++){
                    bestway[i] = chrom[j][i];
                }
            }
        }

        for(i=1;i<=subpros;i++){
            MPI_Send(&chrom[i*aver_popu][0],aver_popu*N,MPI_INT,i,MASTER_TAG,MPI_COMM_WORLD);
        }
    }

    if(rank!=0){
        //MPI_Recv(&bestsolution, 1, MPI_DOUBLE, 0, MASTER_TAG, MPI_COMM_WORLD, &status);
        memset(subchrom, 0, sizeof(subchrom));
        MPI_Recv(&subchrom[0][0], aver_popu*N,MPI_INT,0,MASTER_TAG, MPI_COMM_WORLD,&status);
    }
    
    MPI_Bcast(&bestsolution, 1, MPI_DOUBLE, MASTER_TAG,MPI_COMM_WORLD);
    
    
    double subsolution[aver_popu];
    subbestsolution = bestsolution;
    int num_gen = 0;
    int mig_flag = 0;
    
    while(num_gen<MAX_GEN/GenChunk){
        
        mig_flag = 0;
        int gen = 0;
        while(gen<GenChunk){
            //        printf("%d\n",num_gen);
            
            //initial solution each generation
            memset(subsolution, 0, sizeof(subsolution));
            
            //do GA operation on chrom
            Choice(subchrom);
            Cross(subchrom);
            Mutation(subchrom);
            Reverse(subchrom);
	    int j;
#pragma omp parallel for
            for(j=0;j<aver_popu;j++){
                subsolution[j] = Distance(subchrom[j]);
            }
#pragma omp critical
            for(j=0;j<aver_popu;j++){
                if(subsolution[j] < subbestsolution){
                    subbestsolution = subsolution[j];
                    for(k=0;k<N;k++){
                        subbestway[k] = subchrom[j][k];
                    }
                }
            }

            gen+=1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
	
        if(rank==0){
            double mrate = ((double)rand())/(RAND_MAX+1.0);
            if(mrate<mig_rate)
                mig_flag = 1;
        }
        MPI_Bcast(&mig_flag, 1, MPI_INT, MASTER_TAG, MPI_COMM_WORLD);
        
        if(mig_flag==1){
          printf("!!!!!!!!!!!!!!!!!!!!!!\n");
            MPI_Request req;
            if(rank==0){
                MPI_Send(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 1, num_gen+10, MPI_COMM_WORLD);
                MPI_Recv(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 3, num_gen+13, MPI_COMM_WORLD, &status);
            }
            if(rank==1){
                MPI_Send(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 2, num_gen+11, MPI_COMM_WORLD);
                MPI_Recv(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 0, num_gen+10, MPI_COMM_WORLD, &status);
            }
            if(rank==2){
                MPI_Send(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 3, num_gen+12, MPI_COMM_WORLD);
                MPI_Recv(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 1, num_gen+11, MPI_COMM_WORLD, &status);
            }
            if(rank==3){
                MPI_Send(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 0, num_gen+13, MPI_COMM_WORLD);
                MPI_Recv(&subchrom[0][0],(aver_popu/2)*N,MPI_INT, 2, num_gen+12, MPI_COMM_WORLD, &status);
            }
        }

        
        num_gen+=1;
    }
    

    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Send back sub chrom to master
    if(rank!=0){
        MPI_Send(&subchrom[0][0],aver_popu*N,MPI_INT, 0, MASTER_TAG, MPI_COMM_WORLD);
        MPI_Send(&subbestsolution,1,MPI_DOUBLE, 0, MASTER_TAG, MPI_COMM_WORLD);
    }
    
    int subvis[aver_ant][N];
    int submap[aver_ant][N];
    double subphe[N][N];
    double subyita[N][N];
    double subadd[N][N];
    
    
    if(rank==0){
        for(i=0;i<aver_popu;i++){
            for(j=0;j<N;j++){
                chrom[i][j] = subchrom[i][j];
            }
        }
        
        //
        int recchrom[aver_popu][N];
        double recsubbest;
        
        //        printf("%f\n",bestsolution);
        for(i=1;i<=subpros;i++){
            memset(recchrom,0,sizeof(recchrom));
            
            MPI_Recv(&recchrom[0][0], aver_popu*N, MPI_INT, i, MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&recsubbest, 1, MPI_DOUBLE, i, MASTER_TAG, MPI_COMM_WORLD, &status);
            if(recsubbest<bestsolution){
                bestsolution = recsubbest;
            }
            for(j=0;j<aver_popu;j++){
                for(k=0;k<N;k++){
                    chrom[aver_popu*i+j][k] = recchrom[j][k];
                }
            }
        }
    }
    //
    
    if(rank == 0){
        /*
         Now got final chrom generated by GA
         */
        
        
        /*
         start par ac
         */
        
        /*
         update AC matrix with GA solution
         */
        init_ac();
        //initial pheromone
        memset(yita,0,sizeof(yita));
        for(i=0; i<N; ++i){
            for(j=0; j<N; ++j){
                phe[i][j] = IN_PHE;
                if(i != j)
                    yita[i][j] = 1.0 / graph[i][j];
            }
        }
        
        memset(add,0,sizeof(add));
        
        for(k=0; k<NUM_ANT; k++){
            for(j=0; j<N-1; j++){
                add[ chrom[k][j] ][ chrom[k][j+1] ] += Q/solution[k];
            }
            add[ chrom[k][N-1] ][ chrom[k][0] ] += Q/solution[k];
        }
	int i,j;
#pragma omp parallel for collapse(2)
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                phe[i][j] = phe[i][j]*rou + add[i][j];
                if(phe[i][j] < tao_min)
                    phe[i][j] = tao_min;
                else if(phe[i][j] > tao_max)
                    phe[i][j] = tao_max;
            }
        }
#pragma omp critical
        /*
         Start MPI AC algorithm
         */
        
        for (i=1; i<=subpros; i++) {
            MPI_Send(&phe[0][0],N*N, MPI_DOUBLE, i, MASTER_TAG, MPI_COMM_WORLD);
            MPI_Send(&yita[0][0],N*N, MPI_DOUBLE, i, MASTER_TAG, MPI_COMM_WORLD);
            
        }
        
        /*
         begin ac part
         */
        
        memset(submap, -1, sizeof(submap));
        memset(subvis, 0, sizeof(subvis));
        memset(subphe, 0.0, sizeof(subphe));
        memset(subyita, 0.0, sizeof(subyita));
        
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                subphe[i][j] = phe[i][j];
                subyita[i][j] = yita[i][j];
            }
        }
    }
    
    
    if(rank!=0){
        
        
        //        printf("rank: %d, aver_ant: %d\n",rank,aver_ant);
        
        memset(submap, -1, sizeof(submap));
        memset(subvis, 0, sizeof(subvis));
        memset(subphe, 0.0, sizeof(subphe));
        memset(subyita, 0.0, sizeof(subyita));
        
        MPI_Recv(&subphe[0][0], N*N, MPI_DOUBLE, 0, MASTER_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&subyita[0][0], N*N, MPI_DOUBLE, 0, MASTER_TAG, MPI_COMM_WORLD, &status);
        
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    init_ac();
    int outer_Nc = 0;
    while(outer_Nc++ <= NcMax/NcChunk){
        int inner_Nc = 0;
        int s;
        double drand, pro, psum;
        memset(submap, -1, sizeof(submap));
        memset(subvis, 0, sizeof(subvis));
        
        while(inner_Nc < NcChunk){
            for(k=0;k<aver_ant;k++){
                submap[k][0] = (k+inner_Nc) % N;
                subvis[k][submap[k][0]] = 1;
            }

            s=1;
            
            while(s<N){
                for(k=0; k<aver_ant; k++){
                    psum = 0;
                    for(j=0; j<N; j++){
                        if(subvis[k][j] == 0){
                            psum += pow(subphe[ submap[k][s-1] ][j], alpha) * pow(subyita[ submap[k][s-1] ][j], beta);
                        }
                    }
                    drand = (double)(rand() % 5000);
                    drand /= 5000.0;
                    pro = 0;
                    
                    int y=0;
                    for(y=0; y<N; y++){

                        if(subvis[k][y] == 0){
                            pro += pow(subphe[ submap[k][s-1] ][y], alpha) * pow(subyita[ submap[k][s-1] ][y], beta) / psum;
                        }
                        if(pro > drand)
                            break;
                    }
                    subvis[k][y] = 1;
                    submap[k][s] = y;
                    
                }
                s++;
            }
//            if(rank==0){
//                for(i=0;i<aver_ant;i++){
//                    printf("%d -> ",submap[4][i]);
//                }
//            }
//            printf("\n");

            
            memset(subadd, 0.0, sizeof(subadd));
            memset(subsolution,0,sizeof(subsolution));
            int j;
            #pragma omp parallel for
            for(j=0;j<aver_ant;j++){
                subsolution[j] = Distance(submap[j]);
            }
#pragma omp critical
            for(j=0;j<aver_ant;j++){
                if(subsolution[j] < subbestsolution){
                    subbestsolution = subsolution[j];
                    for(k=0;k<N;k++){
                        subbestway[k] = submap[j][k];
                    }
                }
            }

            
//            if(subbestsolution<1000){
//                for(i=0;i<N;i++){
//                    printf("%d >",subbestway[i]);
//                }
//            }
//            printf("\n");
            
            for(k=0; k<aver_ant; k++){
                for(j=0; j<N-1; j++){
                    subadd[ submap[k][j] ][ submap[k][j+1] ] += Q/subsolution[k];
                }
                subadd[ submap[k][N-1] ][ submap[k][0] ] += Q/subsolution[k];
            }
	    int i;
#pragma omp parallel for collapse(2)
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    subphe[i][j] = subphe[i][j]*rou + subadd[i][j];
                    if(subphe[i][j] < tao_min)
                        subphe[i][j] = tao_min;
                    else if(subphe[i][j] > tao_max)
                        subphe[i][j] = tao_max;
                }
            }
#pragma omp critical
            memset(subvis,0,sizeof(subvis));
            memset(submap,-1,sizeof(submap));
            inner_Nc+=1;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank!=0){
            MPI_Send(&subphe[0][0], N*N, MPI_DOUBLE, 0, MASTER_TAG,MPI_COMM_WORLD);
        }
        if(rank==0){
            double recephe[N][N];
            double finalphe[N][N];
            
            memset(recephe, 0.0, sizeof(recephe));
            memset(finalphe, 0.0, sizeof(finalphe));
            
            for (i=0; i<N; i++) {
                for (j=0; j<N; j++) {
                  finalphe[i][j] = finalphe[i][j] + subphe[i][j]*0.25;
                }
            }
            
            for(i=1; i<=subpros; i++){
                MPI_Recv(&recephe[0][0], N*N, MPI_DOUBLE, i, MASTER_TAG,MPI_COMM_WORLD,&status);
                for (j=0; j<N; j++) {
                    for (k=0; k<N; k++) {
                        finalphe[j][k] = finalphe[j][k] + recephe[j][k]*0.25;
                    }
                }
            }
            

            
            // Send updated finalphe to slaves
            for (i=1; i<=subpros; i++) {
                MPI_Send(&finalphe[0][0],N*N, MPI_DOUBLE, i, MASTER_TAG, MPI_COMM_WORLD);
            }
            
            /*
             update subphe in master
             */
            

            
            for (i=0; i<N; i++) {
                for (j=0; j<N; j++) {
                    subphe[i][j] = finalphe[i][j];
                }
            }
        }
        
        if(rank!=0){
            //        printf("rank: %d, aver_ant: %d\n",rank,aver_ant);
            
            memset(subphe,0.0,sizeof(subphe));
            MPI_Recv(&subphe[0][0], N*N, MPI_DOUBLE, 0, MASTER_TAG, MPI_COMM_WORLD, &status);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if(rank!=0){
        MPI_Send(&subbestsolution,1,MPI_DOUBLE,0,MASTER_TAG,MPI_COMM_WORLD);
    }
    
    //
    //
    //    //        printf("%f\n",subbestsolution);
    //
    //
    //
    //
    //
    if(rank==0){
        bestsolution = subbestsolution;
        double recbest;
        for(i=1;i<=subpros;i++){
            MPI_Recv(&recbest,1,MPI_DOUBLE,i,MASTER_TAG,MPI_COMM_WORLD,&status);
            if(recbest<bestsolution){
                bestsolution = recbest;
            }
        }
	double end = MPI_Wtime();	
	float seconds = (float)(end - start);
	printf("time is: %f\n",seconds);
        FILE *fl;
        fl = fopen("out.txt","a");
        Result(fl);
        printf("result saved in out.txt\n");
    }
    //    //     Save result to file
    //    //     */
    //
    
    MPI_Finalize();
    
    return 0;
}
/*
 global methods
 */
void Inputcoords(FILE *fp){
    if(fp==NULL){
        printf("Sorry the file is not exist\n");
        exit(1);
    }
    else{
        for(i=0;i<N;i++){
            fscanf(fp,"%d",&coords[i].city);
            fscanf(fp,"%f",&coords[i].x);
            fscanf(fp,"%f",&coords[i].y);
        }
    }
}
void CreateGraph(){
    double d;
    for(i=0;i<N-1;i++){
        graph[i][i]=MAX;
        for(j=i+1;j<N;j++){
            d = sqrt((double)((coords[i].x-coords[j].x) * (coords[i].x-coords[j].x) + (coords[i].y-coords[j].y) * (coords[i].y-coords[j].y)));
            graph[j][i] = graph[i][j] = d;
        }
    }
    graph[N-1][N-1] = MAX;
    return;
}
double Distance(int *p){
    double d=0;
    for(i=0;i<N-1;i++){
        d+=graph[(*(p+i))][(*(p+i+1))];
    }
    d+=graph[(*(p+N-1))][(*p)];
    return d;
}
void Result(FILE *fl){
    fprintf(fl,"%s\n","Result for this expreiments");
    fprintf(fl,"alpha=%.3lf, beta=%.3lf, rou=%.3lf, Q=%.3lf\n",alpha,beta,rou,Q);
    fprintf(fl,"%s %.4lf\n","Best solution is: ",bestsolution);
    //    fprintf(fl,"%s \n", "Best path is: ");
    /*
     for(i=0;i<N;i++){
     fprintf(fl,"%d →  ", coords[bestway[i]].city);
     }
     fprintf(fl,"%d", coords[bestway[0]].city);
     */
    fprintf(fl, "\n\n\n");
    fclose(fl);
    return;
    
}
/*
 GA methods
 */
void init_gen(int chrom[NUM_ANT][N]){
    int num=0;
    while(num<NUM_ANT){
        for(i=0;i<NUM_ANT;i++)
            for(j=0;j<N;j++)
                chrom[i][j] = j;
        num++;
        for(i=0;i<N-1;i++){
            for(j=i+1;j<N;j++){
                int temp = chrom[num][i];
                chrom[num][i] = chrom[num][j];
                chrom[num][j] = temp;
                num++;
                if(num >= NUM_ANT)
                    break;
            }
            if(num>=NUM_ANT)
                break;
        }
        while(num<NUM_ANT){
            double r1 = ((double)rand())/(RAND_MAX+1.0);
            double r2 = ((double)rand())/(RAND_MAX+1.0);
            int p1 = (int)(N*r1);
            int p2 = (int)(N*r2);
            int temp = chrom[num][p1];
            chrom[num][p1] = chrom[num][p2];
            chrom[num][p2] = temp;
            num++;
        }
    }
}
void Choice(int subchrom[aver_popu][N]){
    double pick;
    double choice_arr[aver_popu][N];
    double fit_pro[aver_popu];
    double sum = 0;
    double fit[aver_popu];
    for(j=0;j<aver_popu;j++){
        double path = Distance(subchrom[j]);
        double fitness = 1/path;
        fit[j] = fitness;
        sum += fitness;
    }
    for(j=0;j<aver_popu;j++){
        fit_pro[j] = fit[j]/sum;
    }
    for(i=0;i<aver_popu;i++){
        pick = ((double)rand())/RAND_MAX;
        for(j=0;j<aver_popu;j++){
            pick = pick - fit_pro[j];
            if(pick<=0){
                for(k=0;k<N;k++)
                    choice_arr[i][k] = subchrom[j][k];
                break;
            }
        }
    }
    for(i=0;i<aver_popu;i++){
        for(j=0;j<N;j++)
            subchrom[i][j] = choice_arr[i][j];
    }
}
void Cross(int subchrom[aver_popu][N]){
    double pick;
    double pick1,pick2;
    int choice1,choice2;
    int pos1,pos2;
    int temp;
    int conflict1[N];
    int conflict2[N];
    int num1,num2;
    int index1,index2;
    int move = 0;
    while(move<aver_popu-1){
        pick = ((double)rand())/RAND_MAX;
        if(pick > P_CROSS){
            move += 2;
            continue;
        }
        
        choice1 = move;
        choice2 = move+1;
        pick1 = ((double)rand())/(RAND_MAX+1.0);
        pick2 = ((double)rand())/(RAND_MAX+1.0);
        pos1 = (int)(pick1*N);
        pos2 = (int)(pick2*N);
        while(pos1 > N -2 || pos1 < 1){
            pick1 = ((double)rand())/(RAND_MAX+1.0);
            pos1 = (int)(pick1*N);
        }
        while(pos2 > N -2 || pos2 < 1){
            pick2 = ((double)rand())/(RAND_MAX+1.0);
            pos2 = (int)(pick2*N);
        }
        if(pos1 > pos2){
            temp = pos1;
            pos1 = pos2;
            pos2 = temp;
        }
        for(j=pos1;j<=pos2;j++){
            temp = subchrom[choice1][j];
            subchrom[choice1][j] = subchrom[choice2][j];
            subchrom[choice2][j] = temp;
        }
        num1 = 0;
        num2 = 0;
        if(pos1 > 0 && pos2 < N-1){
            for(j =0;j<=pos1-1;j++){
                for(k=pos1;k<=pos2;k++){
                    if(subchrom[choice1][j] == subchrom[choice1][k]){
                        conflict1[num1] = j;
                        num1++;
                    }
                    if(subchrom[choice2][j] == subchrom[choice2][k]){
                        conflict2[num2] = j;
                        num2++;
                    }
                }
            }
            for(j=pos2+1;j<N;j++){
                for( k=pos1;k<=pos2;k++){
                    if(subchrom[choice1][j] == subchrom[choice1][k]){
                        conflict1[num1] = j;
                        num1++;
                    }
                    if(subchrom[choice2][j] == subchrom[choice2][k]){
                        conflict2[num2] = j;
                        num2++;
                    }
                }
            }
        }
        if((num1 == num2) && num1 > 0){
            for(j=0;j<num1;j++){
                index1 = conflict1[j];
                index2 = conflict2[j];
                temp = subchrom[choice1][index1];
                subchrom[choice1][index1] = subchrom[choice2][index2];
                subchrom[choice2][index2] = temp;
            }
        }
        move += 2;
    }
}

void Mutation(int subchrom[aver_popu][N]){
    double pick,pick1,pick2;
    int pos1,pos2,temp;
    for(i=0;i<aver_popu;i++){
        pick = ((double)rand())/RAND_MAX;
        if(pick > P_MUT)
            continue;
        pick1 = ((double)rand())/(RAND_MAX+1.0);
        pick2 = ((double)rand())/(RAND_MAX+1.0);
        pos1 = (int)(pick1*N);
        pos2 = (int)(pick2*N);
        while(pos1 > N-1){
            pick1 = ((double)rand())/(RAND_MAX+1.0);
            pos1 = (int)(pick1*N);
        }
        while(pos2 > N-1){
            pick2 = ((double)rand())/(RAND_MAX+1.0);
            pos2 = (int)(pick2*N);
        }
        temp = subchrom[i][pos1];
        subchrom[i][pos1] = subchrom[i][pos2];
        subchrom[i][pos2] = temp;
    }
}
void Reverse(int subchrom[aver_popu][N]){
    double pick1,pick2;
    double dis,reverse_dis;
    int n;
    int flag,pos1,pos2,temp;
    int reverse_arr[N];
    int i,j;
    for(i=0;i<aver_popu;i++){
        flag = 0;
        while(flag == 0){
            pick1 = ((double)rand())/(RAND_MAX+1.0);
            pick2 = ((double)rand())/(RAND_MAX+1.0);
            pos1 = (int)(pick1*N);
            pos2 = (int)(pick2*N);
            while(pos1 > N-1){
                pick1 = ((double)rand())/(RAND_MAX+1.0);
                pos1 = (int)(pick1*N);
            }
            while(pos2 > N -1){
                pick2 = ((double)rand())/(RAND_MAX+1.0);
                pos2 = (int)(pick2*N);
            }
            if(pos1 > pos2){
                temp = pos1;
                pos1 = pos2;
                pos2 = temp;
            }
            if(pos1 < pos2){
                for(j=0;j<N;j++)
                    reverse_arr[j] = subchrom[i][j];
                n = 0;
                for(j=pos1;j<=pos2;j++){
                    reverse_arr[j] = subchrom[i][pos2-n];
                    n++;
                }
                reverse_dis = Distance(reverse_arr);
                dis = Distance(subchrom[i]);
                if(reverse_dis < dis){
                    for(j=0;j<N;j++)
                        subchrom[i][j] = reverse_arr[j];
                }
            }
            flag = 1;
        }
    }
}
/*
 AC methods
 */
void init_ac(){
    alpha = 1;
    beta = 1;
    rou = 0.7;
    Q = 10;
    NcMax = RC_ITER;
    return;
}
