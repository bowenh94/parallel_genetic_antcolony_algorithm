#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

/*
 Global parameters
 */
#define N 194
#define MAX 0x7fffffff

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
    
    clock_t start = clock();
    /*Do something*/
    
    /*
     initial variables
     */
    memset(solution,0.0,sizeof(solution));
    memset(bestway, 0, sizeof(bestway));
    
    /*
     GA part
     */
    int chrom[NUM_ANT][N];
    memset(chrom, 0, sizeof(chrom));
    init_gen(chrom);
    
    for(j=0;j<NUM_ANT;j++){
        solution[j] = Distance(chrom[j]);
        if(solution[j] < bestsolution){
            bestsolution = solution[j];
            for(i=0;i<N;i++){
                bestway[i] = chrom[j][i];
            }
        }
    }
    
    int num_gen = 0;
    while(num_gen<MAX_GEN){
        //initial solution each generation
        memset(solution, 0, sizeof(solution));
        
        //do GA operation on chrom
        Choice(chrom);
        Cross(chrom);
        Mutation(chrom);
        Reverse(chrom);
        
        for(j=0;j<NUM_ANT;j++){
            solution[j] = Distance(chrom[j]);
            if(solution[j] < bestsolution){
                bestsolution = solution[j];
                for(k=0;k<N;k++){
                    bestway[k] = chrom[j][k];
                }
            }
        }
        num_gen+=1;
    }

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
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            phe[i][j] = phe[i][j]*rou + add[i][j];
            if(phe[i][j] < tao_min)
                phe[i][j] = tao_min;
            else if(phe[i][j] > tao_max)
                phe[i][j] = tao_max;
        }
    }
    
    
    /*
     AC algorithm
     */
    
    int NC = 0;
    int s;
    double drand, pro, psum;
    
    memset(map, -1, sizeof(map));
    memset(vis, 0, sizeof(vis));
    while(NC++ <= NcMax){
        for(k=0;k<NUM_ANT;k++){
            map[k][0] = (k+NC) % N;
            vis[k][map[k][0]] = 1;
        }
        s=1;
        while(s<N){
            for(k=0; k<NUM_ANT; ++k){
                psum = 0;
                for(j=0; j<N; ++j){
                    if(vis[k][j] == 0){
                        psum += pow(phe[ map[k][s-1] ][j], alpha) * pow(yita[ map[k][s-1] ][j], beta);
                    }
                }
                drand = (double)(rand() % 5000);
                drand /= 5000.0;
                pro = 0;
                for(j=0; j<N; ++j){
                    if(vis[k][j] == 0){
                        pro += pow(phe[map[k][s-1]][j], alpha) * pow(yita[map[k][s-1]][j], beta) / psum;
                    }
                    if(pro > drand)
                        break;
                }
                vis[k][j] = 1;
                map[k][s] = j;
            }
            s++;
        }
        
//        for (i=0; i<NUM_ANT; i++) {
//            for (j=0; j<N; j++) {
//                printf("%d ",map[i][j]);
//            }
//            printf("\n");
//        }
        memset(add,0,sizeof(add));
        /*
         for(i=0;i<N;i++){
         printf("%d -> ",bestway[i]);
         }
         */
        memset(solution,0,sizeof(solution));
        for(j=0;j<NUM_ANT;j++){
            solution[j] = Distance(map[j]);
            if(solution[j] < bestsolution){
                bestsolution = solution[j];
                for(k=0;k<N;k++){
                    bestway[k] = map[j][k];
                }
            }
        }
        for(k=0; k<NUM_ANT; k++){
            for(j=0; j<N-1; j++){
                add[ map[k][j] ][ map[k][j+1] ] += Q/solution[k];
            }
            add[ map[k][N-1] ][ map[k][0] ] += Q/solution[k];
        }
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                phe[i][j] = phe[i][j]*rou + add[i][j];
                if(phe[i][j] < tao_min)
                    phe[i][j] = tao_min;
                else if(phe[i][j] > tao_max)
                    phe[i][j] = tao_max;
            }
        }
        memset(vis,0,sizeof(vis));
        memset(map,-1,sizeof(map));
    }
    /*
     for(i=0;i<N;i++){
     printf("%d -> ",bestway[i]);
     }
     */
    
    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
    printf("time is: %f\n",seconds);
    
    /*
     Save result to file
     */
    FILE *fl;
    fl = fopen("out.txt","a");
    Result(fl);
    printf("result saved in out.txt\n");
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
    fprintf(fl,"%s \n", "Best path is: ");
    
    for(i=0;i<N;i++){
        fprintf(fl,"%d →  ", coords[bestway[i]].city);
    }
    fprintf(fl,"%d", coords[bestway[0]].city);
    
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
void Choice(int chrom[NUM_ANT][N]){
    double pick;
    double choice_arr[NUM_ANT][N];
    double fit_pro[NUM_ANT];
    double sum = 0;
    double fit[NUM_ANT];
    for(j=0;j<NUM_ANT;j++){
        double path = Distance(chrom[j]);
        double fitness = 1/path;
        fit[j] = fitness;
        sum += fitness;
    }
    for(j=0;j<NUM_ANT;j++){
        fit_pro[j] = fit[j]/sum;
    }
    for(i=0;i<NUM_ANT;i++){
        pick = ((double)rand())/RAND_MAX;
        for(j=0;j<NUM_ANT;j++){
            pick = pick - fit_pro[j];
            if(pick<=0){
                for(k=0;k<N;k++)
                    choice_arr[i][k] = chrom[j][k];
                break;
            }
        }
    }
    for(i=0;i<NUM_ANT;i++){
        for(j=0;j<N;j++)
            chrom[i][j] = choice_arr[i][j];
    }
}
void Cross(int chrom[NUM_ANT][N]){
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
    while(move<N-1){
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
            temp = chrom[choice1][j];
            chrom[choice1][j] = chrom[choice2][j];
            chrom[choice2][j] = temp;
        }
        num1 = 0;
        num2 = 0;
        if(pos1 > 0 && pos2 < N-1){
            for(j =0;j<=pos1-1;j++){
                for(k=pos1;k<=pos2;k++){
                    if(chrom[choice1][j] == chrom[choice1][k]){
                        conflict1[num1] = j;
                        num1++;
                    }
                    if(chrom[choice2][j] == chrom[choice2][k]){
                        conflict2[num2] = j;
                        num2++;
                    }
                }
            }
            for( j=pos2+1;j<N;j++){
                for(k=pos1;k<=pos2;k++){
                    if(chrom[choice1][j] == chrom[choice1][k]){
                        conflict1[num1] = j;
                        num1++;
                    }
                    if(chrom[choice2][j] == chrom[choice2][k]){
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
                temp = chrom[choice1][index1];
                chrom[choice1][index1] = chrom[choice2][index2];
                chrom[choice2][index2] = temp;
            }
        }
        move += 2;
    }
}
void Mutation(int chrom[NUM_ANT][N]){
    double pick,pick1,pick2;
    int pos1,pos2,temp;
    for(i=0;i<NUM_ANT;i++){
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
        temp = chrom[i][pos1];
        chrom[i][pos1] = chrom[i][pos2];
        chrom[i][pos2] = temp;
    }
}
void Reverse(int chrom[NUM_ANT][N]){
    double pick1,pick2;
    double dis,reverse_dis;
    int n;
    int flag,pos1,pos2,temp;
    int reverse_arr[N];
    int i,j;
    for(i=0;i<NUM_ANT;i++){
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
                    reverse_arr[j] = chrom[i][j];
                n = 0;
                for(j=pos1;j<=pos2;j++){
                    reverse_arr[j] = chrom[i][pos2-n];
                    n++;
                }
                reverse_dis = Distance(reverse_arr);
                dis = Distance(chrom[i]);
                if(reverse_dis < dis){
                    for(j=0;j<N;j++)
                        chrom[i][j] = reverse_arr[j];
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
