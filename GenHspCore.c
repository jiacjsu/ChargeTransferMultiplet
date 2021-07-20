#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int sumEveryBit(unsigned long int tempa){
    long int tempsum=0;
    unsigned long int tempa1=tempa;
    while(tempa1){
        tempsum += tempa1 % 2;
        tempa1 = tempa1 >> 1;
    }
    return tempsum;
}


/* When the number of ligands orbital approches 9, the calculation becomes really slow,
   since 2 ** (2*N) traversing is very slow...

int GenHsp2p1HoleCore(int N_ligand, int nup, int ndn, unsigned int *Hsp){

    int N_core = 3;
    int N_d = 5;
    int N = N_ligand + N_core + N_d;
    unsigned int flagCoreFull = 7 * (1 << (N-N_core)) * (1 + (1 << N));
    int ntemp = 0;
    for(int stat=0; stat<pow(2,2*N); stat++){
        if(sumEveryBit(flagCoreFull & stat) == 5 && \
           sumEveryBit(stat) == (nup+ndn)){
            Hsp[ntemp] = stat;
            ntemp++;     
        }
    }
    return ntemp;
}

 */

int cmpfunc (const void * a, const void * b) {
    //return ( *(long int*)a - *(long int*)b );
    if(*(long int*)a > *(long int*)b) return 1;
    else if(*(long int*)a < *(long int*)b) return -1;
    else return 0;
}
// find this problem, needs to be changed, because compare long int, return int


int GenHsp2p1HoleCore(int N_ligand, int nup, int ndn, unsigned long int *Hsp){

    int N_core = 3;
    int N_d = 5;
    int N = N_ligand + N_core + N_d;
    unsigned long int temp = 7 * (1 << (N-N_core));
    //unsigned long int flagCoreFull = 7 * (1 << (N-N_core)) * (1 + (1 << N));
    unsigned long int flagCoreFull = temp;
    flagCoreFull = (flagCoreFull << N) + flagCoreFull; 
    unsigned long int tlow = (1 << (N-N_core));
    unsigned long int thigh = tlow;
    thigh = thigh << N;
    unsigned long int corestates[6] = {3 * thigh + 7 * tlow,
                                       5 * thigh + 7 * tlow,
                                       6 * thigh + 7 * tlow,
                                       7 * thigh + 3 * tlow,
                                       7 * thigh + 5 * tlow,
                                       7 * thigh + 6 * tlow};

    //printf("%li\n", flagCoreFull);
    //for(int i=0; i<6; i++) printf("%li\n", corestates[i]);
    //printf("\n");  

    int ntemp = 0;
    for(long int stat=0; stat<pow(2,2*(N-N_core)); stat++){
        if(sumEveryBit(stat) == (nup+ndn-5)){
            unsigned long int statshift = ((stat >> (N-N_core)) << N) | (stat % (1 << (N-N_core)));
            for(int i=0; i<6; i++){
                unsigned long int statfinal = statshift | corestates[i];
                Hsp[ntemp] = statfinal;
                ntemp++;                       
                //printf("%d, %li\n", ntemp, statfinal);
            }
        }
    }
    qsort(Hsp, ntemp, sizeof(unsigned long int), cmpfunc);
    return ntemp;
}
