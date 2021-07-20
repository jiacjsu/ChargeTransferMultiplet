#include <stdio.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// C does not have reference, need to pass a pointer


// np.zeros((HsizeEst), dtype = np.uint32)
//    hop = np.zeros((2*ntol, 2*ntol), dtype=np.complex64) 
//    U_ext = np.zeros((ntol, ntol))
//    J_ext = np.zeros((ntol, ntol))
//    U_onsite = np.zeros((ntol))
//        U_pddp = np.zeros((3,5,5,3))
//    U_dpdp = np.zeros((5,3,5,3))


// to implement: sumEverybit()


int ntmax = 28;  // be consistent with ntolmax in SetModelPara.py


int binary_conversion(int num)
{
    if (num == 0)
    {
        return 0;
    }
    else
    {
        return (num % 2) + 10 * binary_conversion(num / 2);
    }
}



int sumEveryBit(unsigned long int tempa){
    int tempsum=0;
    unsigned long int tempa1=tempa;
    while(tempa1){
        tempsum += tempa1 % 2;
        tempa1 = tempa1 >> 1;
    }
    return tempsum;
}     

int BinarySearch(unsigned long int array[], int num, unsigned long int keynum){
    int low = 0;
    int high = num-1;
    int mid;
    do
    {
        mid = (low + high) / 2;
        if ((long int)keynum < (long int)(array[mid]))
            high = mid - 1;
        else if ((long int)keynum > (long int)(array[mid]))
            low = mid + 1;
    } while ((long int)keynum != (long int)(array[mid]) && low <= high);
    if ((long int)keynum == (long int)(array[mid]))
    {
        //printf("SEARCH SUCCESSFUL \n");
        return mid;
    }
    else
    {
        //printf("SEARCH FAILED \n");
        return -1;
    }
}


/*                
# t_{pp,kk} c^\dagger_{kk,\sigma} c_{pp,\sigma}
# U_SO{pp,kk} c^\dagger_{kk,\sigma} c_{pp,-\sigma}
# t_{pp,kk} c^\dagger_{kk,\sigma} c_{pp, \sigma^{\prime}}


#hopping (same spin) and Fe 2p spin-orbital coupling (different spin)
 */

void HopHami(unsigned long int statupdn, double complex *hop, int N, bool corr, 
             double *E_corr, double complex *Values, unsigned long int *States, int *ntag){
    
    /*
    printf("We are in HopHami\n");
    printf(" corr is %d\n", corr);
    printf(" N is %d\n", N);
    printf(" ntmax is %d\n", ntmax);
    printf(" statupdn %u\n", statupdn);
    */

    float hoptol = 0.01;
    int pp,kk,tempsign;
    unsigned long int istatupdn, fstatupdn;
    float tvalue = 0.0;
    for(int pp=0; pp < 2*N; pp++){ 
        unsigned long int opp = 1;
        opp = opp << pp;
        //printf("inside loop pp %d\n", pp);

        for(int kk=0; kk < 2*N; kk++){ 
            unsigned long int okk = 1;
            okk = okk << kk;
            //if (cabs(hop[kk][pp]) < hoptol) continue;
            if (cabs(hop[kk*(2*ntmax)+pp]) < hoptol) continue;
            if (statupdn & opp){
                istatupdn = statupdn & (~opp);

                //printf(" istatupdn %u\n", istatupdn);

                tempsign = pow(-1, sumEveryBit(statupdn >> (pp+1)));

                //printf(" tempsign %d\n", tempsign);

                if(!(istatupdn & okk)){
                    fstatupdn = istatupdn | okk;
                
                    //printf(" fstatupdn %u\n", fstatupdn);

                    tempsign *= pow(-1,sumEveryBit(istatupdn >> (kk+1)));

                    //printf(" tempsign %d\n", tempsign);

                    States[*ntag] = fstatupdn;
                    Values[*ntag] = tempsign * hop[kk*(2*ntmax)+pp];
                    (*ntag)++;
                    //printf(" HopHami *ntag %li %li %li %d\n", statupdn, istatupdn, fstatupdn, *ntag);
                }        
            } 
        }
    }
                       
    if(corr) {
        double tvalue = 0.0;
        for(int pp=0; pp < 2*N; pp++){
            unsigned long int opp = 1;
            opp = opp << pp;
            if(statupdn & opp){
                tvalue += E_corr[pp];
            }
        }
        States[*ntag] = statupdn;
        Values[*ntag] = tvalue;
        (*ntag)++;
    }
}


    
void UonsiteHami(unsigned long int statupdn, double *U_onsite, int N,
                 double complex *Values, unsigned long int *States, int *ntag){

    //printf("We are in UonsiteHami\n");

    double tvalue = 0.0;
    int N_core = 3;
    int N_d = 5;
    for(int pp = 0; pp<N; pp++){
        int kk = pp + N;
        unsigned long int okk = 1;
        okk = okk << kk;
        if((statupdn & (1 << pp)) && (statupdn & okk)){
            tvalue += U_onsite[pp];
        }
    }
    States[*ntag] = statupdn;
    Values[*ntag] = tvalue;
    (*ntag)++;
}
    

void UextHami(unsigned long int statupdn, double *U_ext, int N,\
              double complex *Values, unsigned long int *States, int *ntag){

    //printf("We are in UextHami\n");

    double totalvalue = 0.0;
    for(int pp=0; pp<2*N; pp++){
        unsigned long int opp = 1;
        opp = opp << pp;
        for(int kk=0; kk<2*N; kk++){
            unsigned long int okk = 1;
            okk = okk << kk;
            if (pp % N != kk % N){
                if ((statupdn & opp) && (statupdn & okk)){
                    totalvalue += (0.5 * U_ext[(pp % N)*(ntmax)+(kk % N)]);
                }
            }
        }       
    }
    States[*ntag] = statupdn;
    Values[*ntag] = totalvalue;
    (*ntag)++;
}


/*
# from Fortran script
# Adam note "multiplets.pdf" Equation 9

#1. J_{pp,kk} c^\dagger_{kk,dn} c^\dagger_{pp, dn} c_{kk,dn} c_{pp, dn}
#2. J_{pp,kk} c^\dagger_{kk,up} c^\dagger_{pp, up} c_{kk,up} c_{pp, up}
#3. J_{pp,kk} c^\dagger_{kk,dn} c^\dagger_{pp, up} c_{kk,up} c_{pp, dn}
#4. J_{pp,kk} c^\dagger_{kk,up} c^\dagger_{pp, dn} c_{kk,dn} c_{pp,up}

# J_{pp,kk} c^\dagger_{kk,dn} c^\dagger_{kk,up} c_{pp,up} c_{pp,dn}
# J_{pp,kk} c^\dagger_{kk,up} c^\dagger_{kk,dn} c_{pp,dn} c_{pp,up}

# part 1: from right to left:
# J_{pp,kk} c^\dagger_{kkp,ms} c^\dagger_{ppp, msp} c_{kk,msp} c_{pp, ms}
# part 2: from right to left:
# J_{pp,kk} c^\dagger_{kkp,ms} c^\dagger_{kk, msp} c_{ppp,msp} c_{pp, ms}
*/

void JextHami(unsigned long int statupdn, double *J_ext, int N,
              double complex *Values, unsigned long int *States, int *ntag){

    //printf("We are in JextHami\n");
    /*
    for(int i=0; i<ntmax; i++)
        for(int j=0; j<ntmax; j++)
            printf("    J_ext at %d,%d is %f\n", i, j, J_ext[i*ntmax + j]);
     */


    double Jexttol = 0.01;
    int tempsign;
    unsigned long int istatupdn, fstatupdn;
    for(int pp=0; pp<2*N; pp++){
        unsigned long int opp = 1;
        opp = opp << pp;
        for(int kk=0; kk<2*N; kk++){
            if ((pp % N == kk % N) || (fabs(J_ext[(kk%N) * (ntmax) + (pp%N)]) < Jexttol)) continue;

            unsigned long int okk = 1;
            okk = okk << kk;
            //printf(" pp %d kk %d J_ext %f\n", pp, kk, J_ext[(kk%N) * (ntmax) + (pp%N)]);

            int ms = pp / N;
            int msp = kk / N;
            int ppp = msp * N + pp % N;
            int kkp = ms * N + kk % N;
            unsigned long int okkp = 1, oppp = 1;
            okkp = okkp << kkp;
            oppp = oppp << ppp;

            if ((statupdn & opp) && (statupdn & okk)){
                istatupdn = statupdn & (~opp);
                tempsign = pow(-1, sumEveryBit(statupdn >> (pp+1)));
                istatupdn = istatupdn & (~okk);
                tempsign *= pow(-1, sumEveryBit(istatupdn >> (kk+1)));
                
                if ((!(istatupdn & oppp)) && (!(istatupdn & okkp))){
                    fstatupdn = istatupdn | oppp;
                    tempsign *= pow(-1, sumEveryBit(istatupdn >> (ppp+1)));
                    fstatupdn = fstatupdn | okkp;
                    tempsign *= pow(-1, sumEveryBit(fstatupdn >> (kkp+1))); 
 
                    States[*ntag] = fstatupdn;
                    Values[*ntag] = (0.5 * tempsign * J_ext[(pp % N)*ntmax + (kk % N)]);
                    (*ntag)++;
                }
            }
        }
    }                               

    for(int pp=0; pp<2*N; pp++){
        unsigned long int opp = 1;
        opp = opp << pp;
        for(int kk=0; kk<2*N; kk++){
            if((pp % N == kk % N) || (pp / N == kk / N) || 
                (fabs(J_ext[(kk%N) * ntmax + (pp%N)]) < Jexttol)) continue;
            unsigned long int okk = 1;
            okk = okk << kk;
            int ms = pp / N;
            int msp = kk / N;
            int ppp = msp * N + pp % N;
            int kkp = ms * N + kk % N;
            unsigned long int okkp = 1, oppp = 1;
            okkp = okkp << kkp;
            oppp = oppp << ppp;

            if((statupdn & opp) && (statupdn & oppp)){
                istatupdn = statupdn & (~opp);
                tempsign = pow(-1, sumEveryBit(statupdn >> (pp+1)));
                istatupdn = istatupdn & (~oppp);
                tempsign *= pow(-1, sumEveryBit(istatupdn >> (ppp+1)));
                
                if ((!(istatupdn & okk)) && (!(istatupdn & okkp))){
                    fstatupdn = istatupdn | okk;
                    tempsign *= pow(-1, sumEveryBit(istatupdn >> (kk+1)));
                    fstatupdn = fstatupdn | okkp;
                    tempsign *= pow(-1, sumEveryBit(fstatupdn >> (kkp+1)));
                                
                    States[*ntag] = fstatupdn;
                    Values[*ntag] = (0.5 * tempsign * J_ext[(pp % N) * ntmax + (kk % N)]);
                    (*ntag)++;
                }
            }
        }
    }

}
 

//void ddMulti(statupdn, U_rest, N){
void ddMulti(unsigned long int statupdn, double *U_rest, int N,
             double complex *Values, unsigned long int *States, int *ntag){
                                                          
    float tol = 0.01;
    int N_d = 5;
    
    int appp, app, ap, a, tempsign;
    unsigned long int istatupdn, fstatupdn;
    for(int apppl=0; apppl < N_d; apppl++){ //appp in range(N_d) + range(N, N+N_d):
        for(int appl=0; appl < N_d; appl++){ // app in range(N_d) + range(N, N+N_d):
            //ms = appp / N
            //msp = app / N
            for(int ms = 0; ms < 2; ms++){
                for(int msp = 0; msp < 2; msp++){
                    int appp = apppl + N*ms;
                    int app = appl + N*msp;
                    if(appp == app) continue;
            
                    unsigned long int oappp = 1, oapp = 1;
                    oappp = oappp << appp;
                    oapp = oapp << app;
 
                    for(int ap = msp*N; ap < msp*N + N_d; ap++){
                        for(int a = ms*N; a < ms*N + N_d; a++){
                            if(ap == a) continue;

                            unsigned long int oap = 1, oa = 1;
                            oap = oap << ap;
                            oa = oa << a;

                            if(fabs(U_rest[(a%N) * 125 + (ap%N) * 25 + (app%N) * 5 + (appp%N)]) < tol) continue;
                            
                            if((statupdn & oappp) && (statupdn & oapp)){
                                istatupdn = statupdn & (~oappp);
                                tempsign = pow(-1, sumEveryBit(statupdn >> (appp+1)));
                                istatupdn = istatupdn & (~(oapp));
                                tempsign *= pow(-1, sumEveryBit(istatupdn >> (app+1)));
                                
                                if((!(istatupdn & (oap))) && (!(istatupdn & (oa)))){
                                    fstatupdn = istatupdn | (oap);
                                    tempsign *= pow(-1, sumEveryBit(istatupdn >> (ap+1)));
                                    fstatupdn = fstatupdn | (oa);
                                    tempsign *= pow(-1, sumEveryBit(fstatupdn >> (a+1)));  
                                     
                                    States[*ntag] = fstatupdn;
                                    Values[*ntag] = (0.5 * tempsign * U_rest[(a%N) * 125 + (ap%N) * 25 + (app%N) * 5 + (appp%N)]);
                                    (*ntag)++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
 
/*  
# Adam's note multiplets_2p3d.pdf equation 2
# Upddp_{b,ap,app,bppp} p^\dagger_{b,ms} d^\dagger_{ap, msp} d_{app,msp} p_{bppp, ms}
# Udpdp_{b,bp,app,bppp} d^\dagger_{b,ms} p^\dagger_{bp, msp} d_{app,msp} p_{bppp, ms}  

# this version p orbital is expressed as the hole language
# use Eequation 3 from Adam's note for hole language for p orbitals
# -Upddp_{b,ap,app,bppp} p_{bppp, ms} d^\dagger_{ap, msp} d_{app,msp} p^\dagger_{b,ms}
# -Udpdp_{b,bp,app,bppp} d^\dagger_{b,ms} p_{bppp, ms} d_{app,msp} p^\dagger_{bp, msp}
*/
void dpMulti(unsigned long int statupdn, double *U_pddp, double *U_dpdp, int N,
             double complex *Values, unsigned long int *States, int *ntag){

    //printf("We are in dpMulti\n");

    float tol = 0.01;
    int N_core = 3;
    int N_d = 5;
    unsigned long int flagCoreFull = 7 * (1 << (N-N_core));
    flagCoreFull = flagCoreFull + (flagCoreFull << N);

    if(sumEveryBit(statupdn & flagCoreFull) == (N_core * 2)) return;
        
    int bppp, app, ap, b, bp, a;
    unsigned long int istatupdn, fstatupdn;
    int tempsign;
    for(int b=N-N_core; b<2*N; b++){
        if ( b>=N && b<(2*N-N_core) ) continue;   

        unsigned long int ob = 1;
        ob = ob << b; 

        if (statupdn & ob) continue;

        for(int app=0; app<N+N_d; app++){ 
            if( app >= N_d && app < N) continue;

            unsigned long int oapp = 1;
            oapp = oapp << app;   
            if (! (statupdn & oapp)) continue;

            int ms = b / N;
            int msp = app / N;
            for(int ap=msp*N; ap < msp*N+N_d; ap++){ 
                unsigned long int oap = 1;
                oap = oap << ap;

                for(int bppp=ms*N + N - N_core; bppp < ms*N + N; bppp++){  
                    if (fabs(U_pddp[(b%N - (N-N_core)) * 5*5*3 + \
                        (ap%N) * 5*3 + (app%N) * 3 + (bppp%N-(N-N_core))]) < tol) continue;
                    
                    unsigned long int obppp = 1;
                    obppp = obppp << bppp;

                    istatupdn = statupdn | ob;
                    tempsign = pow(-1, sumEveryBit(statupdn >> (b+1)));
                    istatupdn = istatupdn & (~oapp);
                    tempsign *= pow(-1, sumEveryBit(istatupdn >> (app+1)));
                    
                    if ((!(istatupdn & oap)) && (istatupdn & obppp)){
                        fstatupdn = istatupdn | oap;
                        tempsign *= pow(-1, sumEveryBit(istatupdn >> (ap+1)));
                        fstatupdn = fstatupdn & (~obppp);
                        tempsign *= pow(-1, sumEveryBit(fstatupdn >> (bppp+1)));  
             
                        States[*ntag] = fstatupdn;
                        Values[*ntag] = (-1 * tempsign * U_pddp[(b%N - (N-N_core)) * 5*5*3 + \
                                                (ap%N) * 5*3 + (app%N)*3 + (bppp%N - (N-N_core))]);
                        //printf(" A %2d C %2d A %2d C %2d | %f7.3 pddp\n", bppp, ap, app, b, creal(Values[*ntag]));
                        (*ntag)++;
                    }
                }
            }
        }
    }

    for(int bp=N-N_core; bp<2*N; bp++){
        if( bp>=N && bp < 2*N-N_core) continue;
        unsigned long int obp = 1;
        obp = obp << bp;
        if (statupdn & obp) continue;
        int msp = bp / N;

        for(int app=msp*N; app < msp*N + N_d; app++){

            unsigned long int oapp = 1;
            oapp = oapp << app;
            if (!(statupdn & oapp)) continue;
        
            for(int bppp=N-N_core; bppp < 2*N; bppp++ ){
                if(bppp >= N && bppp < 2*N-N_core) continue; 
                int ms = bppp / N;

                unsigned long int obppp = 1;
                obppp = obppp << bppp;

                for(int a=ms*N; a < ms*N+N_d; a++){
                    
                    if (fabs(U_dpdp[(a%N) * 3*5*3 + (bp%N-(N-N_core))*5*3 + 
                                    (app%N) * 3 + (bppp%N-(N-N_core))]) < tol) continue;
                    
                    unsigned long int oa = 1;
                    oa = oa << a;

                    istatupdn = statupdn | (obp);
                    tempsign = pow(-1, sumEveryBit(statupdn >> (bp+1)));
                    istatupdn = istatupdn & (~(oapp));
                    tempsign *= pow(-1, sumEveryBit(istatupdn >> (app+1)));
                    
                    if ((istatupdn & (obppp)) && (!(istatupdn & (oa)))){
                        fstatupdn = istatupdn & (~(obppp));
                        tempsign *= pow(-1, sumEveryBit(istatupdn >> (bppp+1)));
                        fstatupdn = fstatupdn | (oa);
                        tempsign *= pow(-1, sumEveryBit(fstatupdn >> (a+1)));  
        
                        States[*ntag] = fstatupdn;
                        Values[*ntag] = (-1 * tempsign * U_dpdp[(a%N) * 3*5*3 + 
                                        (bp%N - (N-N_core)) * 5*3 + (app%N) * 3 + 
                                        (bppp%N - (N-N_core))]);
                        //printf(" C %2d C %2d A %2d A %2d | %f7.3\n", a, bp, bppp, app, creal(Values[*ntag]));
                        (*ntag)++;
                    }
                }
            }
        }
    }
    // C*Nd, where C = -6F0 + 6G1 + 63G3
}
                                                                            
    

const char *byte_to_binary(int x)
{
    static char b[9];
    b[0] = '\0';

    int z;
    for (z = 128; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}
      

int GenMatrixCore(int Hsize, int N, double complex *hop, double *U_onsite, 
                    double *U_ext, double *J_ext, double *U_rest, double *U_pddp, double *U_dpdp, 
                    unsigned long int *Hsp, bool corr, double *E_corr, 
                    double complex *Hvalue, int *IndexI, int *IndexJ){

    int nsparse=0;

    /*
    printf(" We are in the C core codes now \n");
    printf(" Hsize value is %d\n", Hsize);


    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            printf("inside GEnMatrixCore J_ext %d,%d, %f\n", i, j, J_ext[i*ntmax + j]);

    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            printf("inside GEnMatrixCore U_ext %d,%d, %f\n", i, j, U_ext[i*ntmax + j]);

    for(int i=0; i<5; i++)
        for(int j=0; j<3; j++)
            for(int k=0; k<5; k++)
                for(int l=0; l<3; l++)
                    printf("inside GEnMatrixCore U_dpdp %d%d%d%d, %f\n", i, j, k, l, \
                           U_dpdp[i*3*5*3 + j*5*3 + k*3 + l]);

    for(int kk=0; kk<N; kk++)
        for(int pp=0; pp<N; pp++)
            printf("inside GenMatrixCore hop[%d,%d] is %f+i%f \n", kk, pp, creal(hop[kk*(2*ntmax)+pp]), \
                   cimag(hop[kk*(2*ntmax)+pp]));
    */

    for(int i=0; i<Hsize; i++){
    //for(int i=0; i<2; i++){

        unsigned long int statupdn = Hsp[i];
        int ntag = 0, nttol = 1000; //ntag: Values, States have been filled up to this index 
        double complex Values[nttol];
        unsigned long int States[nttol];
        double complex ValuesCompact[nttol];
        unsigned long int StatesCompact[nttol];

        HopHami(statupdn, hop, N, corr, E_corr, Values, States, &ntag);
        //printf("after HopHami ntag is %d\n", ntag);
        UonsiteHami(statupdn, U_onsite, N, Values, States, &ntag);
        //printf("after UonsiteHami ntag is %d\n", ntag);
        UextHami(statupdn, U_ext, N, Values, States, &ntag);
        //printf("after UextHami ntag is %d\n", ntag);
        JextHami(statupdn, J_ext, N, Values, States, &ntag);
        //printf("after JextHami ntag is %d\n", ntag);
        ddMulti(statupdn, U_rest, N, Values, States, &ntag);
        //printf("after dpMulti ntag is %d\n", ntag);
        dpMulti(statupdn, U_pddp, U_dpdp, N, Values, States, &ntag);
        //printf("after dpMulti ntag is %d\n", ntag);


        int ntagCompact = 0;
        for(int j=0; j<ntag; j++){

            int index = 0;
            bool stateinCompact = false;
            for(int k=0; k<ntagCompact; k++){
                if(States[j] == StatesCompact[k]) {
                    stateinCompact = true;
                    index = k;
                }
            }

            if(stateinCompact){
                ValuesCompact[index] += Values[j];
            } else {
                StatesCompact[ntagCompact] = States[j];
                ValuesCompact[ntagCompact] = Values[j];
                ntagCompact++;
            }
        }

        // Binary Search and collect to the IndexI, IndexJ, Hvalue now
        for(int j=0; j<ntagCompact; j++){
            if(cabs(ValuesCompact[j]) < 0.0001) continue;
            IndexI[nsparse] = i;
            IndexJ[nsparse] = BinarySearch(Hsp, Hsize, StatesCompact[j]);
            Hvalue[nsparse] = ValuesCompact[j];
            if(IndexI[nsparse] == IndexJ[nsparse]) Hvalue[nsparse] += 0;
            nsparse++;

            //printf("   Hami matrix %d,%d,%f+i%f\n", i, BinarySearch(Hsp, Hsize, States[j]), \
                   creal(Values[j]), cimag(Values[j]));

            if(BinarySearch(Hsp, Hsize, StatesCompact[j]) == -1){
                printf(" StatesCompact[j] print as %s and tag %d\n", \
                         byte_to_binary(StatesCompact[j]), j);
            }
   
        }
        if(i%10000 == 0) printf("   Matrix generates at %d ntag %d nsparse %d\n", i, ntagCompact, nsparse);
    }

    return nsparse;

}
        /* d = ddMulti(statupdn, U_rest, N)
        ! for k in d:
        !     if k in tdict:
         */


double dot_product_vector(double complex v[], int n)
{
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += creal(v[i])*creal(v[i]) + cimag(v[i])*cimag(v[i]);
    return result;
}

double dot_product(double complex v)
{
    return creal(v)*creal(v) + cimag(v)*cimag(v);
}


void averageElectronCore(int N, int Hsize, unsigned long int *Hsp, double complex *v, double *n){
    //n = np.zeros((2*N))
    //Hsize = Hsp.size
    int N_eg = 2;
    int N_t2g = 3;
    int N_d = 5;
    for(int i=0; i<Hsize; i++){
        unsigned long int statupdn = Hsp[i];
        for(int j=0; j<(2*N); j++){
            unsigned long int oj = 1;
            oj = oj << j;
            if(statupdn & (oj)){
                n[j] += dot_product(v[i]);
                //n[j] += np.real(np.vdot(v[i], v[i]))
            }
        }
    }
    //for(int j=0; j<N; j++){
    //    double t = (n[j] + n[j+N])/2.0;
    //    n[j] = t;
    //    n[j+N] = t;
    //}
    //double tsum = 0.0;
    //for(int j=0; j<N_eg; j++) tsum += n[j];
    //tsum = tsum/2.0;
    ////double tsum = sum(n[:N_eg]) / N_eg
    //for(int i=0; i<N_eg; i++){
    //    n[i] = tsum;
    //    n[i+N] = tsum;    
    //}
    ////for i in range(N_eg) + range(N, N+N_eg):
    ////    n[i] = tsum

    //tsum = 0.0;
    //for(int j=N_eg; j<N_d; j++) tsum += n[j];
    //tsum = tsum/3.0;
    ////tsum = sum(n[N_eg:N_d]) / N_t2g
    //for(int i=N_eg; i<N_d; i++){
    //    n[i] = tsum;
    //    n[i+N] = tsum;    
    //}
    ////for i in range(N_eg, N_d) + range(N + N_eg, N+N_d):
    ////    n[i] = tsum

    //double modv = dot_product_vector(v, Hsize);
    //for(int i=0; i<2*N; i++){
    //    n[i] = n[i]/modv;
    //}
    ////n = n/modv
    return;
}

void GenOStatesCoreSpinUp(double complex *H_0, int Hsize_src, unsigned long int *Hsp_src, \
                double complex *H_e, int Hsize_dst, unsigned long int *Hsp_dst, \
                double *Gn, int N){

    //for(int i=0; i<20; i++) printf("   test H_0 %i %f +i%f\n", i, creal(H_0[i]), cimag(H_0[i]));


    int N_core = 3;
    int N_3d = 5;
    //H_e = np.zeros((Hsp_dst.size), dtype = np.complex64)

    for(int jj=0; jj<Hsize_src; jj++){
        unsigned long int statup = Hsp_src[jj] / (1<<N);
        unsigned long int statdn = Hsp_src[jj] % (1<<N);
        for(int popo=N-N_core; popo < N; popo++){
            //#if (.not.BTEST(statdn,popo)): cycle   !c2p, dn
            for(int toto=0; toto < N_3d; toto++){
                unsigned long int maskpo = 1 << popo;
                unsigned long int maskto = 1 << toto;
                if(statup & maskpo){
                    unsigned long int istatup = statup & (~maskpo);
                    int tempsign = pow(-1, sumEveryBit(statup >> (popo+1)));
                    if(!(istatup & maskto)){
                        unsigned long int fstatup = istatup | maskto;
                        tempsign *= pow(-1, sumEveryBit(istatup >> (toto+1)));                    

                        unsigned long int fstatupdn = fstatup;
                        fstatupdn = (fstatupdn << N) + statdn;

                        //int l = BinarySearch(Hsp_dst,Hsize_dst,fstatup*(1<<N)+statdn);
                        int l = BinarySearch(Hsp_dst,Hsize_dst,fstatupdn);

                        if(l == -1) printf("GenOStates is wrong!\n");
  
                        H_e[l] += Gn[toto*3+ popo - (N-N_core)]*tempsign*H_0[jj];
                    }
                }
            }
        }
    }
    return;
}

void GenODaggerStatesCoreSpinUpDown(double complex *H_0, int Hsize_src, unsigned long int *Hsp_src, \
                double complex *H_e, int Hsize_dst, unsigned long int *Hsp_dst, \
                double *Gn, int N){

    //for(int i=0; i<20; i++) printf("   test H_0 %i %f +i%f\n", i, creal(H_0[i]), cimag(H_0[i]));


    int N_core = 3;
    int N_3d = 5;
    //H_e = np.zeros((Hsp_dst.size), dtype = np.complex64)

    for(int jj=0; jj<Hsize_src; jj++){
        unsigned long int statup = Hsp_src[jj] / (1<<N);
        unsigned long int statdn = Hsp_src[jj] % (1<<N);
        for(int toto=N-N_core; toto < N; toto++){
            //#if (.not.BTEST(statdn,popo)): cycle   !c2p, dn
            for(int popo=0; popo < N_3d; popo++){
                unsigned long int maskpo = 1 << popo;
                unsigned long int maskto = 1 << toto;
                if(statup & maskpo){
                    unsigned long int istatup = statup & (~maskpo);
                    int tempsign = pow(-1, sumEveryBit(statup >> (popo+1)));
                    if(!(istatup & maskto)){
                        unsigned long int fstatup = istatup | maskto;
                        tempsign *= pow(-1, sumEveryBit(istatup >> (toto+1)));                    

                        unsigned long int fstatupdn = fstatup;
                        fstatupdn = (fstatupdn << N) + statdn;
                        //int l = BinarySearch(Hsp_dst,Hsize_dst,fstatup*(1<<N)+statdn);
                        int l = BinarySearch(Hsp_dst,Hsize_dst,fstatupdn);

                        if(l == -1) {
                          //  printf("    GenODaggerStates is wrong!\n");
                          //  printf("    %d, %d\n", Hsp_src[jj], fstatup*(1<<N)+statdn);
                            continue;
                        }
  
                        H_e[l] += Gn[popo*3+ toto - (N-N_core)]*tempsign*H_0[jj];
                    }
                }
                if(statdn & maskpo){
                    unsigned long int istatdn = statdn & (~maskpo);
                    int tempsign = pow(-1, sumEveryBit(statdn >> (popo+1)));
                    if(!(istatdn & maskto)){
                        unsigned long int fstatdn = istatdn | maskto;
                        tempsign *= pow(-1, sumEveryBit(istatdn >> (toto+1)));                    

                        unsigned long int fstatupdn = statup;
                        fstatupdn = (fstatupdn << N) + fstatdn;
                        //int l = BinarySearch(Hsp_dst,Hsize_dst,statup*(1<<N)+fstatdn);
                        int l = BinarySearch(Hsp_dst,Hsize_dst,fstatupdn);

                        if(l == -1) {
                          //  printf("    GenODaggerStates is wrong!\n");
                          //  printf("    %d, %d\n", Hsp_src[jj], fstatup*(1<<N)+statdn);
                            continue;
                        }
  
                        H_e[l] += Gn[popo*3+ toto - (N-N_core)]*tempsign*H_0[jj];
                    }
                }
            }
        }
    }
    return;
}



void GenODaggerStatesCoreSpinUp(double complex *H_0, int Hsize_src, unsigned long int *Hsp_src, \
                double complex *H_e, int Hsize_dst, unsigned long int *Hsp_dst, \
                double *Gn, int N){

    //for(int i=0; i<20; i++) printf("   test H_0 %i %f +i%f\n", i, creal(H_0[i]), cimag(H_0[i]));


    int N_core = 3;
    int N_3d = 5;
    //H_e = np.zeros((Hsp_dst.size), dtype = np.complex64)

    for(int jj=0; jj<Hsize_src; jj++){
        unsigned long int statup = Hsp_src[jj] / (1<<N);
        unsigned long int statdn = Hsp_src[jj] % (1<<N);
        for(int toto=N-N_core; toto < N; toto++){
            //#if (.not.BTEST(statdn,popo)): cycle   !c2p, dn
            for(int popo=0; popo < N_3d; popo++){
                unsigned long int maskpo = 1 << popo;
                unsigned long int maskto = 1 << toto;
                if(statup & maskpo){
                    unsigned long int istatup = statup & (~maskpo);
                    int tempsign = pow(-1, sumEveryBit(statup >> (popo+1)));
                    if(!(istatup & maskto)){
                        unsigned long int fstatup = istatup | maskto;
                        tempsign *= pow(-1, sumEveryBit(istatup >> (toto+1)));                    

                        unsigned long int fstatupdn = fstatup;
                        fstatupdn = (fstatupdn << N) + statdn;
                        //int l = BinarySearch(Hsp_dst,Hsize_dst,fstatup*(1<<N)+statdn);
                        int l = BinarySearch(Hsp_dst,Hsize_dst,fstatupdn);

                        if(l == -1) {
                          //  printf("    GenODaggerStates is wrong!\n");
                          //  printf("    %d, %d\n", Hsp_src[jj], fstatup*(1<<N)+statdn);
                            continue;
                        }
  
                        H_e[l] += Gn[popo*3+ toto - (N-N_core)]*tempsign*H_0[jj];
                    }
                }
            }
        }
    }
    return;
}



void GenODaggerStatesCoreSpinDown(double complex *H_0, int Hsize_src, unsigned long int *Hsp_src, \
                double complex *H_e, int Hsize_dst, unsigned long int *Hsp_dst, \
                double *Gn, int N){

    //for(int i=0; i<20; i++) printf("   test H_0 %i %f +i%f\n", i, creal(H_0[i]), cimag(H_0[i]));


    int N_core = 3;
    int N_3d = 5;
    //H_e = np.zeros((Hsp_dst.size), dtype = np.complex64)

    for(int jj=0; jj<Hsize_src; jj++){
        unsigned long int statup = Hsp_src[jj] / (1<<N);
        unsigned long int statdn = Hsp_src[jj] % (1<<N);
        for(int toto=N-N_core; toto < N; toto++){
            //#if (.not.BTEST(statdn,popo)): cycle   !c2p, dn
            for(int popo=0; popo < N_3d; popo++){
                unsigned long int maskpo = 1 << popo;
                unsigned long int maskto = 1 << toto;
                if(statdn & maskpo){
                    unsigned long int istatdn = statdn & (~maskpo);
                    int tempsign = pow(-1, sumEveryBit(statdn >> (popo+1)));
                    if(!(istatdn & maskto)){
                        unsigned long int fstatdn = istatdn | maskto;
                        tempsign *= pow(-1, sumEveryBit(istatdn >> (toto+1)));                    

                        unsigned long int fstatupdn = statup;
                        fstatupdn = (fstatupdn << N) + fstatdn;
                        //int l = BinarySearch(Hsp_dst,Hsize_dst,statup*(1<<N)+fstatdn);
                        int l = BinarySearch(Hsp_dst,Hsize_dst,fstatupdn);

                        if(l == -1) {
                          //  printf("    GenODaggerStates is wrong!\n");
                          //  printf("    %d, %d\n", Hsp_src[jj], fstatup*(1<<N)+statdn);
                            continue;
                        }
  
                        H_e[l] += Gn[popo*3+ toto - (N-N_core)]*tempsign*H_0[jj];
                    }
                }
            }
        }
    }
    return;
}



