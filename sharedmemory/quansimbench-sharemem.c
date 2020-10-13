/////////////////////////////////////////////////////////////////////////////
//  Quantum Factorization Simulation as a Benchmark for HPC
//  Shared memory version
//  Verifies that the area under the peaks of the Quantum Fourier Transform
//  of delta(2^x mod n,1) is larger than 1/2, where n=p*q is an
//  integer that satisfies n^2<=2^QUBITS<2n^2 and maximizes the period r of 2^x mod n with r even and 2^(r/2)~=-1 mod n.
//  It is a simplification of Shor's factorization algorithm
//  (c) Santiago Ignacio Betelu, Denton 2018
//  Thanks Datavortex Technologies, UNT/HPC, TACC and LANL for providing the hardware and research support for developing this benchmark.
//  gcc -Ofast quansimbench-sharemem.c -o quansimbench -lm -Wall -fopenmp
//  export OMP_NUM_THREADS=32 
//  ./quansimbench
//    _______                    ______ _       ______                    _
//   (_______)                  / _____|_)     (____  \                  | |
//    _    _  _   _ _____ ____ ( (____  _ ____  ____)  )_____ ____   ____| |__
//   | |  | || | | (____ |  _ \ \____ \| |    \|  __  (| ___ |  _ \ / ___)  _ |
//   | |__| || |_| / ___ | | | |_____) ) | | | | |__)  ) ____| | | ( (___| | | |
//    \______)____/\_____|_| |_(______/|_|_|_|_|______/|_____)_| |_|\____)_| |_|
//
////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>

#define MINQUBITS 9
#define MAXQUBITS 60

complex float *c=NULL; // quantum amplitudes
int64_t QUBITS,N;
/////////////////////////////////////////////////////////////////////////////
//  Quantum numstates addressing example with 4 nodes.
//  1- The QUBITS-NODEBITS least significant bits can be swapped within each node.
//  2- The NODEBITS most significant digits is node number
//     NODE
//       |  Local bits
//  c0   00 000           c16  10 000
//  c1   00 001           c17  10 001
//  c2   00 010           c18  10 010
//  c3   00 011  N0       c19  10 011  N2
//  c4   00 100           c20  10 100
//  c5   00 101           c21  10 101
//  c6   00 110           c22  10 110
//  c7   00 111           c23  10 111
//       ......                ......
//  c8   01 000           c24  11 000
//  c9   01 001           c25  11 001
//  c10  01 010           c26  11 010
//  c11  01 011  N1       c27  11 011  N3
//  c12  01 100           c28  11 100
//  c13  01 101           c29  11 101
//  c14  01 110           c30  11 110
//  c15  01 111           c31  11 111
//       ......                ......
//////////////////////////////////////////////////////////////////////////////
//  H= | 1  1 |
//     | 1 -1 | /sqrt(2)
static void H(int64_t qubit){  // Hadamard gate acting on qubit
    int64_t x,y,mask1,mask2,q;
    complex float aux;
    mask1= (0xFFFFFFFFFFFFFFFFll<<qubit);  // to avoid branching and half of memory accesses
    mask2=  ~mask1;
    mask1= (mask1<<1);
    #pragma omp parallel for private(x,y,aux)
    for(q=0;q<N/2;q++){
            x= ((q<<1)&mask1) | (q&mask2); // 64 bit index with 0 on the qubit'th position
            y= x|(1ll<<qubit);             //        index with 1 on the qubit'th position
            aux=  (c[x]-c[y])*M_SQRT1_2;
            c[x]= (c[x]+c[y])*M_SQRT1_2;
            c[y]=aux;
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////
static void SWAP(int64_t qubit1, int64_t qubit2){  // SWAP between qubit1 and qubit2, qubit1!=quibit2
    int64_t x,y,b1,b2,q;
    complex float aux;
    //
    if(qubit1>qubit2){ // sort qubit1 < qubit2
        q=qubit1;
        qubit1=qubit2;
        qubit2=q;
    }
    #pragma omp parallel for private(x,y,b1,b2,aux)
    for(x=0;x<N;x++){
            y= (x^(1ll<<qubit1))^(1ll<<qubit2);
            if(y>x){ // to avoid overwriting previously computed
                b1= (x>>qubit1)&1ll;
                b2= (x>>qubit2)&1ll;
                if(b1!=b2){
                    aux= c[x];
                    c[x]=c[y];
                    c[y]=aux;
                }
            }
    }
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
static void init_expphase(int64_t nq,complex float *expphase){   // initialize the phase exponentials
    float phase;
    int64_t k;

    for(k=1;k<=nq;k++){
        phase= M_PI*powf(2.0,-(float)k);
        expphase[k]= cexpf(I*phase);
    }
}

static void CPN(int64_t qubit1, int64_t nq, complex float *expphase){  // PHASE between control qubit1 and qubit+1,2,3,..nq, phase= pi/2^1, pi/2^2,...
    int64_t x,b1,b2,k,qubit2;
    #pragma omp parallel for private(x,b1,b2,k,qubit2)
    for(x=0;x<N;x++){
        b1= ((x>>qubit1)&1ll);
        if( b1 == 1 ){
          for(k=1;k<=nq;k++){
             qubit2=qubit1-k;
             if(qubit2>=0){
                b2= ((x>>qubit2)&1ll);
                if( b2==1 ) c[x]=c[x]*expphase[k];
             }
          }
        }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////
static int64_t min(int64_t x, int64_t y){
    if(x<y) return(x);
    else return(y);
}
//  (a^b) mod n
static int64_t powmod(int64_t a, int64_t b, int64_t n){
    int64_t xq=1,yq=a; // avoid overflow of intermediate results
    while(b>0){
        if(b&1ll) xq=(xq*yq)%n;
        yq=(yq*yq)%n; // squaring the base
        b=b>>1;
    }
    return(xq%n);
}
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){
    int64_t nphase,max_nphase,n,mulperiod,peaknumber,x,l,q,numgates,npeaks,predictedx,threads=1;
    struct timespec tim0,tim1;
    double timeperstate,timeqft,prob,peakspacing,s; // don't change to float
    char texfactors[32];
    complex float *expphase=NULL;

    // largest integers that can be factored with Shor's algoritm with register size 'qubits'
    // n[qubits]= factor1[qubits]*factor2[qubits]   2^qubits <= n^2 < 2^{qubits+1}, qubits>=9
    int64_t factor1[61]={0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 5, 7, 5, 11, 11, 5, 19, 23, 19, 23, 29, 47, 29, 29, 47, 71, 83, 79, 103, 149, 101, 149, 269, 167, 479, 479, 367, 859, 563, 1039, 947, 1307, 2027, 2039, 2357, 2237, 3917, 4127, 4813, 6173, 6029, 7243, 10357, 12757, 11399, 19427, 20771, 24847, 27779};

    int64_t factor2[61]={0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 13, 11, 11, 23, 13, 23, 71, 23, 29, 53, 61, 67, 61, 139, 199, 173, 163, 197, 293, 317, 311, 647, 619, 487, 1109, 547, 773, 1427, 863, 1861, 1427, 2213, 2269, 2069, 2909, 3559, 5303, 4283, 5749, 6971, 7687, 11131, 13103, 12959, 14879, 23549, 19541, 25847, 30557, 38653};

    #if defined(_OPENMP)
        threads= omp_get_max_threads();
        printf("OpenMP threads: %ld\n\n", threads);
    #else
        threads=1;
        printf("Not using OpenMP\n\n");
    #endif
    printf("Qubits      Factors    Probability   Time   Coeffs/s  Coeffs/s/thread  Pass\n");

    // pre-initialize the phase exponentials
    max_nphase= 2+log2(1.0*MAXQUBITS);
    expphase= malloc(max_nphase*sizeof(complex float));
    init_expphase(max_nphase,expphase);

    // iterate over number of qubits
    for(QUBITS=MINQUBITS; QUBITS<=MAXQUBITS; QUBITS++){
        N= (1ll<<QUBITS); // state vector size
        c= realloc(c, N*sizeof(complex float));

        n= factor1[QUBITS]*factor2[QUBITS]; // number to factor
        mulperiod= (factor1[QUBITS]-1)*(factor2[QUBITS]-1); // Euler totient function is multiple of period of (2^x mod n)
        peakspacing= 1.0*N/mulperiod; // so the space between peaks in the spectrum is a multiple of this

        if(n*n>=N){ // n^2<N for Shor's algorithm validity
            printf("Error n*n>=N\n");
            exit(1);
        }

        // initial state is | z, 2^z mod n > collapsed by a measurement of second register with outcome 1
       s=0.0; // for normalization
       #pragma omp parallel for private(l,x) reduction(+:s)
       for(x=0; x<N; x++){
           c[x]=0.0;
           l= powmod(2,x,n);
           if(l==1){
              c[x]=1.0;
              s=s+ 1.0;
           }
       }
       s=1.0/sqrt(s);

       #pragma omp parallel for 
       for(x=0;x<N;x++) c[x]= c[x]*s; // normalize initial condition

       nphase= 1 + (int64_t)log2(1.0*QUBITS);   // number of phases in each step of Approximate Quantum Fourier Transform
       clock_gettime(CLOCK_REALTIME,&tim0);  // only time AQFT
       // the Approximate Quantum Fourier Transform
       numgates=0;
       for(q=QUBITS-1;q>=0; q--){
            H(q);
            CPN(q,nphase,expphase); // all nphase phases folded into a single call
            numgates=numgates+1+min(q,nphase);
       }
       for(q=0;q<QUBITS/2;q++){
            SWAP(q,QUBITS-q-1);
            numgates=numgates+1;
       }
       // end AQFT

       clock_gettime(CLOCK_REALTIME,&tim1);
       timeqft= 1.0*(tim1.tv_sec-tim0.tv_sec)+1.e-9*(tim1.tv_nsec-tim0.tv_nsec); // time of QFT in seconds
       timeperstate= (N*numgates)/timeqft;

       // compute probability that the solution is a multiple of peakspacing
       prob=0.0;
       npeaks= mulperiod;
       for(peaknumber=0 ; peaknumber<=npeaks; peaknumber++){  // note that this lists << N peaks
            if(peaknumber>0) {
                predictedx= peaknumber*peakspacing +0.5; // state number x where a peak may occur, add 0.5 to round to nearest
                if(predictedx>=0 && predictedx<N) prob=prob+creal(c[predictedx]*conjf(c[predictedx]));  // resulting area under theoretical peaknumber
            }
       }
       sprintf(texfactors,"%ld*%ld ",factor1[QUBITS], factor2[QUBITS]);
       printf(" %ld  %12s  %13.6f   %10.4e  %10.4e  %10.4e %4s\n", QUBITS, texfactors, prob, timeqft, timeperstate, timeperstate/threads, prob > 0.5 ? "yes" : "no");
       fflush(stdout);
    }

    free(expphase);
    free(c);
}
////////////////////////////////////////////////////////////////////////////////
