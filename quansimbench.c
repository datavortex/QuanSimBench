/////////////////////////////////////////////////////////////////////////////
//  Quantum Factorization Simulation as a Benchmark for HPC
//  Verifies that the area under the peaks of the Quantum Fourier Transform
//  of delta(2^x mod n,1) is larger than 1/2, where n=p*q is an
//  integer that satisfies n^2<=2^QUBITS<2n^2 and maximizes the period r of 2^x mod n with r even and 2^(r/2)~=-1 mod n.
//  It is a simplification of Shor's factorization algorithm
//  (c) Santiago Ignacio Betelu, Denton 2018
//  Thanks Datavortex Technologies for providing the hardware and research support for developing this benchmark.
//  mpicc -Ofast quansimbench.c -o quansimbench -lm -Wall
//  sbatch quansimbench.batch
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
#include <inttypes.h>
#include <time.h>
#include <mpi.h>
#ifdef _OPENMP
# include <omp.h>
#endif

#define VERSION "1.0"
#ifndef MINQUBITS
# define MINQUBITS 9
#endif
#if MINQUBITS < 9
# error MINQUBITS must be at least 9
#endif
#ifndef MAXQUBITS
# define MAXQUBITS 60
#endif

static complex float *c=NULL, *buffer=NULL; // quantum amplitudes
static int64_t QUBITS,N,NBUFFERS,NODEBITS;
static int BUFFERSIZE,nranks,inode;
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
    int64_t x,y,mask1,mask2,q,chunk;
    int node,b,tag;
    complex float aux;
    static MPI_Request reqsend[1024], reqrecv[1024];
    //
    if(qubit< QUBITS-NODEBITS){
        mask1= (0xFFFFFFFFFFFFFFFFll<<qubit);  // to avoid branching and half of memory accesses
        mask2=  ~mask1;
        mask1= (mask1<<1);
#pragma omp parallel for private(x,y,aux)
        for(q=0;q<N/2/nranks;q++){
            x= ((q<<1)&mask1) | (q&mask2); // 64 bit index with 0 on the qubit'th position
            y= x|(1ll<<qubit);             //        index with 1 on the qubit'th position
            aux=  (c[x]-c[y])*M_SQRT1_2;
            c[x]= (c[x]+c[y])*M_SQRT1_2;
            c[y]=aux;
        }
    }else{
        node= inode^(1ULL<<(qubit-(QUBITS-NODEBITS)));
        tag=0;
        for(chunk=0; chunk<N/nranks; chunk=chunk+NBUFFERS*BUFFERSIZE){
            for(b=0;b<NBUFFERS;b++){
                tag= tag+1;
                MPI_Irecv(  &buffer[b*BUFFERSIZE], BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqrecv[b]);
                MPI_Isend( &c[chunk+b*BUFFERSIZE], BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqsend[b]);
            }
            for(b=0;b<NBUFFERS;b++){
                MPI_Wait(&reqsend[b],MPI_STATUS_IGNORE);
                MPI_Wait(&reqrecv[b],MPI_STATUS_IGNORE);
                if( inode&(1ll<<(qubit-(QUBITS-NODEBITS))) ){
#pragma omp parallel for
                    for(q=0; q<BUFFERSIZE; q++){
                        c[chunk+q+b*BUFFERSIZE]= -(c[chunk+q+b*BUFFERSIZE]-buffer[b*BUFFERSIZE+q])*M_SQRT1_2;
                    }
                }else{
#pragma omp parallel for
                    for(q=0; q<BUFFERSIZE; q++){
                        c[chunk+q+b*BUFFERSIZE]=  (c[chunk+q+b*BUFFERSIZE]+buffer[b*BUFFERSIZE+q])*M_SQRT1_2;
                    }
                }
            }
        }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////
static void SWAP(int64_t qubit1, int64_t qubit2){  // SWAP between qubit1 and qubit2, qubit1!=quibit2
    int64_t x,y,b1,b2,chunk,q;
    int node,b,tag;
    complex float aux;
    static MPI_Request reqsend[1024], reqrecv[1024];
    //
    if(qubit1>qubit2){ // sort qubit1 < qubit2
        q=qubit1;
        qubit1=qubit2;
        qubit2=q;
    }
    if(qubit2<QUBITS-NODEBITS && qubit1<QUBITS-NODEBITS){
#pragma omp parallel for private(x,y,b1,b2,aux)
        for(q=0;q<N/nranks;q++){
            x= q+ 0*inode*(N/nranks);  // 0* because affects only lower qubits
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
    }else if(qubit1 >= QUBITS-NODEBITS && qubit2 >= QUBITS-NODEBITS) { // in this case swap all array alements with another node
        x=  inode*(N/nranks);
        b1= (x>>qubit1)&1ll;
        b2= (x>>qubit2)&1ll;
        if( b1!=b2 ){
            node= inode^(1<<(qubit2-(QUBITS-NODEBITS)));  // here qubit2 >= QUBITS-NODEBITS for sure
            node=  node^(1<<(qubit1-(QUBITS-NODEBITS)));
            tag=0;
            for(chunk=0; chunk<N/nranks; chunk=chunk+NBUFFERS*BUFFERSIZE){
                for(b=0;b<NBUFFERS;b++){
                    tag=tag+1;
                    MPI_Irecv( &buffer[b*BUFFERSIZE],  BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqrecv[b]);
                    MPI_Isend( &c[b*BUFFERSIZE+chunk], BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqsend[b]);
                }
                for(b=0;b<NBUFFERS;b++){
                    MPI_Wait(&reqsend[b],MPI_STATUS_IGNORE);
                    MPI_Wait(&reqrecv[b],MPI_STATUS_IGNORE);
#pragma omp parallel for
                    for(q=0; q<BUFFERSIZE; q++){
                        c[b*BUFFERSIZE+chunk+q]= buffer[b*BUFFERSIZE+q];
                    }
                }
            }
        }
    }else{  // qubit1 inside same node but qubit2 in another node
        node= inode^(1<<(qubit2-(QUBITS-NODEBITS)));  // here qubit2 >= QUBITS-NODEBITS for sure
        x= node*(N/nranks);
        b2= (x>>qubit2)&1ll;
        tag=0;
        for(chunk=0; chunk<N/nranks; chunk=chunk+NBUFFERS*BUFFERSIZE){
            for(b=0;b<NBUFFERS;b++){
                tag=tag+1;
                MPI_Irecv( &buffer[b*BUFFERSIZE],  BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqrecv[b]);
                MPI_Isend( &c[b*BUFFERSIZE+chunk], BUFFERSIZE, MPI_COMPLEX, node, tag, MPI_COMM_WORLD, &reqsend[b]);
            }
            for(b=0;b<NBUFFERS;b++){
                MPI_Wait(&reqsend[b],MPI_STATUS_IGNORE);
                MPI_Wait(&reqrecv[b],MPI_STATUS_IGNORE);
#pragma omp parallel for private(x,y,b1)
                for(q=0; q<BUFFERSIZE; q=q+1){
                    x= b*BUFFERSIZE+chunk+q; // received register
                    b1= (x>>qubit1)&1ll;
                    y= (b*BUFFERSIZE+chunk+q)^(1ll<<qubit1);  // guaranteed y<x
                    if( b1!=b2 ) c[y]= buffer[b*BUFFERSIZE+q];
                }
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
    int64_t x,q,b1,b2,k,qubit2;
#pragma omp parallel for private(x,b1,b2,k,qubit2)
    for(q=0;q<N/nranks;q++){
        x= q+inode*(N/nranks);
        b1= ((x>>qubit1)&1ll);
        if( b1 == 0 )
            continue;
        for(k=1;k<=nq;k++){
            qubit2=qubit1-k;
            if(qubit2>=0){
                b2= ((x>>qubit2)&1ll);
                if( b2 == 0 )
                    continue;
                c[q]=c[q]*expphase[k];
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
static void init_quantum_state(int64_t n){  // initial state is | z, 2^z mod n > collapsed by a measurement of second register with outcome 1
    int64_t x,l,z,zp0,period;
    double s0,s;   // don't change to float

    // find the starting point and period of (2^x mod n)
    x= inode*(N/nranks);
    l= powmod(2,x,n); // l is the value of (2^x mod n)
    zp0= -1;
    period= 1;
    for(z=0;z<N/nranks;z++){
        if (l==1){
            if (zp0==-1)
                zp0= z;
            else {
                period=z-zp0;
                break;
            }
        }
        l= (2*l)%n;  // fast computation of (2^x mod n)
    }

    // set every period'th value to 1.0 and the rest to 0.0
#pragma omp parallel for private(z)
    for(z=0;z<N/nranks;z++){
        c[z]=0.0;
    }
#pragma omp parallel for private(z) firstprivate(zp0,period)
    for(z=zp0;z<N/nranks;z+=period){
        c[z]=1.0;
    }

    // normalize the initial condition
    s0=0.0, s=0.0; // for normalization
#pragma omp parallel for private(z) reduction(+:s0)
    for(z=0;z<N/nranks;z++){
        s0=s0+ cabsf(c[z]*conjf(c[z]));
    }
    MPI_Allreduce(&s0,&s, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    s=1.0/sqrt(s);
#pragma omp parallel for private(z)
    for(z=0;z<N/nranks;z++){
        c[z]= c[z]*s;
    }
}

int main(int argc, char **argv){
    int64_t aux,nphase,max_nphase,n,mulperiod,peaknumber,z,q,numstates,npeaks,predictedx;
    struct timespec tim0,tim1;
    double timeperstate,timeqft,prob,prob0,peakspacing; // don't change to float
    char texfactors[32];
    complex float *expphase=NULL;
    int retval = EXIT_FAILURE;  // assume failure

    // largest integers that can be factored with Shor's algoritm with register size 'qubits'
    // n[qubits]= factor1[qubits]*factor2[qubits]   2^qubits <= n^2 < 2^{qubits+1}, qubits>=9
    int64_t factor1[61]={0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 5, 7, 5, 11, 11, 5, 19, 23, 19, 23, 29, 47, 29, 29, 47, 71, 83, 79, 103, 149, 101, 149, 269, 167, 479, 479, 367, 859, 563, 1039, 947, 1307, 2027, 2039, 2357, 2237, 3917, 4127, 4813, 6173, 6029, 7243, 10357, 12757, 11399, 19427, 20771, 24847, 27779};

    int64_t factor2[61]={0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 13, 11, 11, 23, 13, 23, 71, 23, 29, 53, 61, 67, 61, 139, 199, 173, 163, 197, 293, 317, 311, 647, 619, 487, 1109, 547, 773, 1427, 863, 1861, 1427, 2213, 2269, 2069, 2909, 3559, 5303, 4283, 5749, 6971, 7687, 11131, 13103, 12959, 14879, 23549, 19541, 25847, 30557, 38653};

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&inode);

    NODEBITS=0;
    for(aux=1; aux<nranks; aux<<=1)
        NODEBITS++;
    if(aux!=nranks){
        if(inode==0) fprintf(stderr,"ERROR: Number of nodes has to be a power of 2\n");
        goto fin;
    }
    if(inode==0){
        printf("QuanSimBench version %s\n",VERSION);
#ifdef QSB_MPI_STUBS
        printf("MPI ranks: N/A\n");
#else
        printf("MPI ranks: %d\n", nranks);
#endif
#ifdef _OPENMP
        printf("OpenMP threads per rank: %d\n", omp_get_max_threads());
#else
        printf("OpenMP threads per rank: N/A\n");
#endif
        printf("\n");
        printf("Qubits      Factors    Probability         Time    Coeffs/s  Pass\n");
    }

    // pre-initialize the phase exponentials
    max_nphase= 1 + (int64_t)log2(1.0*MAXQUBITS);
    if( posix_memalign((void **)&expphase, sizeof(complex float), (1 + max_nphase)*sizeof(complex float)) != 0 ){     // re-allocate phase exponentials
        if(inode==0) fprintf(stderr,"Ending due to allocation error\n");
        goto fin;
    }
    init_expphase(max_nphase,expphase);

    // iterate over number of qubits
    for(QUBITS=MINQUBITS; QUBITS<=MAXQUBITS; QUBITS++){
        N= (1ll<<QUBITS); // state vector size
        if( N<nranks ) continue;  // too many nodes for small N

        BUFFERSIZE= (1ll<<18);  // number of complex numbers used in chunk of communication
        NBUFFERS=4; // must be a power of 2 to simplify code, and <=1024 (which is too large)
        if( NBUFFERS> N/nranks/BUFFERSIZE ) NBUFFERS= N/nranks/BUFFERSIZE;
        if( NBUFFERS<1 ) NBUFFERS=1;
        if( BUFFERSIZE>N/nranks/NBUFFERS ) BUFFERSIZE=N/nranks/NBUFFERS;
        if( N%(nranks*BUFFERSIZE*NBUFFERS)!=0){
            if(inode==0) fprintf(stderr,"ERROR: nranks*BUFFERSIZE must divide N %" PRId64 "\n",nranks*BUFFERSIZE*NBUFFERS);
            goto fin;
        }
        free(c);
        if( posix_memalign((void **)&c, sizeof(complex float), (N/nranks)*sizeof(complex float)) != 0 ){     // re-allocate double float amplitudes
            if(inode==0) fprintf(stderr,"Ending due to allocation error\n");
            goto fin;
        }
        free(buffer);
        if( posix_memalign((void **)&buffer, sizeof(complex float), NBUFFERS*BUFFERSIZE*sizeof(complex float)) != 0 ){    // for communication
            if(inode==0) fprintf(stderr,"Ending due to allocation error\n");
            goto fin;
        }

        n= factor1[QUBITS]*factor2[QUBITS]; // number to factor
        mulperiod= (factor1[QUBITS]-1)*(factor2[QUBITS]-1); // Euler totient function is multiple of period of (2^x mod n)
        peakspacing= 1.0*N/mulperiod; // so the space between peaks in the spectrum is a multiple of this

        if(n*n>=N){ // n^2<N for Shor's algorithm validity
            fprintf(stderr,"Error n*n>=N\n");
            goto fin;
        }

        init_quantum_state(n);
        nphase= 1 + (int64_t)log2(1.0*QUBITS);   // number of phases in each step of Approximate Quantum Fourier Transform
        clock_gettime(CLOCK_REALTIME,&tim0);  // only time AQFT
        // the Approximate Quantum Fourier Transform
        numstates=0;
        for(q=QUBITS-1;q>=0; q--){
            H(q);
            CPN(q,nphase,expphase); // all nphase phases folded into a single call
            numstates=numstates+1+min(q,nphase);
        }
        for(q=0;q<QUBITS/2;q++){
            SWAP(q,QUBITS-q-1);
            numstates=numstates+1;
        }
        // end AQFT

        clock_gettime(CLOCK_REALTIME,&tim1);
        timeqft= 1.0*(tim1.tv_sec-tim0.tv_sec)+1.e-9*(tim1.tv_nsec-tim0.tv_nsec); // time of QFT in seconds
        timeperstate= (N*numstates)/timeqft;

        // compute probability that the solution is a multiple of peakspacing
        prob0=0.0;
        npeaks= mulperiod;
        for(peaknumber= inode*npeaks/nranks; peaknumber<=(inode+1)*npeaks/nranks; peaknumber++){  // note that this lists << N peaks
            if(peaknumber>0) {
                predictedx= peaknumber*peakspacing +0.5; // state number x where a peak may occur, add 0.5 to round to nearest
                z= predictedx -N/nranks*inode;           // convert to int and reduce to interval in this node
                if(z>=0 && z<N/nranks) prob0=prob0+cabsf(c[z]*conjf(c[z]));  // resulting area under theoretical peaknumber
            }
        }
        MPI_Allreduce(&prob0,&prob, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(inode==0){
            sprintf(texfactors,"%" PRId64 "*%" PRId64,factor1[QUBITS], factor2[QUBITS]);
            printf("%6" PRId64 " %12s  %13.6f   %10.4e  %10.4e  %4s\n", QUBITS, texfactors, prob, timeqft, timeperstate, prob > 0.5 ? "yes" : "no");
            fflush(stdout);
        }
    }
    retval = EXIT_SUCCESS;

 fin:
    free(expphase);
    free(buffer);
    free(c);
    MPI_Finalize();
    return retval;
}
////////////////////////////////////////////////////////////////////////////////
