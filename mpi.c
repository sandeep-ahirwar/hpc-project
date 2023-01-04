#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2
#define ORGS 10000
#define GENES 100
#define ALLELES 4
#define MUT 1000

char **curG, **nextG, *mod;
int *f, totF, Eval(), Run();
void Mem(), Init(), Gen();

// cluster /parallel code
int main(int argc, char *argv[]) {

    int numtasks, taskid, numworkers, source, dest,
    mtype, rows, averow, extra, offset, i, j, k, rc;
    long double a[ALLELES][ALLELES], b[ALLELES][ALLELES], c[ALLELES][ALLELES];
    MPI_Status status;
    double start, end;
    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    if (numtasks < 2) {
        printf("Need at least two MPI tasks\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }
    numworkers = numtasks - 1;
    char pro_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(pro_name, &name_len);
    if (taskid == MASTER) {
        for (i = 0; i < ALLELES; i++)
            for (j = 0; j < ALLELES; j++)
                a[i][j] = i + j;

        for (i = 0; i < ALLELES; i++)
            for (j = 0; j < ALLELES; j++) {
            b[i][j] = rand()%MUT;
            c[i][j] = rand()%MUT;
        }
        averow = ALLELES / numworkers;
        extra = ALLELES % numworkers;
        offset = 0;
        mtype = FROM_MASTER;
        for (dest = 1; dest <= numworkers; dest++) {
            rows = (dest <= extra) ? averow + 1 : averow;
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&a[offset][0], rows * ALLELES, MPI_LONG_DOUBLE, dest, mtype,
            MPI_COMM_WORLD);
            MPI_Send(&b, ALLELES * ALLELES, MPI_LONG_DOUBLE, dest, mtype,
            MPI_COMM_WORLD);
            offset = offset + rows;
        }

        mtype = FROM_WORKER;
        for (i = 1; i <= numworkers; i++) {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&c[offset][0], rows * ALLELES, MPI_LONG_DOUBLE, source, mtype,
            MPI_COMM_WORLD, &status);
        }
        Mem();
        printf("The final generation was: %d\n", Run());

        end = MPI_Wtime(); 
        printf("\nTime= %f", end - start);
    }

    if (taskid > MASTER) {
        mtype = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, rows * ALLELES, MPI_LONG_DOUBLE, MASTER, mtype,
        MPI_COMM_WORLD, &status);
        MPI_Recv(&b, ALLELES * ALLELES, MPI_LONG_DOUBLE, MASTER, mtype,
        MPI_COMM_WORLD, &status);

        for (k = 0; k < ALLELES; k++)
            for (i = 0; i < rows; i++) {
                c[i][k] = 0.0;
                for (j = 0; j < ALLELES; j++)
                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
            }

        mtype = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&c, rows * ALLELES, MPI_LONG_DOUBLE, MASTER, mtype,
        MPI_COMM_WORLD);
    }
    MPI_Finalize();
}


// Initialize Memory Function
void Mem(){
	int o;
	curG=(char**)malloc(sizeof(char*)*ORGS);
	nextG=(char**)malloc(sizeof(char*)*ORGS);
	mod=(char*)malloc(sizeof(char)*GENES);
	f=(int*)malloc(sizeof(int)*ORGS);

	for(o=0; o<ORGS; ++o){
		curG[o]=(char*)malloc(sizeof(char)*GENES);
		nextG[o]=(char*)malloc(sizeof(char)*GENES);
	}
}

// Run Function
int Run(){
	int gen=0;
	Init();
	while(++gen) {
		if(Eval()) {
			return gen;
		}
		Gen();
	}
}

// Initialization Function
void Init(){
	int o, g;
	for(o=0; o<ORGS; o++) {
		for(g=0; g<GENES; g++) {
			curG[o][g]=rand()%ALLELES;
		}
	}	
	for(g=0; g<GENES; g++){
		mod[g]=rand()%ALLELES;
	}
}

// Evaluation Function
int Eval(){
	int o, g, curF;
	for(totF=0, o=0; o<ORGS; ++o) {
		for(curF=0, g=0; g<GENES; ++g) {
			if(curG[o][g]==mod[g]) {
				++curF;
			}
		}
		if(curF==GENES) { 
			return 1;
		}
		totF += f[o]=curF;
	}
	return 0;
}

// Generation Function
void Gen(){
	int o, g, p1, p2, cp;
	#pragma omp parallel for
	for(o=0;o<ORGS;++o) {
		int tot = 0, pt = rand()%(totF+1);
		
		for(int o1 = 0; o1 < ORGS; ++o1){
			if((tot += f[o1]) >= pt) {
				p1 = o1;
				break;
			}
		}
		tot = 0; pt = rand() % (1+totF);
		for(int o2 = 0; o2 < ORGS; ++o2){
			if((tot += f[o2]) >= pt) {
				p2 = o2;
				break;
			}
		}
		
		for(p1, p2, cp = rand() % GENES, g = 0; g < GENES; ++g) {
			if(rand()%MUT) {
				if(g<cp) {
					nextG[o][g] = curG[p1][g];
				} else {
					nextG[o][g] = curG[p2][g];
				}
			} else {
				nextG[o][g] = rand()%ALLELES;
			}
		} 
	}

	for(o=0; o<ORGS; ++o) {
		for(g=0; g<GENES; ++g) { 
			curG[o][g]=nextG[o][g];
		}
	}
}
