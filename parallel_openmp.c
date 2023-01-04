#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define ORGS 10000
#define GENES 50
#define ALLELES 4
#define MUT 1000

char **curG, **nextG, *mod;
int *f, totF, Eval(), Run();
void Mem(), Init(), Gen();

int main(){
  float pTime = omp_get_wtime();
  Mem();
  printf("The final generation was: %d\n", Run());
  float cTime = omp_get_wtime();
  printf("%f",cTime-pTime);}

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
