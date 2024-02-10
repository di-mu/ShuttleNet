#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define K 20
#define DIM 3
#define ZWIN 10
#define ADRWIN	20
#define ADRMAX	1
#define SF12DOM	1
#define SFSTAT	5

#define TRAIN	1
#define TEST	2
#define PROB	4
#define ADR 	8
#define HIST	16
#define GPSDB	32
#define ERREP	64
#define STAT	128
#define GPSIN	256
#define KNN 	512
#define RAWTHPT	1024
#define ORGTHPT 2048
#define SFREC	4096
#define NOGPS	8192

#define UNDEF	0xDD
#define ALL 	0xFF
#define NAN 	0x7FFF

#define ABS(x) ((x) > 0 ? (x) : -(x))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define IN(x,a,b) ((x) >= (a) && (x) <= (b))

float VOTMAJ[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; //qualified voting majority
const uint8_t SFPKT[6] = {20, 12, 7, 4, 2, 1};
const float SFTIME[6] = {0.103, 0.170, 0.296, 0.520, 1.013, 2.000};
FILE *pof, *pef, *psf;

//#define KNN_GPS
#ifdef KNN_GPS
const int16_t XMIN[3] = {0, 0, 0};
const int16_t XMAX[3] = {150, 150, 0};
#define TS0 151
#define TS1 151
#define TS2 1
#else
const int16_t XMIN[3] = {-110, -110, -110};
const int16_t XMAX[3] = {-80, -80, -80};
#define TS0 31
#define TS1 31
#define TS2 31
#endif

typedef struct datacell {
	int16_t x[3];
	int16_t gps[2];
	uint8_t y;
	uint8_t sf;
	struct datacell* next;
} DataCell;

typedef struct sfselect {
	int16_t optscs;
	int16_t selscs;
} SFSelect;

DataCell traincell[10000], testcell[500000];
DataCell* table[TS0][TS1][TS2] = {0};

int16_t distance(DataCell* c1, DataCell* c2) {
	int16_t d = 0; uint8_t ii;
	for (ii=0; ii<3; ii+=1) d += ABS(c1->x[ii] - c2->x[ii]);
}

void addlink(DataCell* c, uint16_t mode) {
	uint8_t ii;
	DataCell** phead;
	int16_t *x = (mode & GPSIN) ? c->gps : c->x;
	//printf("%d,%d,%d\n", x[0] - XMIN[0], x[1] - XMIN[1], x[2] - XMIN[2]);
	for (ii=0; ii < DIM; ii+=1) if (!IN(x[ii], XMIN[ii], XMAX[ii])) return;
	phead = &table[x[0] - XMIN[0]][x[1] - XMIN[1]][(DIM == 3) ? x[2] - XMIN[2] : XMIN[2]];
	c->next = *phead;
	*phead = c;
}

int32_t readfile(char* fn, DataCell* c, int32_t maxlen, uint16_t mode) {
	int16_t x0, x1, gps0, gps1, x00 = NAN, x01 = NAN;
	int32_t ii = 0, jj, kk;
	uint8_t y = 0;
	int8_t n_col;
	char strgps[20];
	FILE *pfgps;
	FILE *pfr = fopen(fn, "r");
	while ((n_col = fscanf(pfr, "%hd,%hd,%hd,%hd", &x0, &x1, &gps0, &gps1)) > 0) {
		y |= (x0 != NAN) << (ii % 6);
		if (!SF12DOM && x0 != NAN) { x00 = x0; x01 = x1; }
		if (ii % 6 < 5) { ii+=1; continue; }
		if (SF12DOM) { x00 = x0; x01 = x1; }
		if (x00 == NAN) { ii-=5; goto NEXT; }
		jj = ii / 6; ii += 1;
		if (mode & GPSDB) {
			sprintf(strgps, "gps_db/%d,%d", gps0, gps1);
			pfgps = fopen(strgps, "a");
			fprintf(pfgps, "%02X\n", y);
			fclose(pfgps);
			goto NEXT;
		}
		c[jj].x[0] = x00;
		c[jj].x[1] = (DIM > 1) ? x01 : XMIN[1]; //(x1-XMIN[1])/1+XMIN[1]
		c[jj].x[2] = (DIM > 2) ? 0 : XMIN[2];
		//if (jj==0) printf("read gps %d %d\n", gps0, gps1);
		if (n_col == 4) {
			c[jj].gps[0] = (gps0 != NAN) ? gps0 : c[jj-1].gps[0];
			c[jj].gps[1] = (gps1 != NAN) ? gps1 : c[jj-1].gps[1];
		}
		if (DIM > 2) {
			for (kk=0; kk < MIN(jj+1, ZWIN); kk+=1) c[jj].x[2] += c[jj-kk].x[0];
			c[jj].x[2] /= kk;
		}
		c[jj].sf = (mode & TRAIN) ? ALL : 0;
		if (jj > 0) {
			c[jj-1].y = y;
			if (mode & TRAIN) addlink(c+jj-1, mode);
		}
NEXT:
		y = 0;
		x00 = x01 = NAN;
		if (jj == maxlen) break;
	}
	fclose(pfr);
	return jj;
}

void addxyz(int16_t x0, int16_t x1, int16_t x2, int16_t* pk, uint8_t* y, uint8_t* sf) {
	DataCell* c;
	//printf("%d,%d,%d\n", x0 - XMIN[0], x1 - XMIN[1], x2 - XMIN[2]);
	if (!IN(x0, XMIN[0], XMAX[0]) || 
		!IN(x1, XMIN[1], XMAX[1]) || 
		!IN(x2, XMIN[2], XMAX[2])) return;
	c = table[x0 - XMIN[0]][x1 - XMIN[1]][x2 - XMIN[2]];
	while (c) {
		y[*pk] = c->y;
		sf[*pk] = c->sf;
		(*pk) += 1;
		c = c->next;
		//printf("found\n");
	}
}

void search(int16_t *x, int16_t d, int16_t *pk, uint8_t *y, uint8_t *sf) {
	int16_t ii, jj, mm, zz;
	if (d==0) { addxyz(x[0], x[1], (DIM == 3) ? x[2] : XMIN[2], pk, y, sf); return; }
	for (mm=-1; mm<=1; mm+=2) {
		addxyz(x[0], x[1], (DIM == 3) ? x[2]+d*mm : XMIN[2], pk, y, sf);
		for (jj=1; jj<=d-(mm+1)/2; jj+=1) {
			zz = (DIM == 3) ? x[2] + (d-jj) * mm : XMIN[2];
//			for (ii=0; ii<=jj; ii+=1) addxyz(x[0]+ii, x[1]-jj+ii, zz, pk, y, sf); //modified
			for (ii=0; ii<=jj; ii+=1) addxyz(x[0]+ii, x[1]+jj-ii, zz, pk, y, sf);
			for (ii=1; ii< jj; ii+=1) addxyz(x[0]+ii, x[1]-jj+ii, zz, pk, y, sf);
			for (ii=0; ii<=jj; ii+=1) addxyz(x[0]-ii, x[1]-jj+ii, zz, pk, y, sf);
			for (ii=1; ii< jj; ii+=1) addxyz(x[0]-ii, x[1]+jj-ii, zz, pk, y, sf);
		}
	}
}

void addgps(int16_t gps0, int16_t gps1, int16_t *pk, uint8_t *y) {
	char strgps[20];
	FILE *pfgps;
	sprintf(strgps, "gps_db/%d,%d", gps0, gps1);
	pfgps = fopen(strgps, "r");
	if (pfgps == NULL) return;
	while (1 == fscanf(pfgps, "%hhX", y+(*pk))) (*pk) += 1;
	fclose(pfgps);
}

SFSelect search_gps(DataCell *c, int16_t d, uint8_t sfsel) {
	static uint8_t y[25000];
	static int16_t rndidx[50];
	static int16_t lastgps[2];
	int16_t ii, jj, kk, dd, scs;
	int16_t *gps = c->gps;
	uint8_t sfmap = c->y;
	SFSelect sfs;
	kk = 0;
AGAIN:
	addgps(gps[0], gps[1], &kk, y); //dd==0
	for (dd = 1; dd <= d; dd += 1) {
		for (ii=0; ii<=dd; ii+=1) addgps(gps[0]+ii, gps[1]+dd-ii, &kk, y);
		for (ii=1; ii< dd; ii+=1) addgps(gps[0]+ii, gps[1]-dd+ii, &kk, y);
		for (ii=0; ii<=dd; ii+=1) addgps(gps[0]-ii, gps[1]-dd+ii, &kk, y);
		for (ii=1; ii< dd; ii+=1) addgps(gps[0]-ii, gps[1]+dd-ii, &kk, y);
	}
	if (kk == 0) {
		//d *= 2;
		printf("(%d,%d) => (%d,%d)\n", gps[0], gps[1], lastgps[0], lastgps[1]);
		gps = lastgps;
		goto AGAIN;
	}
	lastgps[0] = gps[0]; lastgps[1] = gps[1]; //printf("OK %d %d\n", lastgps[0], lastgps[1]);
	for (ii=0; ii < SFPKT[0]-1; ii+=1) rndidx[ii] = rand() % kk;
	sfs.optscs = -1;
	for (jj=0; jj<6; jj+=1) if ((sfmap >> jj) & 1) {
		scs = 0;
		for (ii=0; ii < SFPKT[jj]-1; ii+=1) scs += (y[rndidx[ii]] >> jj) & 1;
		if (scs > sfs.optscs) sfs.optscs = scs; //optsf=jj
	}
	sfs.selscs = 0;
	for (ii=0; ii < SFPKT[sfsel]-1; ii+=1) sfs.selscs += (y[rndidx[ii]] >> sfsel) & 1;
	return sfs;
}

void test(DataCell* c, int32_t testsize, int16_t adjlen, float relia, uint16_t mode) {
	static uint8_t y[25000], sf[25000];
	uint8_t ypred, ydiff, sfsel = 0, sfscs;
	int16_t kk, mm, positives, votes, sfpkt, maxpkt, steps, snr[ADRWIN], estsnr;
	int32_t ii, jj, pktopt = 0, pktsel = 0, attsel = 0, scssel = 0;
	int32_t tp[6]={0}, tn[6]={0}, pp[6]={0}, nn[6]={0}, scs[6]={0}, sel[6]={0};
	float prr, prrsf, thruput, timecost = 0.0;
	int16_t *px;
	DataCell* cc = c;
	SFSelect sfs = {0,0};
	srand(time(NULL));
	for (ii=0; ii<testsize-1; ii+=1) {
		if (mode & STAT) { sfsel = SFSTAT; goto SELECTED; }
		if (mode & PROB) {
			if ((ii+1) % adjlen) goto SELECTED;
			if (prr < relia && sfsel < 5) sfsel += 1;
			if (prr > relia + 0.05 && sfsel > 0) sfsel -= 1;
			goto SELECTED;
		}
		if (mode & ADR) {
			snr[ii % ADRWIN] = cc->x[0] - cc->x[1];
			if (cc != c + ii) { if (sfsel < 5) sfsel += 1; goto SELECTED; }
			if (ii < ADRWIN) goto SELECTED;
			if (ADRMAX) {
				estsnr = -20;
				for (jj=0; jj < ADRWIN; jj+=1) if (snr[jj] > estsnr) estsnr = snr[jj];
			} else {
				estsnr = 0;
				for (jj=0; jj < ADRWIN; jj+=1) estsnr += snr[jj];
				estsnr /= ADRWIN;
			}
			steps = (estsnr - (-20) - 10) / 3;
			if (steps < 0) sfsel = 5;
			else if (steps > 5) sfsel = 0;
			else sfsel = 5 - steps;
			goto SELECTED;
		}
		if (cc != c + ii) { sfsel = 5; goto SELECTED; }
		//printf("[%hd, %hd] => [", c[ii].x[0], c[ii].x[1]);
		//for (jj=0; jj<6; jj+=1) printf("%hd, ", (c[ii].y >> jj) & 1);
		kk = 0; ypred = 0; sfsel = UNDEF;
		px = (mode & GPSIN) ? cc->gps : cc->x;
		//printf("before search\n");
		for (jj=0; kk<K; jj+=1) search(px, jj, &kk, y, sf);
		//printf("after search\n");
		for (jj=0; jj<6; jj+=1) {
			positives = votes = 0;
			for (mm=0; mm<kk; mm+=1) {
				if (sf[mm]!=ALL && sf[mm]!=jj) continue;
				positives += (y[mm] >> jj) & 1;
				votes += 1;
			}
			//printf("kk=%d, votes=%d\n",kk,votes);
			ypred |= (positives > (float)votes * VOTMAJ[jj]) << jj;
			if (sfsel==UNDEF && ypred) sfsel = jj;
		}
		//printf("ii=%d\n",ii);
SELECTED:
		if (mode & SFREC) fprintf(psf, "%d\n", sfsel+7);
		sfscs = (c[ii].y >> sfsel) & 1;
		if (mode & NOGPS) {
			for (jj=0; jj<6; jj+=1) if ((c[ii].y >> jj) & 1) break;
			if (jj < 6) pktopt += SFPKT[jj];
			if (sfscs) pktsel += SFPKT[sfsel];
		}
		else {
			sfs = search_gps(c + ii, 5, sfsel);
			pktopt += sfs.optscs + 1;
			pktsel += sfs.selscs + sfscs;
		}
		attsel += 1; //SFPKT[sfsel]
		scssel += sfscs;
		sel[sfsel] += 1; //SFPKT[sfsel]
		scs[sfsel] += sfscs; //+ sfs.selscs
		if (mode & ORGTHPT) timecost += SFTIME[sfsel];
		if (sfscs || ((mode & KNN) && sfs.selscs > 0)) cc = c + ii + 1;
		if ((mode & HIST) && sfscs) { //|| sfs.selscs > 0
			c[ii].sf = sfsel;
			addlink(c + ii, mode);
		}
		ydiff = ypred ^ c[ii].y;
		for (jj=0; jj<6; jj+=1) {
			if ((c[ii].y >> jj) & 1) pp[jj]+=1; else nn[jj]+=1;
			if ((ydiff >> jj) & 1) {
				if (mode & ERREP) 
					fprintf(pef, "%d,%d,%d,%d\n", jj, c[ii].x[0], c[ii].x[1], c[ii].x[2]);
				continue;
			}
			if ((ypred >> jj) & 1) tp[jj]+=1; else tn[jj]+=1;
		}
		if ((ii+1) % adjlen) continue;
		prr = (float)scssel / attsel;
		for (jj=0; jj<5; jj+=1) {
			if (sel[jj] < 20) continue;
			prrsf = (float)scs[jj] / sel[jj];
			if (prr < relia && prrsf < relia && VOTMAJ[jj] < 0.89) VOTMAJ[jj] += 0.1;
			if (prr > relia + 0.05 && prrsf > relia && VOTMAJ[jj] > 0.11) VOTMAJ[jj] -= 0.05;
			scs[jj] = sel[jj] = 0;
			if (mode & TRAIN) printf("prr = %f, VOTMAJ[%d] = %f\n", prrsf, jj, VOTMAJ[jj]);
		}
		thruput = (float)pktsel / ((mode & RAWTHPT) ? adjlen : 
								((mode & ORGTHPT) ? timecost : pktopt));
		if ((mode & TEST) && (pktopt > 1000 || adjlen < 250 || (mode & ORGTHPT)))
			fprintf(pof, "%0.3f,%0.3f\n", prr, thruput); //, pktopt, sfsel,c[ii].gps[0], c[ii].gps[1]
		timecost = pktsel = pktopt = attsel = scssel = 0;
	}
	//if (mode & TEST) return;
	for (ii=0; ii<5; ii+=1) {
		printf("-- SF = %d --\n", ii+7);
		printf("TPR = %d/%d = %0.3f\n", tp[ii], pp[ii], (float)tp[ii] / pp[ii]);
		printf("TNR = %d/%d = %0.3f\n", tn[ii], nn[ii], (float)tn[ii] / nn[ii]);
		printf("ACC = %d/%d = %0.3f\n", tp[ii] + tn[ii], pp[ii] + nn[ii],
			(float)(tp[ii] + tn[ii]) / (pp[ii] + nn[ii]));
	}
}

int main(int argn, char* argv[]) {
	int32_t trainsize, testsize, ii;
	if (argn != 4) { printf("usage: %s <init.csv> <test.csv> <result.csv>\n", argv[0]); return 0; }
	trainsize = readfile(argv[1], traincell, 250, TRAIN);
	testsize = readfile(argv[2], testcell, -1, TEST);
	printf("trainsize=%d, testsize=%d\n", trainsize, testsize);
	for (ii=0; ii<10; ii+=1) test(traincell, trainsize, trainsize-1, 0.8, TRAIN | ORGTHPT /*| NOGPS*/);
	pof = fopen(argv[3], "w");
	pef = fopen("error.csv", "w");
	psf = fopen("sfsel.csv", "w");
	printf("---TESTING---\n");
	test(testcell, testsize, 300, 0.8, TEST | KNN | SFREC);
	//test(testcell, testsize, 300, 0.8, TEST | KNN | NOGPS | SFREC);
	//test(testcell, testsize, 300, 0.8, TEST | ADR | NOGPS | SFREC);
	//test(testcell, testsize, 300, 0.8, TEST | PROB | NOGPS | SFREC);
	//test(testcell, testsize, 300, 0.8, TEST | STAT | NOGPS | SFREC);
/*
	system("> ts.csv");
	for (ii=0; ii<testsize-100; ii+=100) {
		system("echo $(($(date +%s%N)/1000)) >> ts.csv");
		test(testcell+ii, 100, 99, 0.9, TEST | KNN | ORGTHPT);
	}
	fclose(pof);
	fclose(pef);
*/
/*
	for (ii=0; ii<trainsize; ii+=1) {
		printf("[%hd, %hd] => [", traincell[ii].x[0], traincell[ii].x[1]);
		for (jj=0; jj<6; jj+=1) printf("%hd, ", (traincell[ii].y >> jj) & 1);
		printf("\b\b]\n");
	}
	for (ii=0; ii<XMAX-XMIN+1; ii+=1) for (jj=0; jj<XMAX-XMIN+1; jj+=1) {
		if (jj==0) printf("\n");
		printf("%02x ", table[jj][ii]);
	}
*/
	return 0;
}





























