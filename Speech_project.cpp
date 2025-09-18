// Speech_project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <Windows.h>
#include <mmsystem.h>
#include <direct.h>  // For creating directories
#pragma comment(lib, "winmm.lib")
using namespace std;

//short int waveIn[LENGTH_WAV];
#define ORDER 12 
#define PI 3.142857142857
#define FRAME_SIZE 320 
#define T_MAX 200
#define double long double

int len[7][101][2];
int lentest[7][26][2];
double universe[60000][12];
double codebook_arr[60000];
double final_codebook[32][12];
double epsilon = 0.03;
double codebookU[32][12];

int N = 5;				//No. of States
int M = 32;				//No. of symbols in a state
int T;				    //length of Observation
double A[6][6];			//Contains transition probability
double B[6][33];		//Contains probability of symbol at respective state
double PIE[6];           //Initial probabitity of states
int O[200]; 
double alpha[T_MAX][6];	//Probability of Forward pass
double Beta[T_MAX][6];     //Probability of Backward pass
double gamma[T_MAX][6];
double Xita[T_MAX][6][6];
double preA[6][6];			//Contains previous transition probability
double preB[6][33];		//Contains previous probability of symbol at respective state
double prePI[6];			//previous initial probabitity of states
const char *Voice[6] = {"M", "F", "MPS", "FPS","MSRS","FSRS"};

double Del[T_MAX][6];		//Contains Probability of optimal state
int Psi[T_MAX][6];		//Contains path
int q[T_MAX];	

int readUniverse(char* filename){
	int length=0;
    FILE* file = fopen(filename, "r");
    if (file==NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
	int count=0;
    while (fscanf(file, "%lf",&universe[length][count++]) != EOF) {
		if(count == 12){
			length++;
			count=0;
		}   
    }
	//printf("%d",*length);

    fclose(file);
	return length;
}

int Signal(char* filename, double* signal) {
	int length=0;
    FILE* file = fopen(filename, "r");
    if (file==NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    while (fscanf(file, "%lf",&signal[length]) != EOF) {
        length++;
    }
	//printf("%d",*length);

    fclose(file);
	return length;
}

void writeAI(char* filename, double* S, int P,int frameNo) {
    FILE* file = fopen(filename, "a");
    if (!file) {
        printf("%s not Found\n", filename);
        exit(1);
    }
	fprintf(file, "Ai's of Frame%d:\n",frameNo);
    for (int i = 1; i <= P; i++) {
        fprintf(file, "%lf \t", S[i]);
    }
    fprintf(file, "\n\n");

    fclose(file);
}

void writeFile(char* filename, double* S, int P,int frameNo) {
    FILE* file = fopen(filename, "a");
    if (!file) {
        printf("%s not Found\n", filename);
        exit(1);
    }
    for (int i = 1; i <= P; i++) {
        fprintf(file, "%lf \t", S[i]);
    }
    fprintf(file, "\n\n");

    fclose(file);
}

void dcShift(double* signal, int length) {
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += signal[i];
    }
    double dcOffset = sum / length;//DC Shift of wave.
    for (int i = 0; i < length; i++) {  
        signal[i] -= dcOffset ;
    }
}

void Normalization(double* signal, int length) {
    double maxAmplitude = abs(signal[0]);
    for (int i = 1; i < length; i++) {
        if (abs(signal[i]) > maxAmplitude) {
            maxAmplitude = abs(signal[i]);
        }
    }
    double normalizationFactor = 5000 / maxAmplitude;
    for (int i = 0; i < length; i++) {
        signal[i] *= normalizationFactor;
    }
}

void steadyFrames(double* signal, int length, double frames[][FRAME_SIZE], int* numFrames) {
    double* energy = (double*)malloc((length / FRAME_SIZE) * sizeof(double));
    int frameIndex = 0;
	
    // Calculate energy for each frame
    for (int i = 0; i < length; i += FRAME_SIZE) {
        energy[frameIndex] = 0.0;
        for (int j = 0; j < FRAME_SIZE; j++) {
            energy[frameIndex] += signal[i + j] * signal[i + j];
        }
        frameIndex++;
    }

    double threshold = 0.8*energy[0];
    *numFrames = 0;

    for (int i = 0; i < frameIndex; i++) {
        if (energy[i] > threshold) {
            for (int j = 0; j < FRAME_SIZE; j++) {
                frames[*numFrames][j] = signal[i * FRAME_SIZE + j];
            }
            (*numFrames)++;
        }
    }
}

void calcRi(double* frame, int frameSize, double* RI, int P) {
    for (int i = 0; i <= P; i++) {
        RI[i] = 0.0;
        for (int j = 0; j < frameSize - i; j++) {
            RI[i] += frame[j] * frame[j + i];
        }
    }
}

void levinsonDurbin(double* RI, double* AI, int order) {
    double E[ORDER+1];
	double Alpha[ORDER+1][ORDER+1];
	double k[ORDER+1];

	E[0] = RI[0];
	for(int i=1;i<=order;i++) {
		double sum = 0.0;

		for(int j=1;j<=i-1;j++) {
			sum += Alpha[i-1][j] * RI[i-j];
		}

		k[i] = (RI[i] - sum) / E[i-1];

		Alpha[i][i] =  k[i];

		for(int j=1;j<=i-1;j++) {
			Alpha[i][j] = Alpha[i-1][j] - k[i] * Alpha[i-1][i-j];
		}

		E[i] = (1 - k[i] * k[i]) * E[i-1];
	}

	for(int i=1;i<=order;i++) {
		AI[i] = Alpha[order][i];
	}
}

// Function to compute Cepstral Coefficients from LPC coefficients
void ComputeCepstralCoefficients(double* CI,  double* Ai, int order) {
    double Wi[ORDER+1];

    for (int n = 1; n <= order; n++) {
        CI[n] = Ai[n];
        double sum = 0.0;
        for (int k = 1; k < n; k++) {
            if (n - k >= 0) {
                sum += k * Ai[n - k] * CI[k];
            }
        }
        CI[n] += sum / n;
        Wi[n] += CI[n] * 6 * sin(PI *( n / 12));
    }
}

int calculateCi(char*universeName,char* inputFile,char* outputFileC,int testNum){
			double signal[50000]; 
			int signalLength=Signal(inputFile, signal);

			// Perform DC shift
			dcShift(signal, signalLength);

			// Perform Normalization
			Normalization(signal, signalLength);

			double steadyFrame[200][FRAME_SIZE];
			int numFrames = 0;

			steadyFrames(signal, signalLength, steadyFrame, &numFrames);

			double* RI = (double*)malloc((ORDER + 1) * sizeof(double));
			double* AI = (double*)malloc((ORDER + 1) * sizeof(double));
			double* CI = (double*)malloc((ORDER + 1) * sizeof(double));
			if (RI == NULL || AI == NULL) {
				printf("Memory allocation failed for Ri's or Ai's array\n");
				exit(1);
			}
			FILE* file2 = fopen(outputFileC, "w");
			fclose(file2);

			for (int i = 0; i < numFrames; i++) {
				calcRi(steadyFrame[i], FRAME_SIZE, RI, ORDER);
				levinsonDurbin(RI, AI, ORDER);
				ComputeCepstralCoefficients(CI, AI, ORDER);
				// Save CI's to file
				writeFile(outputFileC, CI, ORDER,i+1);
				writeFile(universeName, CI, ORDER,i+1);
			}
			printf("Ci's are saved to %s\n", outputFileC);

			free(RI);
			free(AI);
			free(CI);
		return numFrames;
}

//Generating Ci for Test Data 
int calculateCiTest(char*universeName,char* inputFile,char* outputFileC,int testNum){
			double signal[50000]; 
			int signalLength=Signal(inputFile, signal);

			// Perform DC shift
			dcShift(signal, signalLength);

			// Perform Normalization
			Normalization(signal, signalLength);

			double steadyFrame[200][FRAME_SIZE];
			int numFrames = 0;

			steadyFrames(signal, signalLength, steadyFrame, &numFrames);

			double* RI = (double*)malloc((ORDER + 1) * sizeof(double));
			double* AI = (double*)malloc((ORDER + 1) * sizeof(double));
			double* CI = (double*)malloc((ORDER + 1) * sizeof(double));
			if (RI == NULL || AI == NULL) {
				printf("Memory allocation failed for Ri's or Ai's array\n");
				exit(1);
			}
			FILE* file2 = fopen(outputFileC, "w");
			fclose(file2);

			for (int i = 0; i < numFrames; i++) {
				calcRi(steadyFrame[i], FRAME_SIZE, RI, ORDER);
				levinsonDurbin(RI, AI, ORDER);
				ComputeCepstralCoefficients(CI, AI, ORDER);
				// Save CI's to file
				writeFile(outputFileC, CI, ORDER,i+1);
			}
			printf("Ci's are saved to %s\n", outputFileC);

			free(RI);
			free(AI);
			free(CI);
		return numFrames;
}

/*.......................LGB Starts..........................*/
// Generate a random index for initializing the codebook
int generateRandomIndex(int limit) {
    return rand() % (limit + 1);
}
// Calculate the squared Tokhura distance between a vector and the codebook
double* calculateDistance(int vectorIndex,  double** codebook, int k, int dimensions) {
    static double minDistanceResult[2]; 
    minDistanceResult[0] = 1e9; 

    double tokhuraWeights[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    double distance, diff;

    for (int i = 0; i < k; i++) {
        distance = 0;
        for (int j = 0; j < dimensions; j++) {
            diff = universe[vectorIndex][j] - codebook[i][j];
            distance += tokhuraWeights[j] * diff * diff;
        }
        if (distance < minDistanceResult[0]) {
            minDistanceResult[0] = distance;
            minDistanceResult[1] = i;
        }
    }

    return minDistanceResult;
}
// Perform K-means clustering to optimize the codebook
double** performKMeans(int datasetSize, double** codebook, int k, int dimensions) {
    double threshold = 0.000001;  
    double prevDistortion = 1e10, currDistortion = 0;

    // Allocate memory for sum of vectors and region sizes
    double** regionSums = (double**) malloc(k * sizeof(double*));
    double* regionCounts = (double*) calloc(k, sizeof(double));
    if (regionSums == NULL || regionCounts == NULL) {
        // handle allocation failure
        return NULL;
    }

    for (int i = 0; i < k; i++) {
        regionSums[i] = (double*) calloc(dimensions, sizeof(double));
        if (regionSums[i] == NULL) {
            // handle allocation failure
            return NULL;
        }
    }

    double* distanceResult;

    do {
        prevDistortion = currDistortion;
        currDistortion = 0;

        // Reset sums and counts
        for (int i = 0; i < k; i++) {
            regionCounts[i] = 0;
            for (int j = 0; j < dimensions; j++) {
                regionSums[i][j] = 0;
            }
        }

        // Assign each vector to the nearest codebook vector
        for (int i = 0; i < datasetSize; i++) {
            distanceResult = calculateDistance(i, codebook, k, dimensions);
            currDistortion += distanceResult[0];
            
            int closestRegion = (int)distanceResult[1];
            for (int j = 0; j < dimensions; j++) {
                regionSums[closestRegion][j] += universe[i][j];
            }
            regionCounts[closestRegion]++;
        }

        currDistortion /= datasetSize;

        // Update the codebook by averaging vectors in each region
        for (int i = 0; i < k; i++) {
            if (regionCounts[i] > 0) {
                for (int j = 0; j < dimensions; j++) {
                    codebook[i][j] = regionSums[i][j] / regionCounts[i];
                }
            }
        }

    } while (fabs(prevDistortion - currDistortion) > threshold);

    printf("Converged with distortion: %f\n", currDistortion);

    // Free allocated memory
    for (int i = 0; i < k; i++) {
        free(regionSums[i]);
    }
    free(regionSums);
    free(regionCounts);

    return codebook;
}
// Convert universe to array (LBG helper function)
void universeToArray( int size, int p) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < p; j++) {
            codebook_arr[j] += universe[i][j];
        }
    }

    for (int i = 0; i < p; i++) {
        codebook_arr[i] /= size;
    }

}
// LBG algorithm to generate the codebook
double** LBG( int size, int k, int p, double epsilon) {
    double delta = 0.001;
    int codebook_size = 1;

    double** temp_codebook = NULL;
    double** final_codebook = (double**) malloc(sizeof(double*));

    final_codebook[0] = (double*) calloc(p, sizeof(double));
    for (int j = 0; j < p; j++) {
        final_codebook[0][j] = codebook_arr[j];
    }

    while (codebook_size <= k) {
        temp_codebook = (double**) malloc(2 * codebook_size * sizeof(double*));

        for (int i = 0; i < 2 * codebook_size; i++) {
            temp_codebook[i] = (double*) calloc(p, sizeof(double));
        }

        for (int i = 0; i < codebook_size; i++) {
            for (int j = 0; j < p; j++) {
                temp_codebook[2 * i][j] = final_codebook[i][j] * (1 + epsilon);
                temp_codebook[2 * i + 1][j] = final_codebook[i][j] * (1 - epsilon);
            }
        }

		printf("\n%d\n",codebook_size);
        codebook_size *= 2;
        free(final_codebook);
        final_codebook = performKMeans(size, temp_codebook, codebook_size, p);
    }

    return final_codebook;
}

/*............................ HMM strats.............................*/
void openA(char *filename)
{
	FILE *fp = NULL;
	int err = fopen_s(&fp,filename,"r");		//Open File
	if(err != NULL)
	{
		printf("\n Cannot open");
		system("pause");
		exit(1);
	}
	int i = 1, j = 1;
	//Here data is being retrieved in the array.
	while(!feof(fp))
	{
		double num;
		if(fscanf(fp,"%lf",&num) == 1)
		{
			A[i][j] = num;
			j++;
			if(j == 6)
			{
				j = 1;
				i++;
			}
		}
	}
	fclose(fp);	
}
void openB(char *filename)
{
	FILE *fp = NULL;
	int err = fopen_s(&fp,filename,"r");		//Open file
	if(err != NULL)
	{
		printf("\n Cannot open");
		system("pause");
		exit(1);
	}
	int i = 1, j = 1;
	//Here data is being retrieved in the array.
	while(!feof(fp))
	{
		double num;
		if(fscanf(fp,"%lf",&num) == 1)
		{
			B[i][j] = num;
			j++;
			if(j == 33)
			{
				j = 1;
				i++;
			}
		}
	}
	fclose(fp);	
}
void openPI(char *filename)
{
	FILE *fp = NULL;
	int err = fopen_s(&fp,filename,"r");		//Open file
	if(err != NULL)
	{
		printf("\n Cannot open");
		system("pause");
		exit(1);
	}
	int j = 1;
	//Here data is being retrieved in the array.
	while(!feof(fp))
	{
		double num;
		if(fscanf(fp,"%lf",&num) == 1)
		{
			PIE[j] = num;
			j++;
		}
	}
	fclose(fp);	
}
void openObs(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("\n Cannot open %s", filename);
        exit(1);
    }

    // Calculate T and allocate O dynamically
    T = 1;
    int temp;
    while (fscanf(fp, "%d", &O[T]) == 1) {
        T++;
	}
	T--;
    fclose(fp);
}
double forwardPropogation()
{
	//Step 1 : Initialization
	for(int i=1;i<=N;i++)
		alpha[1][i] = PIE[i]*B[i][O[1]];
	
	//Step 2 : Induction
	for(int t = 2;t<=T;t++)
		for(int j=1;j<=N;j++)
		{
			alpha[t][j] = 0;
			for(int i=1;i<=N;i++)
				alpha[t][j] += alpha[t-1][i] * A[i][j];
			alpha[t][j] *= B[j][O[t]];
		}
	double ans = 0;

	//Termination
	for(int i=1;i<=N;i++)
		ans += alpha[T][i];

	//Printing of array alpha
	/*
	printf("Forward Propogaton:-\n");
	for(int i=1;i<=T;i++)
	{
		for(int j=1;j<=N;j++)
			printf("%e ",alpha[i][j]);
		printf("\n");
	}
	*/
	return ans;
}
void backwardPropogation()
{
	for(int i=1;i<=N;i++)
		Beta[T][i] = 1;
	for(int t = T-1;t>=1;t--)
		for(int i=1;i<=N;i++)
		{
			Beta[t][i] = 0;
			for(int j=1;j<=N;j++)
				Beta[t][i] += A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
		}
	/*
	printf("Backward Propogaton:-\n");
	for(int i=1;i<=T;i++)
	{
		for(int j=1;j<=N;j++)
			printf("%e ",Beta[i][j]);
		printf("\n");
	}
	*/
}
void calculateGamma(int t)
{
	//Calculation of Gamma
	for(int i=1;i<=N;i++)
	{
		gamma[t][i] = 0;
		for(int j=1;j<=N;j++)
			gamma[t][i] += Xita[t][i][j];
	}
}
void EitaGammaCalculation()
{
	//Calculation of Eita
	for(int t = 1;t<T;t++)
	{
		double deno = 0;
		for(int i=1;i<=N;i++)
			for(int j=1;j<=N;j++)
				deno += alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
		for(int i=1;i<=N;i++)
			for(int j=1;j<=N;j++)
				Xita[t][i][j] = (alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j])/deno;
		calculateGamma(t);			//Calculation of gamma at time t.
	}
	/*
	printf("\n\n\nGamma:-\n");
	for(int t=1;t<T;t++)
	{
		for(int i=1;i<=N;i++)
			printf("%e ",gamma[t][i]);
		printf("\n");
	}
	*/
}
void newModelCalculation()
{
		for(int i=1;i<N;i++)
			prePI[i] = PIE[i];
		for(int i=1;i<=N;i++)
			for(int j=1;j<=N;j++)
				preA[i][j] = A[i][j];
		for(int i=1;i<=N;i++)
			for(int j=1;j<=M;j++)
				preB[i][j] = B[i][j];
	//Calculation of PI
	for(int i=1;i<=N;i++)
		PIE[i] = gamma[1][i];

	//Calculation of A
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++)
		{
			double num = 0, deno = 0;
			for(int t = 1;t<T;t++)
			{
				num += Xita[t][i][j];
				deno += gamma[t][i];
			}
			A[i][j] = num/deno;
		}

		//Calculation of B
		for(int j=1;j<=N;j++){
			for(int k=1;k<=M;k++)
			{
				double num = 0, deno = 0;
				for(int t=1;t<T;t++)
				{
					if(O[t] == k)
						num += gamma[t][j];
					deno += gamma[t][j];
				}
				if(deno!=0)
					B[j][k] = num/deno;
				if(B[j][k]<1e-170)
					B[j][k]=1e-15;
			}
		}
		double sum=0;
		for(int k = 1; k <= M; k++){
			sum+=B[i][k];
		}
		for(int k = 1; k <= M; k++){
			B[i][k]/=sum;
		}
		/*
		printf("\n\n\nnew PI:-\n");
		for(int i=1;i<N;i++)
			printf("%lf ",PIE[i]);
		printf("\n\n\nnew A:-\n");
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
				printf("%lf ",A[i][j]);
			printf("\n");
		}
		printf("\n\n\nnew B :-\n");
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=M;j++)
				printf("%lf ",B[i][j]);
			printf("\n");
		}
		*/
	}
			
}
double EM(){

	double prob = forwardPropogation();		//This function implement forward propogation and calculates alpha matrix.
	backwardPropogation();					//This function implement backward propogation and calculates beta matrix.
	EitaGammaCalculation();					//This calculates Eita and Gamma
	newModelCalculation();					//This calculates new mode A,B and PI
	return prob;
}
void viterbiAlgo()
{
	//Step 1 : Initialization
	for(int i = 1;i<=N;i++)
	{
		Del[1][i] = PIE[i] * B[i][O[1]];
		Psi[1][i] = 0;
	}

	//Step 2 : Recursion
	for(int t = 2;t<=T;t++)
		for(int j = 1;j<=N;j++)
		{
			double mx = -1;
			int index=0;
			for(int i = 1;i<=N;i++)
				if((Del[t-1][i]*A[i][j])>mx)
				{
					mx = Del[t-1][i]*A[i][j];
					index = i;
				}
			Del[t][j] = mx*B[j][O[t]];
			Psi[t][j] = index;
		}

	//Step 3 : Termination
	double pStar = -1;
	for(int i=1;i<=N;i++)
	{
		if(Del[T][i]>pStar)
		{
			pStar = Del[T][i];
			q[T] = i;
		}
	}

	//Step 4 : State Sequence
	for(int t = T-1;t>=1;t--)
	{
		q[t] = Psi[t+1][q[t+1]];
	}

	//Printing optimal state sequence and probability of it.
	for(int t = 1;t<T;t++)
		printf("%d->",q[t]);	
	printf("%d\n",q[T]);
	printf("Probability of the most optimal state sequence : %e\n",pStar);
}

/*...................Observation file generation.....................*/
int findObs( double currCi[12],  int k, int dimensions) {
    double minDistance =INT_MAX;
    int minDistanceindex = -1 ; 

    double tokhuraWeights[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    double distance, diff;

    for (int i = 0; i < k; i++) {
        distance = 0;
        for (int j = 0; j < 12; j++) {
            diff = currCi[j] - codebookU[i][j];
            distance += tokhuraWeights[j] * diff * diff;
        }
        if (distance < minDistance) {
            minDistance = distance;
            minDistanceindex = i;
        }
    }

    return minDistanceindex;
}
void findObervationFile(char* inputNameC, char* obserFile,int count){
	double currCi[12];
	int obser[200];
	int currN=0;
			FILE *file3=fopen(inputNameC,"r");
			if(file3==NULL){
				printf("file can't be opened");
				exit(1);
			}
			currN=0;
			
			for(int k=0;k<count;k++){
				for(int l=0; l<12; l++){
					fscanf(file3,"%lf",&currCi[l]);
				}
				int obs = findObs(currCi,32,ORDER);
				obser[k]=obs+1;
			}
			FILE *fileO = fopen(obserFile,"w");
			for(int k=0;k<count;k++){
				fprintf(fileO,"%d\t", obser[k]);
			}
			fclose(fileO);
			fclose(file3);
}

/*..................Traning the data................................*/
void trainingData(int iteration) {
	char fileAavg[50], fileBavg[50], filePiAvg[50];
	for(int i=0;i<6;i++){
		if(iteration != 0){
			sprintf(fileAavg, "avg/model_%s/A_%s_%d_avg.txt",Voice[i],Voice[i],iteration-1);
			sprintf(fileBavg, "avg/model_%s/B_%s_%d_avg.txt",Voice[i],Voice[i],iteration-1);
			sprintf(filePiAvg, "avg/model_%s/pi_%s_%d_avg.txt",Voice[i],Voice[i],iteration-1);
		}else{
			sprintf(fileAavg, "initialModel/A.txt");
			sprintf(fileBavg, "initialModel/B.txt");
			sprintf(filePiAvg, "initialModel/pi.txt");
		}

		// Training 
		for(int j=1;j<=100;j++){
			openA(fileAavg);
			openB(fileBavg);
			openPI(filePiAvg);
			char inputNameO[50];
			sprintf(inputNameO, "Observation/data_%s_%d_obs.txt",Voice[i], j);
			openObs(inputNameO);
			double curr = EM(),prev;
			printf("Initial Model Optimal probability = %e\n",curr);
			int iteration=0;
			do
			{
				iteration++;
				prev = curr;
				curr = EM();
				printf("Model probability = %e\n",curr);
			}while(curr>prev && iteration<100);

			//viterbiAlgo();			
			
			//Printing new PI and new A.
			char Anew[50];
			sprintf(Anew, "model/model_%s/A_%d.txt",Voice[i], j);
			char Bnew[50];
			sprintf(Bnew, "model/model_%s/B_%d.txt",Voice[i], j);
			char piNew[50];
			sprintf(piNew, "model/model_%s/pi_%d.txt",Voice[i], j);

			FILE *filePi = fopen(piNew,"w");
			printf("\n\n\nnew PI:-\n");
			for(int i=1;i<N;i++)
				fprintf(filePi,"%e ",prePI[i]);
			fclose(filePi);

			FILE *fileA = fopen(Anew,"w");
			printf("\n\n\nnew A:-\n");
			for(int i=1;i<=N;i++)
			{
				for(int j=1;j<=N;j++)
					fprintf(fileA,"%e ",preA[i][j]);
				fprintf(fileA,"\n");
			}
			fclose(fileA);

			FILE *fileB = fopen(Bnew,"w");
			printf("\n\n\nnew B:-\n");
			for(int i=1;i<=N;i++)
			{
				for(int j=1;j<=M;j++)
					fprintf(fileB,"%e ",preB[i][j]);
				fprintf(fileB,"\n");
			}
			fclose(fileB);
		}
	}
}
void averageModel(int iteration) {
    for (int i = 0; i < 6; i++) {  // Loop through models
		// Arrays to hold cumulative sums for averaging
		double avgA[6][6] = {0}; 
		double avgB[6][33] = {0}; 
		double avgPi[6] = {0};
        for (int j = 1; j <= 100; j++) {  // Loop through each file in the model
            // Construct file paths
            char Aavg[50], Bavg[50], piAvg[50];
            sprintf(Aavg, "model/model_%s/A_%d.txt", Voice[i], j);
            sprintf(Bavg, "model/model_%s/B_%d.txt", Voice[i], j);
            sprintf(piAvg, "model/model_%s/pi_%d.txt",Voice[i], j);

            // Open files
            FILE *fileA = fopen(Aavg, "r");
            FILE *fileB = fopen(Bavg, "r");
            FILE *filePi = fopen(piAvg, "r");
            if (!fileA || !fileB || !filePi) {
                printf("Failed to open one of the files for model_%s, sequence %d.\n",Voice[i], j);
                continue;
            }

            // Accumulate values for avgA
            double num = 0;
            for (int m = 1; m <= N; m++) {
                for (int n = 1; n <= N; n++) {
                    if (fscanf(fileA, "%lf", &num) == 1) {
                        avgA[m][n] += num;
                    }
                }
            }

            // Accumulate values for avgB
            for (int m = 1; m <= N; m++) {
                for (int n = 1; n <= M; n++) {
                    if (fscanf(fileB, "%lf", &num) == 1) {
                        avgB[m][n] += num;
                    }
                }
            }

            // Accumulate values for avgPi
            for (int m = 1; m <= N; m++) {
                if (fscanf(filePi, "%lf", &num) == 1) {
                    avgPi[m] += num;
                }
            }

            // Close files
            fclose(fileA);
            fclose(fileB);
            fclose(filePi);
        }

		// Divide by the total number of models (100) to get the average
		int totalModels = 100;
		for (int m = 1; m <= N; m++) {
			for (int n = 1; n <= N; n++) {
				avgA[m][n] /= totalModels;
			}
			for (int n = 1; n <= M; n++) {
				avgB[m][n] /= totalModels;
			}
			avgPi[m] /= totalModels;
		}


		// Print or return averages as needed
		char fileAavg[50], fileBavg[50], filePiAvg[50];
        sprintf(fileAavg, "avg/model_%s/A_%s_%d_avg.txt",Voice[i],Voice[i],iteration);
        sprintf(fileBavg, "avg/model_%s/B_%s_%d_avg.txt",Voice[i],Voice[i],iteration);
        sprintf(filePiAvg, "avg/model_%s/Pi_%s_%d_avg.txt",Voice[i],Voice[i],iteration);
		FILE *fileA = fopen(fileAavg, "w");
		FILE *fileB = fopen(fileBavg, "w");
		FILE *filePi = fopen(filePiAvg, "w");
		printf("Averaged Model A:\n");
		for (int m = 1; m <= N; m++) {
			for (int n = 1; n <= N; n++) {
				fprintf(fileA,"%e\t ", avgA[m][n]);
			}
			fprintf(fileA,"\n");
		}
		printf("\nAveraged Model B:\n");
		for (int m = 1; m <= N; m++) {
			for (int n = 1; n <= M; n++) {
				fprintf(fileB,"%e\t", avgB[m][n]);
			}
			fprintf(fileB,"\n");
		}
		printf("\nAveraged Initial State Probabilities (PI):\n");
		for (int m = 1; m <= N; m++) {
			fprintf(filePi,"%e\t ", avgPi[m]);
		}
		fprintf(filePi,"\n");

		fclose(fileA);
		fclose(fileB);
		fclose(filePi);
	}
}
int testData(char* testFile,int it){
	openObs(testFile);
	double maxprob=0;
	int digit=-1;
	double prob=0;
	for(int k=0;k<6;k++){
		char fileAavg[50], fileBavg[50], filePiAvg[50];
		sprintf(fileAavg, "avg/model_%s/A_%s_%d_avg.txt",Voice[k],Voice[k],it-1);
		sprintf(fileBavg, "avg/model_%s/B_%s_%d_avg.txt",Voice[k],Voice[k],it-1);
		sprintf(filePiAvg, "avg/model_%s/pi_%s_%d_avg.txt",Voice[k],Voice[k],it-1);
		openA(fileAavg);
		openB(fileBavg);
		openPI(filePiAvg);
		prob = forwardPropogation();
		if(prob > maxprob){
			maxprob = prob;
			digit = k;
		}
	}
			
	return digit;
}

int _tmain(int argc, _TCHAR* argv[])
{
	//Training 
	char universeName[50];
	int count;
	sprintf(universeName, "universe.txt");
	FILE* fileU = fopen(universeName, "w");
	fclose(fileU);
	for(int s=0;s<6;s++){
		int testNum=1;
		while(testNum<=100){
			char inputName[50];
			char outputNameC[50];
			sprintf(inputName, "train_input/data_%s_ (%d).txt",Voice[s], testNum);
			sprintf(outputNameC, "Ci/data_%s_%d_ci.txt",Voice[s], testNum);
		
			int numFrame = calculateCi(universeName,inputName,outputNameC,testNum);
			len[s+1][testNum][1]= numFrame;
			testNum++;
		}
	 }

	//Read Universe
	count = readUniverse(universeName);

	//codebook Generation
	double** codebook = LBG( count ,32, ORDER, epsilon);
	FILE *fileC = fopen("codebook.txt","w");
    printf("\nFinal codebook:\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < ORDER; j++) {
			codebookU[i][j] = codebook[i][j];
            fprintf(fileC,"%f\t", codebook[i][j]);
        }
        printf("\n");
    }
	fclose(fileC);

	// creating Obesrvation for each file
	char inputNameC[50];
	char obserFile[50];
	for(int s=0;s<6;s++){
		int testNum=1;
		while(testNum<=100){
			sprintf(inputNameC, "Ci/data_%s_%d_ci.txt",Voice[s], testNum);
			sprintf(obserFile, "Observation/data_%s_%d_obs.txt",Voice[s], testNum);
			findObervationFile(inputNameC,obserFile,len[s+1][testNum][1]);
			testNum++;
		}
	 }
	 
	// create model sequence

	//Training
	int it=3; 
	for (it =0;it<2;it++){
		trainingData(it);
		averageModel(it);
	}

	//Testing
	for(int s=0;s<6;s++){
		int testNum=1;
		while(testNum<=25){
			char inputName[50];
			char outputNameC[50];
			sprintf(inputName, "Test_input/data_%s_ (%d).txt",Voice[s], testNum);
			sprintf(outputNameC, "Ci_test/data_%s_%d_ci.txt",Voice[s], testNum);
		
			int numFrame = calculateCiTest(universeName,inputName,outputNameC,testNum);
			lentest[s+1][testNum][1]= numFrame;
			testNum++;
		}
	 }
	// creating Obesrvation for each file
	for(int s=0;s<6;s++){
		int testNum=1;
		while(testNum<=25){
			sprintf(inputNameC, "Ci_test/data_%s_%d_ci.txt",Voice[s], testNum);
			sprintf(obserFile, "Observation_test/data_%s_%d_obs.txt",Voice[s], testNum);
			findObervationFile(inputNameC,obserFile,lentest[s+1][testNum][1]);
			testNum++;
		}
	 }

	//Start Testing of dataset
	float accuracy=0;
	int Digit;
	for(int i=0;i<6;i++){
		char testFile[50];
		for(int j=1;j<=25;j++){
			sprintf(testFile, "Observation_test/data_%s_%d_obs.txt",Voice[i], j);
			Digit = testData(testFile,it);
			printf("\nFor given voice of %s Estimated voice of %s\n",Voice[i],Voice[Digit]);
			if(Voice[i]==Voice[Digit])
				accuracy++;
		}
	}
	accuracy = (accuracy/125)*100;
	printf("Accuracy of testing is %f%%",accuracy);
 	return 0;
}

