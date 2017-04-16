#include <iostream>
#include <omp.h>
#include "CPSO_OPENMP_1_Header.h"
using namespace std;
int main(int argc, char* argv[])
{
	clock_t t1 = clock();	//the moment the program starts
	int test_number=0;		//the times of tests
	int failed_number=0;	//the times of failed tests
	ParticleSwarm p[TESTTIMES];		//defining an objective of particle population
	int j,k;

	FILE *PSOData;			//a file pointer
	double  *resultStore;	//an array to store the gbest
	double Best_Value[TESTTIMES],Best_Position[TESTTIMES][GANVARS]={0.0};//To store the best value and best position

	////////////////////////////allocate space for resultStore and initialize it/////////////////////////////
	resultStore = new double[PSO_maxgen];
/*	for (int i = 0; i < PSO_maxgen; i++)
	{
		resultStore[i] = 0.0;
	}  */
	////////////////////////////allocate space for resultStore and initialize it/////////////////////////////

	///////////////////////////////open the file//////////////////////////////////
	PSOData = fopen("CPSO_1.dat","w+");//clear the file and write new data into it。 
	if (PSOData == NULL)
	{
		printf("The file is not exit!\n");
		exit(1);
	}
	///////////////////////////////open the file//////////////////////////////////

	else
	{
		
		#pragma omp parallel for
		for (test_number=0;test_number<TESTTIMES;test_number++)
		{
			Best_Value[test_number]=numeric_limits<double>::max();
			fprintf(PSOData,"Test number:  %d\n",test_number);
			//the program will run RUNTIMES times
		//	#pragma omp parallel for
			for (int i = 0; i < RUNTIMES; i++)
			{
				p[test_number].init();//to initilaze
				p[test_number].search(resultStore);	

				fprintf(PSOData,"Run number:  %d, Best values:  %lf  \n",i,p[test_number].gbest);
				fprintf(PSOData,"\nBest position: ");
				for (int j=0;j<GANVARS;j++)
				{
					fprintf(PSOData,"%lf ",p[test_number].COOP_pos[j]);
				}
				fprintf(PSOData,"\n");
				fprintf(PSOData,"-----------------------------------------------------------------------------------------\n");

				//保存每一代的gbest,可以屏蔽掉下面两行。使文件简洁易读。 
				for (int row = 0; row <PSO_maxgen; row++)				
					fprintf(PSOData,"%d   %lf\n",row,resultStore[row]);	

				if (p[test_number].gbest<=Best_Value[test_number]
					&&Is_feasible(p[test_number].COOP_pos)
					)
				{
					Best_Value[test_number]=p[test_number].gbest;
					for (int k = 0; k < GANVARS; k++)
					{
						Best_Position[test_number][k]=p[test_number].COOP_pos[k];
					}
				}

			}
			fprintf(PSOData,"*********************************************************************************************\n");
			if (Best_Value[test_number]==numeric_limits<double>::max())
			{
				failed_number++;
				fprintf(PSOData,"Failed!\n");
				printf("Failed\n");
			}
			else
			{
				fprintf(PSOData,"Final Best values:  %lf  \n",Best_Value[test_number]);
				fprintf(PSOData,"Final Best position: ");
				for (int k = 0; k < GANVARS; k++)
				{
					fprintf(PSOData,"%lf ",Best_Position[test_number][k]);
				}
				fprintf(PSOData,"\n");
				printf("Final Best values:  %lf  \n",Best_Value[test_number]);
				printf("Final Best position: ");
				for (int k = 0; k < GANVARS; k++)
				{
					printf("%lf ",Best_Position[test_number][k]);
				}
				printf("\n");
			}

			fprintf(PSOData,"*********************************************************************************************\n");

		}

		//扫尾工作 
		delete  []resultStore;
		fclose(PSOData);


	}
	cout<<"Failed times: "<<failed_number<<endl;
	clock_t t2 = clock();
	cout<<"It takes: "<<t2-t1<<"mS"<<endl;

	system("pause");
	return 0;
}