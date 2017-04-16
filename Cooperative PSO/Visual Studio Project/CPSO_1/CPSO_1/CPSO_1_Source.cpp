#include <iostream>
#include "CPSO_1_ Header.h"
using namespace std;
int main(int argc, char* argv[])
{
	clock_t t1 = clock();	//the moment the program starts
	int test_number=0;		//the times of tests
	int failed_number=0;	//the times of failed tests
	ParticleSwarm p;		//defining an objective of particle population
	int j,k;

	FILE *PSOData;			//a file pointer
	double  *resultStore;	//an array to store the gbest
	double Best_Value,Best_Position[GANVARS]={0.0};//To store the best value and best position

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
		for (test_number=0;test_number<TESTTIMES;test_number++)
		{
			Best_Value=numeric_limits<double>::max();
			fprintf(PSOData,"Test number:  %d\n",test_number);
			//the program will run RUNTIMES times
			for (int i = 0; i < RUNTIMES; i++)
			{
				p.init();//to initilaze
				p.search(resultStore);	

				fprintf(PSOData,"Run number:  %d, Best values:  %lf  \n",i,p.gbest);
				fprintf(PSOData,"\nBest position: ");
				for (int j=0;j<GANVARS;j++)
				{
					fprintf(PSOData,"%lf ",p.COOP_pos[j]);
				}
				fprintf(PSOData,"\n");
				fprintf(PSOData,"-----------------------------------------------------------------------------------------\n");

				//保存每一代的gbest,可以屏蔽掉下面两行。使文件简洁易读。 
				for (int row = 0; row <PSO_maxgen; row++)				
					fprintf(PSOData,"%d   %lf\n",row,resultStore[row]);	

				if (p.gbest<=Best_Value
					&&Is_feasible(COOP_pos_global)
					)
				{
					Best_Value=p.gbest;
					for (int k = 0; k < GANVARS; k++)
					{
						Best_Position[k]=COOP_pos_global[k];
					}
				}

			}
			fprintf(PSOData,"*********************************************************************************************\n");
			if (Best_Value==numeric_limits<double>::max())
			{
				failed_number++;
				fprintf(PSOData,"Failed!\n");
				printf("Failed\n");
			}
			else
			{
				fprintf(PSOData,"Final Best values:  %lf  \n",Best_Value);
				fprintf(PSOData,"Final Best position: ");
				for (int k = 0; k < GANVARS; k++)
				{
					fprintf(PSOData,"%lf ",Best_Position[k]);
				}
				fprintf(PSOData,"\n");
				printf("Final Best values:  %lf  \n",Best_Value);
				printf("Final Best position: ");
				for (int k = 0; k < GANVARS; k++)
				{
					printf("%lf ",Best_Position[k]);
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