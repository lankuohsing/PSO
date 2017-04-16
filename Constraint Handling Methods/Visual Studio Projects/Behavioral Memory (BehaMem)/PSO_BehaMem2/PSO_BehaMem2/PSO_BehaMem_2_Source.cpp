#include <iostream>
#include "PSO_BehaMem_2_Header.h"
using namespace std;
int main(int argc, char* argv[])
{
	clock_t t1 = clock();	//The moment the program starts
	int test_number=0;		//number of tests
	int failed_number=0;	//number of failed tests
	ParticleSwarm p;		//definite a population object 
	int i,j,k;				//loop variables
	
	//array of which the data is to be analyzed in Excel
	double run_best[TESTTIMES][RUNTIMES];
	double final_best[TESTTIMES];
	double convergence[PSO_maxgen];
	bool is_convergent=true;

	FILE *PSOData;			//file pointer
	FILE *EXCELData;
	double  *gbestStore;	//to store "gbest"
	double Best_Value,Best_Position[GANVARS]={0.0};//To store the best value and best position
	
	//在堆上分配gbestStore，保存每一代的gbest。是数组。每运行一次，保存一次。 
	gbestStore = new double[PSO_maxgen];
	
	//initialize the array with 0.0
	for (i = 0; i < PSO_maxgen; i++)
	{
			gbestStore[i] = 0.0;
	}
	
	PSOData = fopen("PSO_BehaMem2.dat","w+");//clear the file and write new data in  
	if (PSOData == NULL)
	{
		printf("The file is not exit!\n");
		exit(1);
	}

	EXCELData = fopen("PSO_BehaMem2_EXCEL.dat","w+");//clear the file and write new data in  
	if (PSOData == NULL)
	{
		printf("The file is not exit!\n");
		exit(1);
	}

	else
	{
		for (test_number=0;test_number<TESTTIMES;test_number++)
		{
			Best_Value=numeric_limits<double>::max();  //maximal value of double data type
			fprintf(PSOData,"Test number:  %d\n",test_number);
			//the program runs RUNTIME times
			for (i = 0; i < RUNTIMES; i++)
			{
				p.init();//to initilaze	the population
				p.search(gbestStore);	
				run_best[test_number][i]=p.gbest;
				fprintf(PSOData,"Run number:  %d, Best values:  %lf  \n",i,p.gbest);
				fprintf(PSOData,"\nBest position: ");
				for (int j=0;j<GANVARS;j++)
				{
					fprintf(PSOData,"%lf ",p.gbest_pos[j]);
				}
				fprintf(PSOData,"\n");
				fprintf(PSOData,"-----------------------------------------------------------------------------------------\n");
			
				//保存每一代的gbest,可以屏蔽掉下面两行。使文件简洁易读。 
				for (int row = 0; row <PSO_maxgen; row++)				
					fprintf(PSOData,"%d   %lf\n",row,gbestStore[row]);	

				if (p.gbest<=Best_Value	
			//		&&Is_feasible(p.gbest_pos)
				)
				{
				   
				   if (is_convergent)
				   {
					   for (int con_gen=0;con_gen<PSO_maxgen;con_gen++)
					   {
						   convergence[con_gen]=gbestStore[con_gen];
					   }
					   is_convergent=false;
				   }
				   

				   Best_Value=p.gbest;
				   for (int k = 0; k < GANVARS; k++)
				   {
					   Best_Position[k]=p.gbest_pos[k];
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
				final_best[test_number]=Best_Value;

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
   		delete  []gbestStore;
		fclose(PSOData);
		for (int i=0;i<TESTTIMES;i++)
		{
			 fprintf(EXCELData,"%lf\n",final_best[i]);

		}
		fprintf(EXCELData,"---------\n");
		for (int i=0;i<TESTTIMES;i++)
		{
			for (int j=0;j<RUNTIMES;j++)
			{
				 fprintf(EXCELData,"%lf\n",run_best[i][j]);
			}

		}
		fprintf(EXCELData,"---------\n");
		for (int i=0;i<PSO_maxgen;i++)
		{
			fprintf(EXCELData,"%lf\n",convergence[i]);
		}
		fclose(EXCELData);

		
	}
//	cout<<"Failed times: "<<failed_number<<endl;
	clock_t t2 = clock();
	cout<<"It takes: "<<t2-t1<<"mS"<<endl;

	system("pause");
	return 0;
}