#include <iostream>
#include "PSO_SF_4_Header.h"
using namespace std;

int main(int argc, char* argv[])
{
	clock_t t1 = clock();	//the time when the program starts
	int test_number=0;		//test number
	int failed_number=0;	//failed number	**mark
	ParticleSwarm p;		//an objective of particle swarm
	int i,j,k;				//loop variables

	FILE *PSOData;			//file pointer
	FILE *PSO_EXCEL;			//file pointer of the excel data

	double  *resultStore;	//to store the gbests
	double Best_Value,Best_Position[GANVARS]={0.0};//To store the best value and best position

	double *excel_final_best=new double[TESTTIMES];
	double **excel_run_best=new double *[TESTTIMES];
	for (i=0;i<TESTTIMES;i++)
	{
		excel_run_best[i]=new double [RUNTIMES];
	}
	
	double *excel_con_series=new double[PSO_maxgen];
	bool excel_series=true;
	//在堆上分配resultStore，保存每一代的gbest。是数组。每运行一次，保存一次。 
	resultStore = new double[PSO_maxgen]; //There is going to be PSO_maxgen iterations 

	//to initialize the array with 0.0
	for (i = 0; i < PSO_maxgen; i++)
	{
		resultStore[i] = 0.0; 
	}

	PSOData = fopen("PSO_SF_4.dat","w+");//to clear the file and write into new data。 
	if (PSOData == NULL)
	{
		printf("The file is not exit!\n");
		exit(1);
	}

	PSO_EXCEL = fopen("PSO_SF_4_EXCEL.dat","w+");//to clear the file and write into new data。 
	if (PSO_EXCEL == NULL)
	{
		printf("The file PSO_SF_1_EXCEL is not exit!\n");
		exit(1);
	}

	else
	{
		for (test_number=0;test_number<TESTTIMES;test_number++)
		{
			Best_Value=numeric_limits<double>::max();  //to initialize the Best_Value
			fprintf(PSOData,"Test number:  %d\n",test_number);

			//in every test, the programruns RUNTIMES times
			for (i = 0; i < RUNTIMES; i++)
			{
				do 
				{
					p.init();//to initilaze
					//	printf("Number of feasible particles: %d\n",p.num_feasible);
				} while (false
					||p.num_feasible<r*PSO_popsize
					);



				//	printf("Bingo\n");
				p.search(resultStore);	

				fprintf(PSOData,"Run number:  %d, Best values:  %lf  \n",i,p.gbest);
				fprintf(PSOData,"\Best position: ");
				for (int j=0;j<GANVARS;j++)
				{
					fprintf(PSOData,"%lf ",p.gbest_pos[j]);
				}
				fprintf(PSOData,"\n");
				fprintf(PSOData,"-----------------------------------------------------------------------------------------\n");
  
				//保存每一代的gbest,可以屏蔽掉下面两行。使文件简洁易读。 
				//	for (int row = 0; row <PSO_maxgen; row++)				
				//		fprintf(PSOData,"%d   %lf\n",row,resultStore[row]);	

				if (p.gbest<=Best_Value
					&&Is_feasible(p.gbest_pos)
					)
				{
					Best_Value=p.gbest;
					for (int k = 0; k < GANVARS; k++)
					{
						Best_Position[k]=p.gbest_pos[k];
					}
					excel_run_best[test_number][i]=Best_Value;
					if (excel_series)
					{
						for (int generation=0;generation<PSO_maxgen;generation++)
						{
							excel_con_series[generation]=resultStore[generation];
						}
						excel_series=false;
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
			excel_final_best[test_number]=Best_Value;

			fprintf(PSOData,"*********************************************************************************************\n");

		}

		//扫尾工作 
		delete  []resultStore;
		fclose(PSOData);


	}
	fprintf(PSO_EXCEL,"Final Best values:\n");
	for (int i=0;i<TESTTIMES;i++)
	{
		fprintf(PSO_EXCEL,"%lf\n",excel_final_best[i]);
	}
  
	fprintf(PSO_EXCEL,"RUN Best values:\n");
	for (int i=0;i<TESTTIMES;i++)
	{
		for (int j=0;j<RUNTIMES;j++)
		{
			fprintf(PSO_EXCEL,"%lf\n",excel_run_best[i][j]);
		}
		
	}
 
	fprintf(PSO_EXCEL,"Convergence series:\n");
	for (int i=0;i<PSO_maxgen;i++)
	{
		fprintf(PSO_EXCEL,"%lf\n",excel_con_series[i]);
	}
	for (int i=0;i<TESTTIMES;i++)
	{
		delete []excel_run_best[i];
	}
	delete[]excel_run_best;
	
	delete [] excel_con_series;
	delete [] excel_final_best;
	
	
	fclose(PSO_EXCEL);
	cout<<"Failed times: "<<failed_number<<endl;
	clock_t t2 = clock();
	cout<<"It takes: "<<t2-t1<<"mS"<<endl;

	system("pause");
	return 0;
}





/*
int main(void)
{
	int i=0;
	int **a=new int *[3];
	for (int i=0;i<3;i++)
	{
		a[i]=new int [5];
	}
	for (i=0;i<3;i++)
	{
		for (int j=0;j<5;j++)
		{
			a[i][j]=i+j;
		}
	}

	for (i=0;i<3;i++)
	{
		for (int j=0;j<5;j++)
		{
			cout<<a[i][j]<<' ';
		}
		cout<<endl;
	}
	for (i=0;i<3;i++)
	{
		delete []a[i];
	}
	delete[] a;
	system("pause");
	return 0;
}
 */

