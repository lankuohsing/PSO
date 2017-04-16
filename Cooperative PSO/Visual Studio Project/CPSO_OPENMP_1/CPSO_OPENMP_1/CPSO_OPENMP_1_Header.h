#include <omp.h>											  #include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //包含rand()和srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//使用double 的最大值
using namespace std;
using std::numeric_limits;		//使用double 的最大值,numeric_limits<double>::max()
#define GANVARS 12	//参数的维数             *****需要更改******
#define GROUPE_MEMBER 20
#define PSO_popsize GANVARS*GROUPE_MEMBER //粒子个数                    
#define PSO_maxgen 1000
#define Vmax 2.5 //速度极值，如果取值过大，不利于个体的优化,一般是upbound-lowbound

//生成[low,uper]的随机double值 
#define rnd(low,uper) ((rand()/(double)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //程序运行次数 
#define TESTTIMES 10 //测试次数
double paramater_w; //历史系数（惯性权重）， 由程序随机生成 
double paramater_c1=1.49445; //认知系数，一般取2.0 / 1.49445
double paramater_c2=1.49445; //社会系数, 一般取2.0 / 1.49445
double paramate_cy = 1.0e-10;
double lowBound[GANVARS],upperBound[GANVARS]; //种群中个体的范围，   *****需要更改******


//全局极值的坐标,注意由ParticleSwarm::getGBest()更新 
//粒子更新速度需要使用，用于全局通信用。 
//double gbest_pos_global[GANVARS]; 
double COOP_gbest_pos_global[GANVARS][GANVARS];
double COOP_pos_global[GANVARS];

/************************************************************************/
/* 单个粒子类                                                           */
/************************************************************************/
class Particle {
public:
	double pos[GANVARS]; 		//粒子位置,position
	double v[GANVARS]; 			//粒子速度,velocity
	double pbest;				//个体极值的适应值
	double pbest_pos[GANVARS];	//个体最优解的坐标，对应每种证券组合比例的值
	double fitness;				//当前算出的一个适应值
public:
	double calcFitness(double pos[],int k);	//计算适应值的函数
	void updatePosition(int i,int gen);				//位置更新函数
	void updatePBest();					//个体极值更新函数
};


double f_function(double *p);//Original objective function
double h_function(int k);	//Penalty value 
double H_function(double *p);//Penalty factor	
double *q_funciton(double *p);//violated function of the constraints
double *g_function(double *p);//constraints
double xita_function(double q);//assignment function
int gama_function(double q);//power of the penalty function
bool Is_feasible(double *p);//To judge whether the constraints are violated
void COOP_UpdatePosition(double *p);//cooperative update



double f_function(double *p)   //Original objective function
{
	double f_value;
	f_value=(p[0]-1)*(p[0]-1)+(p[1]-2)*(p[1]-2)+(p[2]-3)*(p[2]-3)
		+(p[3]-4)*(p[3]-4)+(p[4]-5)*(p[4]-5)+(p[5]-6)*(p[5]-6)
		+(p[6]-7)*(p[6]-7)+(p[7]-8)*(p[7]-8)+(p[8]-9)*(p[8]-9)
		+(p[9]-10)*(p[9]-10)+(p[10]-11)*(p[10]-11)+(p[11]-12)*(p[11]-12);
	return f_value;

}
double h_function(int k)	//Penalty value
{
	double h_value;
	h_value=pow(1.5,k);
	return h_value;

}
double H_function(double *p)//Penalty factor
{
	int i;
	double *g,*q,gama_value,xita_value,pow_value;
	g=g_function(p);
	q=q_funciton(g);

	double H_value=0.0;
	for (i=0;i<GANVARS;i++)
	{
		xita_value=xita_function(q[i]);
		gama_value=gama_function(q[i]);
		pow_value=pow(q[i],gama_value);
		H_value+=xita_value*pow_value; 
	}
	delete []q;

	return H_value;
}
double *q_funciton(double *g)//violated function of the constraints
{
	double *q_value=new double[GANVARS];
	int i=0;
	for (i=0;i<GANVARS;i++)
	{
		q_value[i]=(g[i]>=0?g[i]:0);
	}

	//	 q_value[0]=abs(g[0]);
	//	 q_value[1]=(g[1]>=0?g[1]:0);
	delete []g;
	return q_value;
}
double *g_function(double *p)//constraints
{
	double *g_value=new double[GANVARS];
	g_value[0]=100-pow(p[0]-5,2)-pow(p[1]-5,2);
	g_value[1]=pow(p[0]-6,2)+pow(p[1]-5,3)-82.81;
	return	g_value;
}
double xita_function(double q)//assignment function
{
	double xita_value;
	if (q<0.001)
	{
		xita_value=10;
	}
	else
	{
		if (q<=0.1)
		{
			xita_value=20;
		}
		else
		{
			if (q<=1)
			{
				xita_value=100;
			}
			else
			{
				xita_value=300;
			}
		}
	}
	return xita_value;
}
int gama_function(double q)//power of the penalty function
{
	int gama_value;
	if (q<1)
	{
		gama_value=1;
	}
	else
	{
		gama_value=2;
	}
	return gama_value;
}
bool Is_feasible(double *p)//To judge whether the constraints are violated
{
/*	double constraint[GANVARS];
	constraint[0]=100-(p[0]-5)*(p[0]-5)-(p[1]-5)*(p[1]-5);
	constraint[1]=(p[0]-6)*(p[0]-6)+(p[1]-5)*(p[1]-5)-82.81;
	if (constraint[0]<1.0&&constraint[1]<1.0)
	{
		return true;
	}
	else
	{
		return false;
	}
	*/
	return true;
}


double Particle::calcFitness(double *p,int k)//计算适应值的函数,评估个体，如果有约束，需要加罚函数呀！ ****需要更改*****
{

	//int k;
	double serr;
	serr=f_function(p);
	return serr;
}

//更新位移和速度
void Particle::updatePosition(int j,int gen)
{
	double r1,r2;
	int i=0;
	#pragma omp parallel for
	for(i=0;i<GANVARS;i++)
	{	
		paramater_w = 1-(double)gen/RUNTIMES;

		//更新速度，利用粒子的pbest_pos，和全局的  gbest_pos_global[]
		v[i] = paramater_w  * v[i] +
			paramater_c1 * rnd(0,1) * (pbest_pos[i]        -  pos[i]) +
			paramater_c2 * rnd(0,1) * (COOP_pos_global[i] -  pos[i]); 

		//判断超出最大速度和最小速度。
		if (v[i]<-Vmax)
			v[i] = -Vmax; 
		if (v[i]>Vmax)
			v[i]=Vmax;

		//极值扰动吗?
		r1 = rnd(0,1);
		if (r1 < paramate_cy)
		{
			r2 = rnd(0,1);
			v[i]=r2*v[i];
		}

		//更新位移
		pos[i]+=v[i]; 

		//出界则拉回 
		if(pos[i]<lowBound[i])
			pos[i]=lowBound[i];
		if(pos[i]>upperBound[i])
			pos[i]=upperBound[i];
	}
}

//更新个体极值	
void Particle::updatePBest(){

	if(this->fitness<pbest)
	{
		pbest=this->fitness;
		for(int i=0;i<GANVARS;i++)
		{
			pbest_pos[i]=pos[i];//更新个体最优解的坐标	
		}
	}
}

/************************************************************************/
/*        粒子群类                                                      */
/************************************************************************/
class ParticleSwarm{
public:
	double gbest; //全局极值的适应值
	double gbest_pos[GANVARS]; //全局极值的坐标
	double COOP_gbest[GANVARS];//global best fitness of each swarm 
	double COOP_gbest_pos[GANVARS][GANVARS];//global best position of each swarm 
	Particle PSO_pop[PSO_popsize];//单个粒子定义为粒子群类的属性
	Particle ***b;//context vector
	double COOP_pos[GANVARS]; //vector b

public: 
	void init();  //初始化种群
	void getGBest(); //获取全局极值
	void search(double *Array); //
	void Coop_search(double *Array);//Cooperative PSO search function
	void init_b();//To update vector b
	~ParticleSwarm(){
		for (int i=0;i<GANVARS;i++)
		{
		delete []b[i];
		}
		delete []b;
	}
};

//初始化种群函数
void ParticleSwarm::init(){

	////////////////////////initialize gbests with maximum of type double ///////////////////////////////
	gbest=numeric_limits<double>::max();
	for (int i=0;i<GANVARS;i++)
	{
		COOP_gbest[i]=numeric_limits<double>::max();
	}
	////////////////////////initialize gbests with maximum of type double ///////////////////////////////

	////////////////////To allocate space for vector b///////////////////////////////	
	b=new Particle **[GANVARS];
	for (int i=0;i<GANVARS;i++)
	{
		b[i]=new Particle *[GROUPE_MEMBER];
	}
	////////////////////To allocate space for vector b///////////////////////////////////

	////////////////////initialize the state of each particle////////////////////////////
	srand((unsigned)time(NULL));//get the seed of random number
	//initialize the boundary ****需要更改*****设定粒子边界
	for(int t=0;t<GANVARS;t++)
	{
		lowBound[t]=-100;
		upperBound[t]=100;
	}

	//to initialize the particle swarm
//	#pragma omp parallel for
	for(int i=0;i<PSO_popsize;i++)
	{   
		//x[] store the positions。y[] store the velocity 
		double x[GANVARS];
		double y[GANVARS];

		for(int j=0;j<GANVARS;j++)
		{		
			x[j] = rnd(lowBound[j],upperBound[j]);//[lowBound[j],upperBound[j]]之间的随机数
			y[j] = rnd(-Vmax,Vmax);//[-Vmax,Vmax]之间的随机数 
		}

		//初始化位置和速度	
		for(int j=0;j<GANVARS;j++){
			PSO_pop[i].pos[j]=x[j];
			PSO_pop[i].v[j]=y[j];
		}

		//calculate the fitness of each particle 
		PSO_pop[i].fitness= PSO_pop[i].calcFitness(PSO_pop[i].pos,0);

		/////////////////initialize pbest of each particle///////////////////////////////////
		PSO_pop[i].pbest=PSO_pop[i].fitness;  
		for(int m=0;m<GANVARS;m++)
		{
			PSO_pop[i].pbest_pos[m]=PSO_pop[i].pos[m];
		}
		/////////////////initialize pbest of each particle///////////////////////////////////

		//to get vector b
		b[i/GROUPE_MEMBER][i%GROUPE_MEMBER]=&(PSO_pop[i]);

	}	
	////////////////////initialize the state of each particle////////////////////////////	

	init_b();
	


}
//获取全局极值
void ParticleSwarm::getGBest()
{
	
	for (int i=0;i<GANVARS;i++)
	{
		//
	}
}

//Toe update vector b
void ParticleSwarm::init_b()
{
	for (int i=0;i<GANVARS;i++)
	{
		for (int j=0;j<GROUPE_MEMBER;j++)
		{
			if (b[i][j]->fitness<COOP_gbest[i])
			{
				COOP_gbest[i]=b[i][j]->fitness;//global best of each swarm
				for (int k=0;k<GANVARS;k++)
				{
					COOP_gbest_pos[i][k]=b[i][j]->pos[k];//global best position of each swarm
					COOP_gbest_pos_global[i][k]=b[i][j]->pos[k]; //
				}
			}
		}
		COOP_pos[i]=COOP_gbest_pos[i][i];//best component
		COOP_pos_global[i]=COOP_pos[i];

	}

}

//迭代
void ParticleSwarm::search(double *Array)
{
	int gen=0;//number of iterations
	double temp1;
	double temp_COOP=numeric_limits<double>::max();

	////////////////////////////Copy best component///////////////////////////////////////
	double *temp_COOP_pos=new double [GANVARS];
	for (int i=0;i<GANVARS;i++)
	{
		temp_COOP_pos[i]=COOP_pos[i];
	}
	////////////////////////////Copy best component///////////////////////////////////////

	
	//代数范围是:[0, PSO_maxgen-1],总共PSO_maxgen代。 
	
	while(gen<PSO_maxgen)
	{
	//	#pragma omp parallel for
		for (int i=0;i<GANVARS;i++)
		{
		//	#pragma omp parallel for
			for (int j=0;j<GROUPE_MEMBER;j++)
			{
				temp_COOP_pos[i]=b[i][j]->pos[i];//substitute each particle of swarm j into b
				
				temp1=b[i][j]->calcFitness(temp_COOP_pos,gen);//

//////////////////////////////update cooperative gbest and position//////////////////////////////////////
				if (temp1<temp_COOP)
				{
					temp_COOP=temp1;
					COOP_pos[i]=temp_COOP_pos[i];
					
				}
//////////////////////////////update cooperative gbest and position//////////////////////////////////////

				
//////////////////////////////update personal best//////////////////////////////////////
				double temp2=b[i][j]->calcFitness(b[i][j]->pos,gen);
				if(temp2<b[i][j]->pbest)
				{
					b[i][j]->pbest=temp2;
					for(int k=0;k<GANVARS;k++)
					{
						b[i][j]->pbest_pos[k]=b[i][j]->pos[k];//更新个体最优解的坐标
						COOP_gbest_pos_global[i][k]=b[i][j]->pos[k];
					}
				}
//////////////////////////////update personal best//////////////////////////////////////
				
				
	
			}
		COOP_pos_global[i]=COOP_pos[i];
		for (int L=0;L<GROUPE_MEMBER;L++)
		{
			b[i][L]->updatePosition(i,gen);
		}
		
		}
	gbest=Array[gen]=b[0][0]->calcFitness(COOP_pos_global,gen);
	gen++;		

	}
	delete []temp_COOP_pos;
}

