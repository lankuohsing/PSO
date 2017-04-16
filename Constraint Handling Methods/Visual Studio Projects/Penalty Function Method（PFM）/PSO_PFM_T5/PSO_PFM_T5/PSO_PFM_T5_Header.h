#include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //包含rand()和srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//使用int 的最大值
using namespace std;
using std::numeric_limits;		//使用int 的最大值,numeric_limits<int>::max()
#define GANVARS 2	//参数的维数             *****需要更改******
#define PSO_popsize 30 //粒子个数                    
#define PSO_maxgen 500
#define Vmax 0.5 //速度极值，如果取值过大，不利于个体的优化,一般是upbound-lowbound

//生成[low,uper]的随机int值 
#define rnd(low,uper) ((rand()/(int)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //程序运行次数 
#define TESTTIMES 20 //测试次数
int paramater_w; //历史系数（惯性权重）， 由程序随机生成 
int paramater_c1=1.49445; //认知系数，一般取2.0 
int paramater_c2=1.49445; //社会系数, 一般取2.0 
int paramate_cy = 1.0e-10;
int lowBound[GANVARS],upperBound[GANVARS]; //种群中个体的范围，   *****需要更改****** 

//全局极值的坐标,注意由ParticleSwarm::getGBest()更新 
//粒子更新速度需要使用，用于全局通信用。 
int gbest_pos_global[GANVARS]; 

/************************************************************************/
/* 单个粒子类                                                           */
/************************************************************************/
class Particle {
public:
	int pos[GANVARS]; 		//粒子位置,position
	int v[GANVARS]; 			//粒子速度,velocity
	int pbest;				//个体极值的适应值
	int pbest_pos[GANVARS];	//个体最优解的坐标，对应每种证券组合比例的值
	int fitness;				//当前算出的一个适应值
public:
	int calcFitness(int pos[],int k);	//计算适应值的函数
	void updatePosition();				//位置更新函数
	void updatePBest();					//个体极值更新函数
};


int f_function(int *p);//Original objective function
int h_function(int k);	//Penalty value 
int H_function(int *p);//Penalty factor	
int *q_funciton(int *p);//violated function of the constraints
int *g_function(int *p);//constraints
int xita_function(int q);//assignment function
int gama_function(int q);//power of the penalty function
bool Is_feasible(int *p);//To judge whether the constraints are violated


int f_function(int *p)   //Original objective function
{
	int f_value;
	f_value=(p[0]-2)*(p[0]-2)+(p[1]-2)*(p[1]-2);
	return f_value;

}
int h_function(int k)	//Penalty value
{
	int h_value;
	h_value=pow(2,k);
	return h_value;

}
int H_function(int *p)//Penalty factor
{
	int i;
	int *g,*q,gama_value,xita_value,pow_value;
	g=g_function(p);
	q=q_funciton(g);

	int H_value=0.0;
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
int *q_funciton(int *g)//violated function of the constraints
{
	int *q_value=new int[GANVARS];
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
int *g_function(int *p)//constraints
{
	int *g_value=new int[GANVARS];
	g_value[0]=(p[0]-2)*(p[0]-2)/4+(p[1]-0.4)*(p[1]-0.4)-1;
	g_value[1]=(p[0]-2)*(p[0]-2)+(p[1]-0.4)*(p[1]-0.4)-2;
//	g_value[1]=(p[0]-2)*(p[0]-2)/4+p[1]*p[1]-1;
	return	g_value;
}
int xita_function(int q)//assignment function
{
	int xita_value;
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
	xita_value*=10;
	return xita_value;
}
int gama_function(int q)//power of the penalty function
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
bool Is_feasible(int *p)//To judge whether the constraints are violated
{
	int constraint[GANVARS];
	constraint[0]=(p[0]-2)*(p[0]-2)/4+p[1]*p[1]-1;
	constraint[1]=(p[0]-2)*(p[0]-2)+p[1]*p[1]-2;
//	constraint[1]=(p[0]-2)*(p[0]-2)/4+p[1]*p[1]-1;
	if ((constraint[0])<0.1&&constraint[1]<0.1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int Particle::calcFitness(int *p,int k)//计算适应值的函数,评估个体，如果有约束，需要加罚函数呀！ ****需要更改*****
{

	//int k;
	int serr;
	serr=f_function(p)+h_function(k)*H_function(p);
	return serr;
}

//更新位移和速度
void Particle::updatePosition()
{
	int r1,r2;
	int i=0,j=0,k=0;

	for(i=0;i<GANVARS;i++)
	{	
		paramater_w = rnd(0,1);

		//更新速度，利用粒子的pbest_pos，和全局的  gbest_pos_global[]
		v[i] = paramater_w  * v[i] +
			paramater_c1 * rnd(0,1) * (pbest_pos[i]        -  pos[i]) +
			paramater_c2 * rnd(0,1) * (gbest_pos_global[i] -  pos[i]); 

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
	int gbest; //全局极值的适应值
	int gbest_pos[GANVARS]; //全局极值的坐标
	Particle PSO_pop[PSO_popsize];//单个粒子定义为粒子群类的属性
public: 
	void init();  //初始化种群
	void getGBest(); //获取全局极值
	void search(int *Array); //迭代,col参数是[0,RUNTIME-1]，指示是第几次运行 	  
};

//初始化种群函数
void ParticleSwarm::init(){

	gbest=numeric_limits<int>::max();//initialize gbest with maximal value of int type
	srand((unsigned)time(NULL));//get the seed of random number

	//initialize the boundary ****需要更改*****设定粒子边界
	for(int t=0;t<GANVARS;t++)
	{
		lowBound[t]=-10;
		upperBound[t]=10;
	}

	//only for the case of test problem 1
	lowBound[0]=0,upperBound[0]=4;
	lowBound[1]=-2,upperBound[1]=2;

	//initialize the particle swarm
	for(int i=0;i<PSO_popsize;i++)
	{   
		//x[] store the positions。y[] store the velocity 
		int x[GANVARS];
		int y[GANVARS];

		for(int j=0;j<GANVARS;j++)
		{		
			x[j] = rnd(lowBound[j],upperBound[j]);//[lowBound[j],upperBound[j]]之间的随机数
			y[j] = rnd(-Vmax,Vmax);//[-Vmax,Vmax]之间的随机数 
		}

		//only for the case of test problem 1
		//	x[1]=rnd(lowBound[1],upperBound[1]);
		//	x[0]=2*x[1]-1;



		//初始化位置和速度	
		for(int j=0;j<GANVARS;j++){
			PSO_pop[i].pos[j]=x[j];
			PSO_pop[i].v[j]=y[j];
		}

		//calculate the fitness 
		PSO_pop[i].fitness= PSO_pop[i].calcFitness(PSO_pop[i].pos,0);

		//初始化，将当前fitness赋给pbest
		PSO_pop[i].pbest=PSO_pop[i].fitness;  
		for(int m=0;m<GANVARS;m++)
		{
			PSO_pop[i].pbest_pos[m]=PSO_pop[i].pos[m];
		}
	}	

	getGBest();	//get gbest


}
//获取全局极值
void ParticleSwarm::getGBest()
{
	for(int i=0;i<PSO_popsize;i++)
	{
		if(PSO_pop[i].fitness<gbest)
		{
			gbest=PSO_pop[i].fitness;
			for(int j=0;j<GANVARS;j++)
			{
				gbest_pos[j]=PSO_pop[i].pos[j];	
				//也要更新全局的 gbest_pos_global,用于全体粒子通讯用。 
				gbest_pos_global[j]=PSO_pop[i].pos[j];		  
			}

		}
	}
}

//迭代
void ParticleSwarm::search(int *Array)
{
	int gen=0;//number of iterations

	//代数范围是:[0, PSO_maxgen-1],总共PSO_maxgen代。 
	while(gen<PSO_maxgen)
	{
		Array[gen] = gbest;//Array数组储存全局极值的的适应值

		//每个粒子进行运动，求值，更行pbest 
		for(int k=0;k<PSO_popsize;k++)
		{
			PSO_pop[k].updatePosition();
			PSO_pop[k].fitness=PSO_pop[k].calcFitness(PSO_pop[k].pos,gen);
			PSO_pop[k].updatePBest();
		}

		//update gbest.注意包含更新用于通信的全局变量 
		getGBest();

		gen++;		

	}
}

