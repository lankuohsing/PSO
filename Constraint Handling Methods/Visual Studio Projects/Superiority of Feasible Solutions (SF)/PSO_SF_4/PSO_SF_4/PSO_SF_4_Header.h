#include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //包含rand()和srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//使用double 的最大值
using namespace std;
using std::numeric_limits;		//使用double 的最大值,numeric_limits<double>::max()
#define GANVARS 2	//参数的维数             *****需要更改******
#define PSO_popsize 20 //粒子个数                    
#define PSO_maxgen 2000
#define Vmax 0.1 //速度极值，如果取值过大，不利于个体的优化,一般是upbound-lowbound
#define  r 0.1
//生成[low,uper]的随机double值 
#define rnd(low,uper) ((rand()/(double)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //程序运行次数 
#define TESTTIMES 20 //测试次数
double paramater_w; //历史系数（惯性权重）， 由程序随机生成 
double paramater_c1=1.49445; //认知系数，一般取2.0 
double paramater_c2=1.49445; //社会系数, 一般取2.0 
double paramate_cy = 1.0e-10;
double lowBound[GANVARS],upperBound[GANVARS]; //种群中个体的范围，   *****需要更改****** 

//全局极值的坐标,注意由ParticleSwarm::getGBest()更新 
//粒子更新速度需要使用，用于全局通信用。 
double gbest_pos_global[GANVARS]; 

/************************************************************************/
/* 单个粒子类                                                           */
/************************************************************************/



double f_function(double *p);//Original objective function
double h_function(int k);	//Penalty value 
double H_function(double *p);//Penalty factor	
double *q_funciton(double *p);//violated function of the constraints
double *g_function(double *p);//constraints
double xita_function(double q);//assignment function
int gama_function(double q);//power of the penalty function
bool Is_feasible(double *p);//To judge whether the constraints are violated



 double f_function(double *p)   //Original objective function
 {
	 double f_value;
	 f_value=pow(p[0]-2,2)+pow(p[1]-2,2);
	 return f_value;

 }
 double h_function(int k)	//Penalty value
 {
	 double h_value;
	 h_value=sqrt(k);
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
	 g_value[0]=(p[0]-2)*(p[0]-2)/4+p[1]*p[1]-1;
	 g_value[1]=(p[0]-2)*(p[0]-2)+p[1]*p[1]-2;
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
	 double *constraint;
	 constraint=g_function(p);
	 
	 if ((constraint[0])<=1&&constraint[1]<=1)
	 {
		 delete[] constraint;
		 return true;
	 }
	 else
	 {
		 delete[] constraint;
		 return false;
	 }
 }
 bool Is_feasible2(double *p)//To judge whether the constraints are violated
 {
	 double *constraint;
	 constraint=g_function(p);

	 if ((constraint[0])<=5000&&constraint[1]<=5000)
	 {
		 delete[] constraint;
		 return true;
	 }
	 else
	 {
		 delete[] constraint;
		 return false;
	 }
 }
 

 class Particle {
 public:
	 double pos[GANVARS]; 		//粒子位置,position
	 double v[GANVARS]; 		//粒子速度,velocity
	 double pbest;				//个体极值的适应值
	 double pbest_pos[GANVARS];	//个体最优解的坐标，对应每种证券组合比例的值
	 double fitness;			//当前算出的一个适应值
	 double violation;			//constraint violation value
 public:
	 double calcFitness(double pos[],int k);	//计算适应值的函数
	 void updatePosition(int gen);				//位置更新函数
	 void updatePBest();					//个体极值更新函数
	 bool Is_Replace(double *p1,double *p2);//to judge whether to replace the position
	 double constraint_violation_value(double *p);//to calculate the constraint violation value of an infeasible solution
 };
 double Particle::constraint_violation_value(double *p)//to calculate the constraint violation value of an infeasible solution
 {
	 double violation_value=0.0;
	 double *g;
	 g=g_function(p);
	 violation_value+=abs(g[0]);
	 violation_value+=(g[1]>0)?g[1]:0;
	 delete []g;
	 return violation_value;
 }

 bool Particle::Is_Replace(double *p1,double *p2)//to judge whether to replace the position
 {
	if (!Is_feasible(p1)&&Is_feasible(p2))
	{
		return true;
	}
	else
	{
		if (Is_feasible(p1)&&Is_feasible(p2))
		{
			if (f_function(p1)>f_function(p2))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			if (!Is_feasible(p1)&&!Is_feasible(p2))
			{
				if (constraint_violation_value(p1)>constraint_violation_value(p2))
				{
					return true;
				}
				else
				{
					 return false;
				}
			}
			else
			{
				return false;
			}
		}
	}
 }

double Particle::calcFitness(double *p,int k)//计算适应值的函数,评估个体，如果有约束，需要加罚函数呀！ ****需要更改*****
{

	//int k;
	double serr;
	serr=f_function(p);
	return serr;
}

//更新位移和速度
void Particle::updatePosition(int gen)
{
	double r1,r2;
	int i=0,j=0,k=0;
	double temp_pos[GANVARS],temp_vol[GANVARS];
	for(i=0;i<GANVARS;i++)
	{	
		temp_pos[i]=pos[i];
		temp_vol[i]=pos[i];
		paramater_w = 1.0-(double)gen/(double)PSO_maxgen;

		//更新速度，利用粒子的pbest_pos，和全局的  gbest_pos_global[]
		temp_vol[i] = paramater_w  * v[i] +
			paramater_c1 * rnd(0,1) * (pbest_pos[i]        -  pos[i]) +
			paramater_c2 * rnd(0,1) * (gbest_pos_global[i] -  pos[i]); 

		//判断超出最大速度和最小速度。
		if (temp_vol[i]<-Vmax)
			temp_vol[i] = -Vmax; 
		if (temp_vol[i]>Vmax)
			temp_vol[i]=Vmax;

		//极值扰动吗?
		r1 = rnd(0,1);
		if (r1 < paramate_cy)
		{
			r2 = rnd(0,1);
			temp_vol[i]=r2*temp_vol[i];
		}

		//更新位移
		temp_pos[i]+=v[i]; 

		//出界则拉回 
		if(temp_pos[i]<lowBound[i])
			temp_pos[i]=lowBound[i];
		if(temp_pos[i]>upperBound[i])
			temp_pos[i]=upperBound[i];
	}
	
	if (Is_Replace(pos,temp_pos))

	{
		for (i = 0; i < GANVARS; i++)
		{
			v[i]=temp_vol[i];
			pos[i]=temp_pos[i];
		}
		this->fitness=this->calcFitness(this->pos,gen);
		this->pbest=this->fitness;
		for(i=0;i<GANVARS;i++)
		{
			this->pbest_pos[i]=this->pos[i];//更新个体最优解的坐标	
		}
		
		
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
	Particle PSO_pop[PSO_popsize];//单个粒子定义为粒子群类的属性
	int FPS[PSO_popsize];//feasible	particle set
	int num_feasible;
public: 
	void init();  //初始化种群
	void getGBest(); //获取全局极值
	void search(double *Array); //迭代,col参数是[0,RUNTIME-1]，指示是第几次运行 
	
};

//初始化种群函数
void ParticleSwarm::init(){
	int i,j,k,t;
	num_feasible=0;	//initialize the number of feasible particles with 0
	gbest=numeric_limits<double>::max();//初始化全局极值的适应值为无穷大
	srand((unsigned)time(NULL));//获取随机数种子

	//初始化Bound ****需要更改*****设定粒子边界
	
	lowBound[0]=-100,upperBound[0]=100;
	lowBound[1]=-100,upperBound[1]=100;

	//初始化种群,包括位置，速度，fitness和pbest 
	for(i=0;i<PSO_popsize;i++)
	{   
		FPS[i]=-1;//to initialize FPS
		//x[] 保存位置。y[]保存速度 
		double x[GANVARS];
		double y[GANVARS];
		
		
		for(int j=0;j<GANVARS;j++)
		{		
			x[j] = rnd(lowBound[j],upperBound[j]);//[lowBound[j],upperBound[j]]之间的随机数
			y[j] = rnd(-Vmax,Vmax);//[-Vmax,Vmax]之间的随机数 
		}

		//初始化位置和速度	
		for(j=0;j<GANVARS;j++){
			PSO_pop[i].pos[j]=x[j];
			PSO_pop[i].v[j]=y[j];
		}

		//to calculate the constraint violation value
		PSO_pop[i].violation=PSO_pop[i].constraint_violation_value(x);

		//计算fitness 
		PSO_pop[i].fitness= PSO_pop[i].calcFitness(PSO_pop[i].pos,0);

		//初始化，将当前fitness赋给pbest								//mark   
		PSO_pop[i].pbest=PSO_pop[i].fitness;
		//initialize pbest_pos with current pos
		for(int m=0;m<GANVARS;m++)
		{
			PSO_pop[i].pbest_pos[m]=PSO_pop[i].pos[m];
		}
		if (Is_feasible2(PSO_pop[i].pos))
		{
			FPS[num_feasible]=i;
			num_feasible++;
		}
	}	
	
	
	getGBest();

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
void ParticleSwarm::search(double *Array)
{
	int gen=0;//迭代次数循环变量

	//代数范围是:[0, PSO_maxgen-1],总共PSO_maxgen代。 
	while(gen<PSO_maxgen)
	{
		Array[gen] = gbest;//Array数组储存全局极值的的适应值

		//每个粒子进行运动，求值，更行pbest 
		for(int k=0;k<PSO_popsize;k++)
		{
			PSO_pop[k].updatePosition(gen);
			
			
		}

		//更新gbest.注意包含更新用于通信的全局变量 
		getGBest();

		gen++;		
		
	}
}

