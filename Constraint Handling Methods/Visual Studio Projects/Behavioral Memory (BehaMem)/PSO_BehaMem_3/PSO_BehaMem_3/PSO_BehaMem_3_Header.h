#include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //包含rand()和srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//使用double 的最大值
using namespace std;
using std::numeric_limits;		//使用double 的最大值,numeric_limits<double>::max()
#define GANVARS 2	//参数的维数             *****需要更改******
#define PSO_popsize 30 //粒子个数                    
#define PSO_maxgen 1500
#define Vmax 0.5 //速度极值，如果取值过大，不利于个体的优化,一般是upbound-lowbound
#define num_constraints 2
//生成[low,uper]的随机double值 
#define rnd(low,uper) ((rand()/(double)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //程序运行次数 
#define TESTTIMES 20 //测试次数
#define flip_threshold 0.9
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
double g2(double *p);      //the 1-th constraint
double g1(double *p);      //the 2-th constraint
class Particle {
public:
	double pos[GANVARS]; 		//粒子位置,position
	double v[GANVARS]; 			//粒子速度,velocity
	double pbest;				//个体极值的适应值
	double pbest_pos[GANVARS];	//个体最优解的坐标，对应每种证券组合比例的值
	double fitness;				//当前算出的一个适应值
	double (*g[2])(double *p);

	Particle(){g[0]=g1;
	g[1]=g2;}
	double calcFitness(double pos[],int k);	//计算适应值的函数
	void updatePosition();				//位置更新函数
	void updatePBest();					//个体极值更新函数
	
	
	double f_function(double *p);//Original objective function
};



double h_function(int k);	//Penalty value 
double H_function(double *p);//Penalty factor	
double *q_funciton(double *p);//violated function of the constraints
double *g_function(double *p);//constraints
double xita_function(double q);//assignment function
int gama_function(double q);//power of the penalty function
bool Is_feasible(double *p);//To judge whether the constraints are violated


 double Particle:: f_function(double *p)   //Original objective function
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
	 g_value[0]=p[0]-2*p[1]+1;
	 g_value[1]=(p[0]*p[0]/4+p[1]*p[1]-1);
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
	 double constraint[GANVARS];
	 constraint[0]=g1(p);
	 constraint[1]=g2(p);
	 if (abs(constraint[0])<1&&constraint[1]<1)
	 {
		 return true;
	 }
	 else
	 {
		 return false;
	 }
 }
 double g1(double *p)      //the 1-th constraint
 {
	 double g1_value;
	 g1_value=(p[0]-2)*(p[0]-2)/4+p[1]*p[1]-1;
	 return g1_value;
 }
 double g2(double *p)      //the 2-th constraint
 {
	 double g2_value;
	 g2_value=(p[0]-2)*(p[0]-2)+p[1]*p[1]-2;
	 return g2_value;
 }

double Particle::calcFitness(double *p,int k)//计算适应值的函数,评估个体，如果有约束，需要加罚函数呀！ ****需要更改*****
{

	//int k;
	double serr;
	serr=f_function(p);
	return serr;
}

//更新位移和速度
void Particle::updatePosition()
{
	double r1,r2;
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
	double gbest; //全局极值的适应值
	double gbest_pos[GANVARS]; //全局极值的坐标
	Particle PSO_pop[PSO_popsize];//单个粒子定义为粒子群类的属性
public: 
	void init();  //初始化种群
	void getGBest(); //获取全局极值
	void search(double *Array); //迭代,col参数是[0,RUNTIME-1]，指示是第几次运行 	  
};

//初始化种群函数
void ParticleSwarm::init(){
	int num_freasible=0;
	gbest=numeric_limits<double>::max();//初始化全局极值的适应值为无穷大
	srand((unsigned)time(NULL));//获取随机数种子

	//初始化Bound ****需要更改*****设定粒子边界
	
	lowBound[0]=0,upperBound[0]=4;
	lowBound[1]=-2,upperBound[1]=2;
	while (true)
	{
		num_freasible=0;
		//初始化种群,包括位置，速度，fitness和pbest 
		for(int i=0;i<PSO_popsize;i++)
		{   
			//x[] store the positions,y[] store the velocity 
			double x[GANVARS];
			double y[GANVARS];

			for(int j=0;j<GANVARS;j++)
			{		
				x[j] = rnd(lowBound[j],upperBound[j]);//[lowBound[j],upperBound[j]]之间的随机数
				y[j] = rnd(-Vmax,Vmax);//[-Vmax,Vmax]之间的随机数 
			}
			x[0]=2*x[1]-1;

			//初始化位置和速度	
			for(int j=0;j<GANVARS;j++){
				PSO_pop[i].pos[j]=x[j];
				PSO_pop[i].v[j]=y[j];
			}

			//计算fitness 
			PSO_pop[i].fitness= PSO_pop[i].g[0](PSO_pop[i].pos);
			if ((PSO_pop[i].fitness)<=2.6)
			{
				num_freasible++;
			}



			//初始化，将当前fitness赋给pbest
			PSO_pop[i].pbest=PSO_pop[i].fitness;  
			for(int m=0;m<GANVARS;m++)
			{
				PSO_pop[i].pbest_pos[m]=PSO_pop[i].pos[m];
			}
		}
		
		if (num_freasible>=PSO_popsize*flip_threshold)
		{
			num_freasible=0;
			break;
		}
		else
		{
			num_freasible=0;
		}
		
	}

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
	int i,j;
	int num_freasible=0;
	
	for (i=0;i<(num_constraints-1);i++)
	{
		for ( j = 0; j < PSO_popsize; j++)
		{
			PSO_pop[j].fitness=PSO_pop[j].g[i+1](PSO_pop[j].pos);
			PSO_pop[j].pbest=PSO_pop[j].fitness;  
			for(int m=0;m<GANVARS;m++)
			{
				PSO_pop[j].pbest_pos[m]=PSO_pop[j].pos[m];
			}
		}
		gbest=numeric_limits<double>::max();
		getGBest();

		while (true)
		{
			
			for (j=0;j<PSO_popsize;j++)
			{
				PSO_pop[j].updatePosition();
				PSO_pop[j].fitness=PSO_pop[j].g[i+1](PSO_pop[j].pos);
				PSO_pop[j].updatePBest();
				if ((PSO_pop[j].fitness)<=-1)
				{
					num_freasible++;
				}

			}
			getGBest();
			
			if (num_freasible>=PSO_popsize*flip_threshold)
			{
				
				num_freasible=0;
				break;
			}
			else
			{
				num_freasible=0;
			}
			
		}

	}
	for ( j = 0; j < PSO_popsize; j++)
	{
		PSO_pop[j].fitness=PSO_pop[j].f_function(PSO_pop[j].pos);
		PSO_pop[j].pbest=PSO_pop[j].fitness;  
		for(int m=0;m<GANVARS;m++)
		{
			PSO_pop[j].pbest_pos[m]=PSO_pop[j].pos[m];
		}
	}
	gbest=numeric_limits<double>::max();
	getGBest();
	int gen=0;
	for (gen=0;gen<PSO_maxgen;gen++)
	{
		Array[gen] = gbest;//Array数组储存全局极值的的适应值
		for (j=0;j<PSO_popsize;j++)
		{
			PSO_pop[j].updatePosition();
			PSO_pop[j].fitness=PSO_pop[j].f_function(PSO_pop[j].pos);
			PSO_pop[j].updatePBest();
		}
	}
	getGBest(); 
	
}

