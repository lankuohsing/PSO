#include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //����rand()��srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//ʹ��double �����ֵ
using namespace std;
using std::numeric_limits;		//ʹ��double �����ֵ,numeric_limits<double>::max()
#define GANVARS 2	//������ά��             *****��Ҫ����******
#define PSO_popsize 30 //���Ӹ���                    
#define PSO_maxgen 1000
#define Vmax 5 //�ٶȼ�ֵ�����ȡֵ���󣬲����ڸ�����Ż�,һ����upbound-lowbound

//����[low,uper]�����doubleֵ 
#define rnd(low,uper) ((rand()/(double)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //�������д��� 
#define TESTTIMES 20 //���Դ���
double paramater_w; //��ʷϵ��������Ȩ�أ��� �ɳ���������� 
double paramater_c1=1.49445; //��֪ϵ����һ��ȡ2.0 
double paramater_c2=1.49445; //���ϵ��, һ��ȡ2.0 
double paramate_cy = 1.0e-10;
double lowBound[GANVARS],upperBound[GANVARS]; //��Ⱥ�и���ķ�Χ��   *****��Ҫ����****** 

//ȫ�ּ�ֵ������,ע����ParticleSwarm::getGBest()���� 
//���Ӹ����ٶ���Ҫʹ�ã�����ȫ��ͨ���á� 
double gbest_pos_global[GANVARS]; 

/************************************************************************/
/* ����������                                                           */
/************************************************************************/
class Particle {
public:
	double pos[GANVARS]; 		//����λ��,position
	double v[GANVARS]; 			//�����ٶ�,velocity
	double pbest;				//���弫ֵ����Ӧֵ
	double pbest_pos[GANVARS];	//�������Ž�����꣬��Ӧÿ��֤ȯ��ϱ�����ֵ
	double fitness;				//��ǰ�����һ����Ӧֵ
public:
	double calcFitness(double pos[],int k);	//������Ӧֵ�ĺ���
	void updatePosition(int gen);				//λ�ø��º���
	void updatePBest();					//���弫ֵ���º���
};


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
	 f_value=pow(p[0]-10,3)+pow(p[1]-20,3);
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
	 double constraint[GANVARS];
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
 }

double Particle::calcFitness(double *p,int k)//������Ӧֵ�ĺ���,�������壬�����Լ������Ҫ�ӷ�����ѽ�� ****��Ҫ����*****
{

	//int k;
	double serr;
	serr=f_function(p)+h_function(k)*H_function(p);
	return serr;
}

//����λ�ƺ��ٶ�
void Particle::updatePosition(int gen)
{
	double r1,r2;
	int i=0,j=0,k=0;

	for(i=0;i<GANVARS;i++)
	{	
		paramater_w = 1.2-(double)gen/RUNTIMES;

		//�����ٶȣ��������ӵ�pbest_pos����ȫ�ֵ�  gbest_pos_global[]
		v[i] = paramater_w  * v[i] +
			paramater_c1 * rnd(0,1) * (pbest_pos[i]        -  pos[i]) +
			paramater_c2 * rnd(0,1) * (gbest_pos_global[i] -  pos[i]); 

		//�жϳ�������ٶȺ���С�ٶȡ�
		if (v[i]<-Vmax)
			v[i] = -Vmax; 
		if (v[i]>Vmax)
			v[i]=Vmax;

		//��ֵ�Ŷ���?
		r1 = rnd(0,1);
		if (r1 < paramate_cy)
		{
			r2 = rnd(0,1);
			v[i]=r2*v[i];
		}

		//����λ��
		pos[i]+=v[i]; 

		//���������� 
		if(pos[i]<lowBound[i])
			pos[i]=lowBound[i];
		if(pos[i]>upperBound[i])
			pos[i]=upperBound[i];
	}
}

//���¸��弫ֵ	
void Particle::updatePBest(){

	if(this->fitness<pbest)
	{
		pbest=this->fitness;
		for(int i=0;i<GANVARS;i++)
		{
			pbest_pos[i]=pos[i];//���¸������Ž������	
		}
	}
}

/************************************************************************/
/*        ����Ⱥ��                                                      */
/************************************************************************/
class ParticleSwarm{
public:
	double gbest; //ȫ�ּ�ֵ����Ӧֵ
	double gbest_pos[GANVARS]; //ȫ�ּ�ֵ������
	Particle PSO_pop[PSO_popsize];//�������Ӷ���Ϊ����Ⱥ�������
public: 
	void init();  //��ʼ����Ⱥ
	void getGBest(); //��ȡȫ�ּ�ֵ
	void search(double *Array); //����,col������[0,RUNTIME-1]��ָʾ�ǵڼ������� 	  
};

//��ʼ����Ⱥ����
void ParticleSwarm::init(){

	gbest=numeric_limits<double>::max();//initialize gbest with maximal value of double type
	srand((unsigned)time(NULL));//get the seed of random number

	//initialize the boundary ****��Ҫ����*****�趨���ӱ߽�
	for(int t=0;t<GANVARS;t++)
	{
		lowBound[t]=-10;
		upperBound[t]=10;
	}

	//only for the case of test problem 1
	lowBound[0]=13,upperBound[0]=100;
	lowBound[1]=0,upperBound[1]=100;

	//initialize the particle swarm
	for(int i=0;i<PSO_popsize;i++)
	{   
		//x[] store the positions��y[] store the velocity 
		double x[GANVARS];
		double y[GANVARS];

		for(int j=0;j<GANVARS;j++)
		{		
			x[j] = rnd(lowBound[j],upperBound[j]);//[lowBound[j],upperBound[j]]֮��������
			y[j] = rnd(-Vmax,Vmax);//[-Vmax,Vmax]֮�������� 
		}

		//only for the case of test problem 1
	//	x[1]=rnd(lowBound[1],upperBound[1]);
	//	x[0]=2*x[1]-1;
	


		//��ʼ��λ�ú��ٶ�	
		for(int j=0;j<GANVARS;j++){
			PSO_pop[i].pos[j]=x[j];
			PSO_pop[i].v[j]=y[j];
		}

		//calculate the fitness 
		PSO_pop[i].fitness= PSO_pop[i].calcFitness(PSO_pop[i].pos,0);

		//��ʼ��������ǰfitness����pbest
		PSO_pop[i].pbest=PSO_pop[i].fitness;  
		for(int m=0;m<GANVARS;m++)
		{
			PSO_pop[i].pbest_pos[m]=PSO_pop[i].pos[m];
		}
	}	

	getGBest();	//get gbest

	
}
//��ȡȫ�ּ�ֵ
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
				//ҲҪ����ȫ�ֵ� gbest_pos_global,����ȫ������ͨѶ�á� 
				gbest_pos_global[j]=PSO_pop[i].pos[j];		  
			}

		}
	}
}

//����
void ParticleSwarm::search(double *Array)
{
	int gen=0;//number of iterations

	//������Χ��:[0, PSO_maxgen-1],�ܹ�PSO_maxgen���� 
	while(gen<PSO_maxgen)
	{
		Array[gen] = gbest;//Array���鴢��ȫ�ּ�ֵ�ĵ���Ӧֵ

		//ÿ�����ӽ����˶�����ֵ������pbest 
		for(int k=0;k<PSO_popsize;k++)
		{
			PSO_pop[k].updatePosition(gen);
			PSO_pop[k].fitness=PSO_pop[k].calcFitness(PSO_pop[k].pos,gen);
			PSO_pop[k].updatePBest();
		}

		//update gbest.ע�������������ͨ�ŵ�ȫ�ֱ��� 
		getGBest();

		gen++;		
		
	}
}

