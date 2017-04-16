#include <iostream> //include cout cin 
#include <cstdio>   //include printf() 
#include <cstdlib>  //����rand()��srand()
#include <cmath>    //include cos(),sin() 
#include <ctime>    //include time()
#include <limits>	//ʹ��int �����ֵ
using namespace std;
using std::numeric_limits;		//ʹ��int �����ֵ,numeric_limits<int>::max()
#define GANVARS 2	//������ά��             *****��Ҫ����******
#define PSO_popsize 30 //���Ӹ���                    
#define PSO_maxgen 500
#define Vmax 0.5 //�ٶȼ�ֵ�����ȡֵ���󣬲����ڸ�����Ż�,һ����upbound-lowbound

//����[low,uper]�����intֵ 
#define rnd(low,uper) ((rand()/(int)RAND_MAX)*((uper)-(low))+(low))
#define RUNTIMES 30   //�������д��� 
#define TESTTIMES 20 //���Դ���
int paramater_w; //��ʷϵ��������Ȩ�أ��� �ɳ���������� 
int paramater_c1=1.49445; //��֪ϵ����һ��ȡ2.0 
int paramater_c2=1.49445; //���ϵ��, һ��ȡ2.0 
int paramate_cy = 1.0e-10;
int lowBound[GANVARS],upperBound[GANVARS]; //��Ⱥ�и���ķ�Χ��   *****��Ҫ����****** 

//ȫ�ּ�ֵ������,ע����ParticleSwarm::getGBest()���� 
//���Ӹ����ٶ���Ҫʹ�ã�����ȫ��ͨ���á� 
int gbest_pos_global[GANVARS]; 

/************************************************************************/
/* ����������                                                           */
/************************************************************************/
class Particle {
public:
	int pos[GANVARS]; 		//����λ��,position
	int v[GANVARS]; 			//�����ٶ�,velocity
	int pbest;				//���弫ֵ����Ӧֵ
	int pbest_pos[GANVARS];	//�������Ž�����꣬��Ӧÿ��֤ȯ��ϱ�����ֵ
	int fitness;				//��ǰ�����һ����Ӧֵ
public:
	int calcFitness(int pos[],int k);	//������Ӧֵ�ĺ���
	void updatePosition();				//λ�ø��º���
	void updatePBest();					//���弫ֵ���º���
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

int Particle::calcFitness(int *p,int k)//������Ӧֵ�ĺ���,�������壬�����Լ������Ҫ�ӷ�����ѽ�� ****��Ҫ����*****
{

	//int k;
	int serr;
	serr=f_function(p)+h_function(k)*H_function(p);
	return serr;
}

//����λ�ƺ��ٶ�
void Particle::updatePosition()
{
	int r1,r2;
	int i=0,j=0,k=0;

	for(i=0;i<GANVARS;i++)
	{	
		paramater_w = rnd(0,1);

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
	int gbest; //ȫ�ּ�ֵ����Ӧֵ
	int gbest_pos[GANVARS]; //ȫ�ּ�ֵ������
	Particle PSO_pop[PSO_popsize];//�������Ӷ���Ϊ����Ⱥ�������
public: 
	void init();  //��ʼ����Ⱥ
	void getGBest(); //��ȡȫ�ּ�ֵ
	void search(int *Array); //����,col������[0,RUNTIME-1]��ָʾ�ǵڼ������� 	  
};

//��ʼ����Ⱥ����
void ParticleSwarm::init(){

	gbest=numeric_limits<int>::max();//initialize gbest with maximal value of int type
	srand((unsigned)time(NULL));//get the seed of random number

	//initialize the boundary ****��Ҫ����*****�趨���ӱ߽�
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
		//x[] store the positions��y[] store the velocity 
		int x[GANVARS];
		int y[GANVARS];

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
void ParticleSwarm::search(int *Array)
{
	int gen=0;//number of iterations

	//������Χ��:[0, PSO_maxgen-1],�ܹ�PSO_maxgen���� 
	while(gen<PSO_maxgen)
	{
		Array[gen] = gbest;//Array���鴢��ȫ�ּ�ֵ�ĵ���Ӧֵ

		//ÿ�����ӽ����˶�����ֵ������pbest 
		for(int k=0;k<PSO_popsize;k++)
		{
			PSO_pop[k].updatePosition();
			PSO_pop[k].fitness=PSO_pop[k].calcFitness(PSO_pop[k].pos,gen);
			PSO_pop[k].updatePBest();
		}

		//update gbest.ע�������������ͨ�ŵ�ȫ�ֱ��� 
		getGBest();

		gen++;		

	}
}

