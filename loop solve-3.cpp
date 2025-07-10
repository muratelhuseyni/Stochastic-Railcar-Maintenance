#include "Common.h"
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>  
#include <thread>
//#include <array>
#include <vector>
#include <numeric> 
//#include <map>
#include <math.h>
#include <random>
#include <chrono>

//added for parallel_for
#include <windows.h>
#include <ppl.h>

using namespace concurrency;
using namespace std;
//added for parallel_for
//globalvars
//paralel for

const int tt = 28; //7 g�n - 4 time period/day
const int K = 1000;   //normalde 1000, senaryo say�s� //real

//const int kk = 50; //real
const int kk = 100; //real

//const int kk = 5;  //test
//const int K = 10;   //normalde 1000, senaryo say�s�     //test
//const int n = 2; // total number of jobs  //test

const int seedsayisi = 10;
const int perfkriteri = 5;

struct Job
{
	int ID;
	int rj;
	int Ymax;
	int dj;
	int age;

	Job(int IDin, int rjIn, int YmaxIn, int djIn, int ageIn)
	{
		ID = IDin;
		rj = rjIn;
		Ymax = YmaxIn;
		dj = djIn;
		age = ageIn;
	}
};

double* FindCI(int size, double* data, double alpha)
{
	double* CI = new double[2];

	double sum = 0.0;
	for (int i = 0; i < size; ++i)
		sum += data[i];

	double avgcost = sum / size;

	//cost sonu�lar� arraya gelmeli
	// https://www.programiz.com/cpp-programming/examples/standard-deviation
	//buradan sample var
	double sumssquaredif = 0.0;

	for (int i = 0; i < size; ++i)
		sumssquaredif += pow(data[i] - avgcost, 2);

	double stdev = sqrt(sumssquaredif / (size - 1));
	//double alpha = 0.05;
	double Z_alfaover2 = 1.96;

	double lb = avgcost - Z_alfaover2 * stdev / sqrt(size);
	double ub = avgcost + Z_alfaover2 * stdev / sqrt(size);

	CI[0] = lb;
	CI[1] = ub;

	return CI;
}


void GenerateJobInterval(vector<Job>& Jobs, int cdftarget, int n)
{
	double beta = 5; //shape //>1 olsun ki increasing failure rate olsun
	double alpha = 50; //scale - mtbf
	double u = 0.90;
	if (cdftarget == 0)
		u = 0.80;

	int targetage = ceil(alpha * pow(-log(1 - u), 1 / beta));
	int halflength = ceil(targetage * 0.10);

	for (int j = 0; j < n; j++)
	{
		int Ymax = 0;
		int Ymin = 0;
		int dj = 0;
		int middle = 0; //intervalde targeta denk gelir
		double rldcdf = 0;

		middle = targetage - Jobs[j].age;
		Jobs[j].rj = middle - halflength;
		Jobs[j].Ymax = middle + halflength;
		Jobs[j].dj = ceil(middle + 0.5 * halflength);

		//Jobs.emplace_back(j, Ymin, Ymax, dj, age, rldcdf);
	}
}

double SetupModel2(vector<vector<double>>& z2, int n, double* pi, int opercost, int Cadd, int pensla, int Cc, int Cp, int penymaxplus, int Yp, int Yc, int L, int* SLA, vector<Job>& Jobs)
{
	//Elhuseyni and ATU model
	IloEnv env2;
	IloModel model2(env2);

	IloNumVarArray E = CreateNumVarArray(env2, n, "E", 0, IloInfinity); //earliness
	IloNumVarArray T = CreateNumVarArray(env2, n, "T", 0, IloInfinity);//tardiness
	IloNumVarArray S = CreateNumVarArray(env2, n, "S", 0, IloInfinity); //Job start
	IloNumVarArray deltaYmaxplus = CreateNumVarArray(env2, n, "deltaYmax+", 0, IloInfinity); //Job start
	IloNumVarArray deltaSLA = CreateNumVarArray(env2, tt, "deltaSLA-", 0, IloInfinity); //sla_deltad_undersatisfaction

	IloNumVarArray gamma = CreateNumVarArray(env2, tt, "deltacap", 0, IloInfinity); //extra capacity

	BoolVarArray2 ZZ = CreateBoolVarArray2(env2, n, tt, "Z");   //maintenance begins for job j at time 1; else 0 -yukar�daki b_jt
	BoolVarArray2 etaa = CreateBoolVarArray2(env2, n, tt, "eta");   //operational variable
	IloBoolVarArray np = CreateBoolVarArray(env2, n, "n");   //operational variable

	//const1
	for (int j = 0; j < n; ++j) //3-job size
	{
		IloExpr exp(env2);

		for (int t = 0; t < tt; ++t)
			exp += ZZ[j][t];

		model2.add(exp + np[j] == 1);
		exp.end();
	}


	//capacity
	for (int t = 0; t < tt; ++t)
	{
		IloExpr exp(env2);
		for (int j = 0; j < n; ++j)
			for (int e = 0; e < Yp; ++e)
			{
				if (t - e < 0)
					break;
				exp += ZZ[j][t - e];
			}

		model2.add(exp <= L + gamma[t]);
		exp.end();
	}

	//availability
	for (int j = 0; j < n; ++j)
		for (int t = 0; t < tt; ++t)
		{
			IloExpr exp(env2);
			for (int e = 0; e < Yp; ++e)
			{
				if (t - e < 0)
					break;
				exp += ZZ[j][t - e];
			}

			model2.add(etaa[j][t] <= 1 - exp);
			exp.end();
		}

	//curtailment
	for (int t = 0; t < tt; ++t)
	{
		IloExpr exp(env2);
		for (int j = 0; j < n; ++j)
			exp += etaa[j][t];

		model2.add(exp + deltaSLA[t] >= SLA[t]);
		exp.end();
	}

	//earliness-tardiness
	for (Job& job : Jobs)
	{
		int j = job.ID;

		IloExpr exp(env2);

		for (int t = 0; t < tt; ++t)
			exp += t * ZZ[j][t];

		model2.add(S[j] == exp);
		exp.end();

		//model2.add(T[j] >= S[j] - (job.dj));

		model2.add(T[j] >= S[j] + tt * np[j] - job.dj);

		model2.add(E[j] >= job.dj - S[j] - tt * np[j]);

		model2.add(T[j] - deltaYmaxplus[j] <= job.Ymax - job.dj);
	}

	//violation upper bound
	/*for (int t = 0; t < tt; ++t)
	{
		model2.add(gamma[t] <= ceil(L * 0.50));
		model2.add(deltaSLA[t] <= ceil(0.50 * SLA[t]));
	}*/

	//write the objective
	IloExpr objExp2(env2);

	for (Job& job : Jobs)
	{
		int j = job.ID;
		objExp2 += T[j] + E[j] + penymaxplus * deltaYmaxplus[j];
	}

	for (int t = 0; t < tt; ++t)
		objExp2 += pensla * deltaSLA[t];

	for (int t = 0; t < tt; ++t)
		objExp2 += Cadd * gamma[t];

	for (Job& job : Jobs)
		for (int t = 0; t < tt; ++t)
			objExp2 += Cp * ZZ[job.ID][t];

	//oper cost
	for (int j = 0; j < n; ++j)
		for (int t = 0; t < tt; ++t)
			objExp2 += etaa[j][t] * opercost;

	IloObjective obj2 = IloMinimize(env2, objExp2, "obj");
	model2.add(obj2);
	objExp2.end();

	IloCplex cplex2(model2);
	cplex2.setParam(IloCplex::TiLim, 60 * 60);
	double time = cplex2.getCplexTime();
	IloBool success = cplex2.solve();

	double earl = 0;
	double tard = 0;
	double curtailcost = 0;
	double aditcost = 0;
	double opcost = 0;
	double prevtotal = 0;

	if (success && cplex2.isPrimalFeasible())
	{
		//Earl ve tard cost alal�m bunlar�n yerine

		for (int j = 0; j < n; ++j)
		{
			earl += cplex2.getValue(E[j]);
			tard += cplex2.getValue(T[j]);
		}

		for (int t = 0; t < tt; ++t)
			curtailcost += pensla * cplex2.getValue(deltaSLA[t]);

		for (int j = 0; j < n; ++j)
			for (int t = 0; t < tt; ++t)
			{
				z2[j][t] = cplex2.getValue(ZZ[j][t]);
				prevtotal += z2[j][t] * Cp;
			}

		//getadditcost
		for (int t = 0; t < tt; ++t)
			aditcost += Cadd * cplex2.getValue(gamma[t]);

		//getopercost

		//oper cost
		for (int j = 0; j < n; ++j)
			for (int t = 0; t < tt; ++t)
			{
				opcost += cplex2.getValue(etaa[j][t]) * opercost;
			}
	}

	double totalcost = earl + tard + curtailcost + aditcost + opcost + prevtotal;

	cplex2.end();
	model2.end();
	env2.end();

	return totalcost;
}


double SetupModelTest(vector<vector<double>>& z2, int n, int* tau, int* gamma, int* sla_v, double& costestcor, double& costestprev, int opercost, int Cadd, int pensla, int Cc, int Cp, int Yp, int Yc, int L, int* SLA)
{

	double objresult = 0;
	int* cornum = new int[tt];
	int* prevnum = new int[tt];

	for (int t = 0; t < tt; ++t)
	{
		cornum[t] = 0;
		prevnum[t] = 0;
	}

	for (int j = 0; j < n; ++j)
	{
		bool nomaint = true;
		for (int t = 0; t < tt; ++t)
		{
			if (z2[j][t] > 0.5)
			{
				if (t >= tau[j])
				{
					objresult += Cc;
					costestcor += Cc;
					nomaint = false;
					for (int cm = tau[j]; cm < min(tau[j] + Yc, tt); ++cm)
						cornum[cm] += 1;
				}
				else
				{
					objresult += Cp;
					costestprev += Cp;
					nomaint = false;
					for (int pm = t; pm < min(t + Yp, tt); ++pm)
						prevnum[pm] += 1;
				}
			}

		}

		if (nomaint && tt > tau[j])
		{
			objresult += Cc;
			costestcor += Cc;
			for (int cm = tau[j]; cm < min(tau[j] + Yc, tt); ++cm)
				cornum[cm] += 1;
		}

	}


	for (int t = 0; t < tt; ++t)
	{
		int maintjobs = prevnum[t] + cornum[t];
		int available = n - maintjobs;
		sla_v[t] = max(SLA[t] - available, 0);
		int exploit = min(available, SLA[t]);
		gamma[t] = max(maintjobs - L, 0);

		objresult += opercost * exploit + pensla * sla_v[t] + Cadd * gamma[t];
	}

	delete[] cornum;
	delete[] prevnum;
	return objresult;
}

void SetupModel1(int senaryo, int n, vector<vector<double>>& z1, double traincost[seedsayisi], double* pi, vector<vector<int>>& tau, int opercost, int Cadd, int pensla, int Cc, int Cp, int Yp, int Yc, int L, int* SLA)
{
	IloEnv env;
	IloModel model(env);

	BoolVarArray3 eta = CreateBoolVarArray3(env, n, tt, kk, "eta");   //operational recourse variable
	NumVarArray2 gamma = CreateNumVarArray2(env, tt, kk, "gamma");   //maintenance recourse variable (additional to the current capacity)
	BoolVarArray2 Z = CreateBoolVarArray2(env, n, tt, "Z");   //maintenance begins for job j at time 1; else 0 -yukar�daki b_jt
	NumVarArray2 deltaminusSLA = CreateNumVarArray2(env, tt, kk, "deltaSLA-", 0, IloInfinity); //sla_deltad_undersatisfaction
	BoolVarArray3 m = CreateBoolVarArray3(env, n, tt, kk, "m");
	IloBoolVarArray np = CreateBoolVarArray(env, n, "n");   //operational variable

	IloExpr objExp(env);

	for (int j = 0; j < n; ++j) //3-job size
	{
		IloExpr exp(env);

		for (int t = 0; t < tt; ++t)
			exp += Z[j][t];

		for (int k = 0; k < kk; ++k)
			if (tt > tau[j][k])
				model.add(exp + np[j] == 1);

		exp.end();
	}

	for (int j = 0; j < n; ++j) //3-job size
	{
		IloExpr exp(env);

		for (int t = 0; t < tt; ++t)
			exp += Z[j][t];

		for (int k = 0; k < kk; ++k)
			if (tau[j][k] >= tt)
				model.add(exp <= 1);

		exp.end();
	}

	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = 0; t < min(tau[j][k], tt); ++t)
			{

				IloExpr exp(env);
				for (int e = 0; e < Yp; ++e)
				{
					if (t - e < 0)
						break;
					exp += Z[j][t - e];
				}

				model.add(m[j][t][k] == exp);
				exp.end();
			}

	//24.4.25 - removal of it makes deterministic feasible!!
	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = tau[j][k]; t < min(tau[j][k] + Yc, tt); ++t)
			{
				IloExpr exp(env);
				int end = min(t - Yp + 1, tau[j][k]);
				for (int tprime = 0; tprime < end; ++tprime)
					exp += Z[j][tprime];
				model.add(m[j][t][k] == 1 - exp);
				exp.end();
			}

	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = tau[j][k] + Yc; t < tt; ++t)
				model.add(m[j][t][k] == 0);

	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = 0; t < tt; ++t)
				model.add(eta[j][t][k] <= 1 - m[j][t][k]);

	//extra capacity
	for (int k = 0; k < kk; ++k)
		for (int t = 0; t < tt; ++t)
		{
			IloExpr exp(env);
			for (int j = 0; j < n; ++j)
				exp += m[j][t][k];

			model.add(exp <= L + gamma[t][k]);
			exp.end();
		}

	//violation upper bound
	/*for (int k = 0; k < kk; ++k)
		for (int t = 0; t < tt; ++t)
		{
			model.add(gamma[t][k] <= ceil(L * 0.50));
			model.add(deltaminusSLA[t][k] <= ceil(0.50 * SLA[t]));
		}*/

	//demand curtailment
	for (int k = 0; k < kk; ++k)
		for (int t = 0; t < tt; ++t)
		{
			IloExpr exp(env);
			for (int j = 0; j < n; ++j)
				exp += eta[j][t][k];

			model.add(exp + deltaminusSLA[t][k] >= SLA[t]);
			//15.4.25
			exp.end();
		}

	//write the objective
	for (int k = 0; k < kk; ++k)
		for (int j = 0; j < n; ++j)
		{
			for (int t = 0; t < tau[j][k]; ++t)
				objExp += pi[k] * Cp * Z[j][t];
			for (int t = tau[j][k]; t < tt; ++t)
				objExp += Cc * pi[k] * Z[j][t];
			//28.4.25 new term
			if (tt - 1 >= tau[j][k])
				objExp += Cc * pi[k] * np[j];
		}

	for (int k = 0; k < kk; ++k)
		for (int t = 0; t < tt; ++t)
			objExp += pi[k] * pensla * deltaminusSLA[t][k];


	for (int k = 0; k < kk; ++k)
		for (int t = 0; t < tt; ++t)
			objExp += pi[k] * Cadd * gamma[t][k];


	//oper cost
	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = 0; t < tt; ++t)
				objExp += eta[j][t][k] * pi[k] * opercost;

	IloObjective obj = IloMinimize(env, objExp, "obj");
	model.add(obj);

	IloCplex cplex(model);
	cplex.setParam(IloCplex::TiLim, 60 * 60); //10 mins new instead of 60
	cplex.setParam(IloCplex::EpGap, 0.01);

	//cplex.exportModel("deneme.lp");
	IloBool success = cplex.solve();

	double earlcost = 0;
	double corcost = 0;
	double curtailcost = 0;
	double aditcost = 0;
	double opcost = 0;

	if (success && cplex.isPrimalFeasible())
	{
		for (int k = 0; k < kk; ++k)
			for (int j = 0; j < n; ++j)
			{
				for (int t = 0; t < tau[j][k]; ++t)
					earlcost += Cp * pi[k] * cplex.getValue(Z[j][t]);


				for (int t = tau[j][k]; t < tt; ++t)
					corcost += Cc * pi[k] * cplex.getValue(Z[j][t]);

				if (tt - 1 >= tau[j][k])
					corcost += Cc * pi[k] * cplex.getValue(np[j]);
			}

		for (int k = 0; k < kk; ++k)
			for (int t = 0; t < tt; ++t)
				curtailcost += pi[k] * pensla * cplex.getValue(deltaminusSLA[t][k]);
	}

	//getadditcost
	for (int k = 0; k < kk; ++k)
	{
		IloExpr objexp(env);
		for (int t = 0; t < tt; ++t)
			aditcost += pi[k] * Cadd * cplex.getValue(gamma[t][k]);
	}

	//getopercost

	//oper cost
	for (int j = 0; j < n; ++j)
		for (int k = 0; k < kk; ++k)
			for (int t = 0; t < tt; ++t)
			{
				double coef = pi[k] * opercost;
				opcost += cplex.getValue(eta[j][t][k]) * coef;
			}

	double totalcost = earlcost + corcost + curtailcost + aditcost + opcost;

	traincost[senaryo] = totalcost;

	//first stage decisions
	for (int j = 0; j < n; ++j)
		for (int t = 0; t < tt; ++t)
			z1[j][t] = cplex.getValue(Z[j][t]);

	cplex.end();
	model.end();
	env.end();

}



void GenerateJobs(vector<Job>& Jobs, int n, mt19937& rngmt, int agetype)
{
	int Ymax = 0;
	int Ymin = 0;
	int dj = 0;
	int middle = 0; //intervalde targeta denk gelir

	for (int j = 0; j < n; j++)
	{

		int age = 0;
		if (agetype == 0)
		{
			//uniform_int_distribution<int> ageunif(10, 20); //genc
			uniform_int_distribution<int> ageunif(30, 40); //genc
			age = ageunif(rngmt);
		}

		else if (agetype == 1)
		{
			//uniform_int_distribution<int> ageunif(50, 60); //ya�l�
			uniform_int_distribution<int> ageunif(70, 80); //ya�l�
			age = ageunif(rngmt);
		}

		else
		{
			//uniform_int_distribution<int> ageunif(10, 60); // kar���k
			uniform_int_distribution<int> ageunif(30, 80); // kar���k
			age = ageunif(rngmt);
		}




		//Jobs.emplace_back(j, age);
		Jobs.emplace_back(j, Ymin, Ymax, dj, age);
	}
}

void GenerateScenarios(vector<std::vector<int>>& scenK, vector<std::vector<std::vector<int>>>& tau, vector<Job>& Jobs, mt19937& rngmt, int type, int trial)
{
	double beta = 5; //shape //>1 olsun ki increasing failure rate olsun
	double alpha = 70; //scale - mtbf

	int scenario;
	if (type == 0)
		scenario = K;
	else
		scenario = kk;

	for (int k = 0; k < scenario; ++k)
		for (Job& job : Jobs)
		{
			double x;
			//int begweek = 70;
			uniform_real_distribution<double> zeroneunif(0.0, 1.0); //generate unirandom number
			double fx = zeroneunif(rngmt);

			x = ceil(pow(pow(job.age / alpha, beta) - log(1 - fx), 1 / beta) * alpha); //log burada ln, pow:power, ceil for discrete

			double additime = x - job.age;

			if (type == 0) //test
			{
				if (additime >= tt)
					scenK[job.ID][k] = tt;
				else
					scenK[job.ID][k] = additime;
			}

			else //train
			{
				if (additime >= tt)
					tau[trial][job.ID][k] = tt;
				else
					tau[trial][job.ID][k] = additime;
			}

		}

}

int main()
{
	int seed = 0; //20 �imdi
	//int initseed = 20; //20 �imdi

	int* SLAdaily;
	int* SLAsat;
	int* SLAsun;
	string pathparent;

	//const int kk = 50; //senaryo say�s�
	int day = 7;
	int period = 4;
	int realbegweek = 0;
	vector<Job> Jobs; //bu jobs seti tan�mlan�p doldurulmal�
	double pi[kk];

	for (int i = 0; i < kk; ++i)
		pi[i] = 1.00 / kk;

	int Yp = 3; //prev maint duration
	int Yc = 5; //cor maint duration

	//10.4.23 - main experimental settings

	int opercost = 2;
	//Preventive Job generation 
	int penymaxplus = opercost * 2;

	int Cc = opercost * 20; //failure penalty //yeni
	int Cp = Cc / 5;
	double pensla = 0;      //1.5, 1, 0.5

	int Cadd = Cp; //additional maintenance cost
	//double* testcost = new double[5];
	double model2cost;
	double alpha = 0.05;

	double epsilon = 0.10; // %90double epsilon = 0.04; // %96
	int rho = 15; //dney1 - 4 oldu, ba�ar�l�

	ofstream resultfile;

	for (int sla = 0; sla < 2; ++sla)
		for (int age = 0; age < 3; ++age)  //0:stcoh; 1:det
		{
			Jobs.clear();

			switch (age)
			{
				case 0:
					seed = 0;
					break;
				case 1:
					seed = 1;
					break;
				case 2:
					seed = 2;
					break;
			}

			int n = 57; // total number of jobs    //real
			if (sla == 1)
				n = 74;
			int L = ceil(n * 0.20); // track capacity //9.6.25

			mt19937 rngmt;
			rngmt.seed(seed); // to ensure the same solution when target changes 
			GenerateJobs(Jobs, n, rngmt, age);

			/*int scenK[n][K];
			int tau[seedsayisi][n][kk];*/

			std::vector<std::vector<int>> scenK(n, std::vector<int>(K, 0));

			// 3D vector for tau[seedsayisi][n][kk]
			std::vector<std::vector<std::vector<int>>> tau(
				seedsayisi, std::vector<std::vector<int>>(n, std::vector<int>(kk, 0)));

			GenerateScenarios(scenK, tau, Jobs, rngmt, 0, -1);

			for (int typ = 0; typ < 2; ++typ)
			{

				if (typ == 0)
					resultfile.open("stochastic_" + to_string(sla) + "_" + to_string(age) + "_expers.txt");
				else
					resultfile.open("detexpers" + to_string(sla) + "_" + to_string(age) + "_expers.txt");


				if (typ == 0)
					for (int i = 0; i < seedsayisi; ++i)
						GenerateScenarios(scenK, tau, Jobs, rngmt, 1, i); //train tau by i seed

				for (int cas = 0; cas < 3; ++cas)
				{
					if (cas == 0)
						pensla = Cc * 1.5;
					else if (cas == 1)
						pensla = Cc;
					else
						pensla = Cc * 2 / 3;

					for (int met = 0; met < 2; met++)
					{
						if (typ == 0 && met == 1)
							continue;

						double testdata[K][perfkriteri];
						double traincost[seedsayisi];
						//double model2testdata[K][perfkriteri];

						double totaltestime = 0;
						string method, casee, slatype, agetype;

						if (typ == 1)
						{
							GenerateJobInterval(Jobs, met, n); //modify det job interval
							if (met == 0)
								method = "Strict";
							else
								method = "Relaxed";
						}

						else
						{
							if (met == 0)
								method = "Stochastic";
							else
								method = "Stochastic+CC";
						}

						if (cas == 0)
							casee = "Case1";
						else if (cas == 1)
							casee = "Case2";
						else
							casee = "Case3";

						if (sla == 0)
							slatype = "M1B";
						else
							slatype = "M4";

						if (age == 0)
							agetype = "Young";
						else if (age == 1)
							agetype = "Old";
						else
							agetype = "Mixed";

						string filename = casee + "_" + slatype + "_" + agetype +
							"_" + method;

						resultfile << filename << "\n";

						//sla settings
						if (sla == 0) //m1
						{
							SLAdaily = new int[4] {0, 45, 50, 49};
							SLAsat = new int[4] {12, 41, 51, 47 };
							SLAsun = new int[4] {12, 40, 43, 41};
						}

						else //m4
						{
							SLAdaily = new int[4] {0, 68, 66, 59 };
							SLAsat = new int[4] {11, 55, 66, 56 };
							SLAsun = new int[4] {11, 48, 66, 54 };
						}


						//Generate SLA
						int* SLA = new int[tt];
						for (int i = 0; i < tt; ++i)
						{
							if (i < 20)
							{
								int k = i % 4;
								SLA[i] = SLAdaily[k];
							}

							else if (i < 24)
							{
								int k = i % 4;
								SLA[i] = SLAsat[k];
							}

							else
							{
								int k = i % 4;
								SLA[i] = SLAsun[k];
							}
						}


						//first stage decisions
						double totcost2 = 0;

						/*double z1_all[seedsayisi][n][tt];
						double z1best[n][tt];
						double z2[n][tt];*/

						std::vector<std::vector<std::vector<double>>> z1_all(
							seedsayisi, std::vector<std::vector<double>>(n, std::vector<double>(tt, 0.0))
						);

						std::vector<std::vector<double>> z1best(
							n, std::vector<double>(tt, 0.0)
						);

						std::vector<std::vector<double>> z2(
							n, std::vector<double>(tt, 0.0)
						);


						auto start = chrono::steady_clock::now();

						if (typ == 1)
							model2cost = SetupModel2(z2, n, pi, opercost, Cadd, pensla, Cc, Cp, penymaxplus, Yp, Yc, L, SLA, Jobs);
						else
							for (int i = 0; i < seedsayisi; ++i) //MODEL1 TRAIN 5 TIMES
							{
								//double z1[n][tt];
								std::vector<std::vector<double>> z1(n, std::vector<double>(tt, 0.0));

								SetupModel1(i, n, z1, traincost, pi, tau[i], opercost, Cadd, pensla, Cc, Cp, Yp, Yc, L, SLA);

								for (int job = 0; job < n; ++job)
									for (int t = 0; t < tt; ++t)
										z1_all[i][job][t] = z1[job][t];
							}

						auto end = chrono::steady_clock::now();
						auto diff = end - start;
						double tottraintime = chrono::duration_cast<chrono::milliseconds>(diff).count() / 1000;

						//Get the best performing total cost seed on average
						double min = 10000000000000;
						int index = -1;
						double avgtrain = 0;
						double* CI_train = new double[2];

						if (typ == 0)
						{
							for (int k = 0; k < seedsayisi; ++k)
							{
								avgtrain += traincost[k];
								if (min > traincost[k])
								{
									min = traincost[k];
									index = k;
								}

							}

							avgtrain = avgtrain / seedsayisi; //14.4.25
							CI_train = FindCI(seedsayisi, traincost, alpha);

							//Get the best performing train results
							for (int job = 0; job < n; ++job)
								for (int t = 0; t < tt; ++t)
									z1best[job][t] = z1_all[0][job][t];
						}
						//model2 Test
						auto start2 = chrono::steady_clock::now();

						//Test the best random seed
						parallel_for(0, K, [&](int k)
							{
								int* sla_v = new int[tt];
								int* gamma = new int[tt];

								for (int t = 0; t < tt; ++t)
								{
									sla_v[t] = 0;
									gamma[t] = 0;
								}

								double costestcor = 0;
								double costestprev = 0;

								int* testscenarios = new int[n];
								for (int j = 0; j < n; ++j)
									testscenarios[j] = scenK[j][k];

								double objresult;

								if (typ == 0)
									objresult = SetupModelTest(z1best, n, testscenarios, gamma, sla_v, costestcor, costestprev, opercost, Cadd, pensla, Cc, Cp, Yp, Yc, L, SLA);
								else
									objresult = SetupModelTest(z2, n, testscenarios, gamma, sla_v, costestcor, costestprev, opercost, Cadd, pensla, Cc, Cp, Yp, Yc, L, SLA);

								int slav = 0;
								int gammasum = 0;

								for (int t = 0; t < tt; ++t)
								{
									slav += sla_v[t];
									gammasum += gamma[t];
								}


								testdata[k][0] = costestprev / Cp;
								testdata[k][1] = costestcor / Cc;
								testdata[k][2] = slav / tt;
								testdata[k][3] = gammasum / tt;
								testdata[k][4] = objresult;

								delete[] testscenarios;
								delete[] sla_v;
								delete[] gamma;
							});

						auto end2 = chrono::steady_clock::now();
						auto diff2 = end2 - start2;
						totaltestime = chrono::duration_cast<chrono::seconds>(diff2).count();

						double avgcorjobs, avgprevjobs, avgslaviol, avgamma, avgtot;

						for (int j = 0; j < perfkriteri; ++j)
						{
							double avg = 0.0;
							for (int k = 0; k < K; ++k)
								avg += testdata[k][j];

							avg = avg / K;

							if (j == 0)
								avgprevjobs = avg;
							else if (j == 1)
								avgcorjobs = avg;
							else if (j == 2)
								avgslaviol = avg;
							else if (j == 3)
								avgamma = avg;
							else if (j == 4)
								avgtot = avg;
						}

						double* testcost = new double[K];
						for (int k = 0; k < K; ++k)
							testcost[k] = testdata[k][4];

						//model1 test interval
						double* CI_test = new double[2];
						CI_test = FindCI(K, testcost, alpha);

						double trueLB = min(CI_train[0], CI_test[0]);
						double trueUB = max(CI_train[1], CI_test[1]);
						double gapercent = (trueUB - trueLB) / trueUB * 100;

						if (typ == 0)
						{
							resultfile << "\n #prev\t" << "#cor \t" << "#sla_v \t" << "track_v\t" << "testcost\t" << "traincost\t" << "traintime\t" << "Testime\t"
								<< "CI for LB \t" << "CI for UB \t" << "Gap% \t" << "age gen seed \n";

							resultfile << avgprevjobs << "\t" << avgcorjobs << "\t" << avgslaviol << "\t" << avgamma << "\t" <<
								avgtot << "\t" << avgtrain << "\t" << tottraintime << "\t" << totaltestime << "\t"
								<< "[" << CI_train[0] << ";" << CI_train[1] << "]\t"
								<< "[" << CI_test[0] << ";" << CI_test[1] << "]\t" << gapercent << "\t" << seed << "\n\n";
						}

						else
						{
							resultfile << "\n #prev\t" << "#cor \t" << "#sla_v \t" << "track_v\t" << "Testcost\t" << "Modelcost\t" << "Modeltime\t" << "Testime\t"
								<< "CI for UB \t" << "Age gen seed \t" << "\n";

							resultfile << avgprevjobs << "\t" << avgcorjobs << "\t" << avgslaviol << "\t" << avgamma << "\t" <<
								avgtot << "\t" << model2cost << "\t" << tottraintime << "\t" << totaltestime << "\t"
								<< "[" << CI_test[0] << ";" << CI_test[1] << "]\t" << "\t" << seed << "\n\n";
						}


						//14.4.25
						delete[] SLA;
						delete[] testcost;
						delete[] CI_train;
						delete[] CI_test;
						//delete[] totraincost;
					}
				}

				resultfile.close();
			}

		}


}