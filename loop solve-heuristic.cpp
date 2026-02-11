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
//paralel for


const int tt = 28; //7 g�n - 4 time period/day
const int K = 1000;
const int kk = 150; //real

//const int tt = 28; //7 g�n - 4 time period/day
//const int K = 1000;
//const int kk = 5; //real


//const int n = 2; // total number of jobs  //test

const int seedsayisi = 5;
const int perfkriteri = 5;

struct Job
{
	int ID;
	int Ymax;
	int dj;
	int age;
	int Sj;
	int Cj;

	Job(int IDin, int YmaxIn, int djIn, int ageIn, int SjIn, int CjIn)
	{
		ID = IDin;
		Sj = SjIn;
		Cj = CjIn;
		Ymax = YmaxIn;
		dj = djIn;
		age = ageIn;
	}
};

double* FindCI(int size, double* data, double train)
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

	double stdev = sqrt(sumssquaredif / (size*(size - 1)) );
	double table;

	if (train == 1)
	{
		if (size == 5)
			table = 2.77637; //t_0.05/2, 5-1
		else if (size == 10)
			table = 2.26211; //t_0.05/2, 10-1
	}
		
	else
		table = 1.96; //Z_0.05/2

	double lb = avgcost - table * stdev;
	double ub = avgcost + table * stdev;

	CI[0] = lb;
	CI[1] = ub;

	return CI;
}


void GenerateJobInterval(vector<Job>& Jobs, int cdftarget, int n)
{
	double beta = 5; //shape //>1 olsun ki increasing failure rate olsun
	double alpha = 50; //scale - mtbf
	//double alpha = 70; //scale - mtbf
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
		//Jobs[j].rj = middle - halflength;
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

	for (int j = 0; j < n; ++j)
		for (int t = 0; t < tt; ++t)
			model2.add(ZZ[j][t] == z2[j][t]);

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
	cplex2.setParam(IloCplex::TiLim, 60 * 60*3);

	double time = cplex2.getCplexTime();
	IloBool success = cplex2.solve();

	double earl = 0;
	double tard = 0;
	double curtailcost = 0;
	double aditcost = 0;
	double opcost = 0;
	double prevtotal = 0;
	double ymaxviolcost = 0; //1.02.26

	if (success && cplex2.isPrimalFeasible())
	{
		//Earl ve tard cost alal�m bunlar�n yerine

		for (int j = 0; j < n; ++j)
		{
			earl += cplex2.getValue(E[j]);
			tard += cplex2.getValue(T[j]);
			ymaxviolcost += penymaxplus*cplex2.getValue(deltaYmaxplus[j]);
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

	double totalcost = earl + tard + curtailcost + aditcost + opcost + prevtotal+ymaxviolcost; //1.2.26

	cplex2.end();
	model2.end();
	env2.end();

	return totalcost;
}

//double DispatchJobs(int* SLA,int capacity, vector<Job>& Jobs, int T, int Yp, vector<vector<double>>& z, int Cp
// int Cadd, int opercost, int pensla, int penymaxplus, int neigh)
//{
//	sort(Jobs.begin(), Jobs.end(),
//		[](const Job& a, const Job& b) {
//			return a.dj < b.dj;
//		});
//
//	vector<int> Hangar;   // store indices of Jobs
//	int n = Jobs.size();
//	int nextJob = 0;      // pointer to first unscheduled job
//
//	int totalcost = 0;
//
//	// for every t in tt
//	for (int t = 0; t < T; ++t)
//	{
//		int availablejobs = n - Hangar.size();
//		
//		// Remove finished jobs from Hangar
//		for (auto it = Hangar.begin(); it != Hangar.end(); )
//		{
//			Job& job = Jobs[*it];
//
//			if (t >= job.Sj && t < job.Cj)
//				++it; // still processing
//			else
//				it = Hangar.erase(it); // job finished
//		}
//
//		// Assign jobs if there is capacity
//		while (Hangar.size() < capacity && nextJob < n)
//		{
//			if (availablejobs < SLA[t])
//				break;
//
//			Job& job = Jobs[nextJob];
//
//			job.Sj = t;
//			job.Cj = t + Yp;
//
//			Hangar.push_back(nextJob);
//			nextJob++;
//			availablejobs--;
//
//			int earliness = max(0, job.dj - job.Sj);
//			int tardiness = max(0, job.Sj - job.dj);
//			int ymaxplus = max(0, job.Sj - job.Ymax);
//
//			totalcost += earliness;
//			totalcost += tardiness;
//			totalcost += ymaxplus * penymaxplus;
//			totalcost += Cp;
//		}
//
//
//		int slaviol = max(SLA[t] - availablejobs, 0);
//		int exploit = min(availablejobs, SLA[t]);
//		int addittrack = max(Hangar.size() - L, 0);
//
//		totalcost += opercost * exploit + pensla * slaviol + Cadd * addittrack;
//	}
//
//	//fill z
//	for (Job& job : Jobs)
//		for (int t = 0; t < tt; ++t)
//		{
//			z[job.ID][t] = 0;
//			if(job.Sj == t)
//				z[job.ID][t] = 1;
//		}
//}

double DispatchJobs(int* SLA, int capacity, vector<Job>& Jobs, int T, int Yp, vector<vector<double>>& z, int Cp, int Cadd, int opercost,
	int pensla, int penymaxplus, int neigh)
{
	// Sort jobs by EDD
	sort(Jobs.begin(), Jobs.end(),
		[](const Job& a, const Job& b) {
			return a.dj < b.dj;
		}); 

	vector<int> Hangar;   // indices of Jobs currently being processed
	int n = Jobs.size();
	int nextJob = 0;      // first not-yet-considered job in EDD order
	double totalcost = 0;

	/*double opcost = 0;
	double slaviolcost = 0;
	double addcost = 0;*/

	// For every time period
	for (int t = 0; t < T; ++t)
	{
		/*-----------------------------------
		  Remove finished jobs from Hangar
		-----------------------------------*/
		for (auto it = Hangar.begin(); it != Hangar.end(); )
		{
			Job& job = Jobs[*it];

			if (t >= job.Sj && t < job.Cj)
				++it;              // still processing
			else
				it = Hangar.erase(it); // finished
		}

		/*-----------------------------------
		  Build schedulable job set Jcons
		  (EDD-respecting + neigh constraint)
		-----------------------------------*/
		vector<int> Jcons;
		int availablejobs = n - Hangar.size();

		for (int j = nextJob; j < n; ++j)
		{
			Job& job = Jobs[j];
			int earliness = max(0, job.dj - t);

			if (earliness <= neigh)
				Jcons.push_back(j);
		}

		/*-----------------------------------
		  Assign jobs from Jcons
		-----------------------------------*/
		//while (Hangar.size() < capacity && !Jcons.empty() && availablejobs >= SLA[t]) //makalede raporlanan
		// while (!Jcons.empty() && availablejobs > SLA[t]) //berbat! en kötü!
		 while (Hangar.size() < capacity && !Jcons.empty() && availablejobs > SLA[t]) //makaledeki pseudo - daha iyi çıktı!
		{
			int jidx = Jcons.front();   // earliest feasible job
			Job& job = Jobs[jidx];

			job.Sj = t;
			job.Cj = t + Yp;

			Hangar.push_back(jidx);
			availablejobs--;			

			// Advance pointer only if we scheduled the front job
			if (jidx == nextJob)
				nextJob++;

			Jcons.erase(Jcons.begin());
		}

		/*-----------------------------------
		  SLA and operational costs
		-----------------------------------*/

		int slaviol = max(SLA[t] - availablejobs, 0);
		int exploit = min(availablejobs, SLA[t]);
		int addittrack = max((int)Hangar.size() - capacity, 0);

		totalcost += opercost * exploit
			+ pensla * slaviol
			+ Cadd * addittrack;

		/*opcost += opercost * exploit;
		slaviolcost += pensla * slaviol;
		addcost += Cadd * addittrack;*/
	}

	//cout << "opercost:" << opcost << endl;
	//cout << "slaviolcost:" << slaviolcost << endl;
	//cout << "addcost:" << addcost << endl;

	/*-----------------------------------
	  Fill z decision matrix
	-----------------------------------*/
	for (Job& job : Jobs)
	{
		bool nomaint = true;
		for (int t = 0; t < T; ++t)
		{
			z[job.ID][t] = 0;
			if (job.Sj == t && job.Cj > 0)
			{
				z[job.ID][t] = 1;
				nomaint = false;
				int earliness = max(0, job.dj - job.Sj);
				int tardiness = max(0, job.Sj - job.dj);
				int ymaxplus = max(0, tardiness - job.Ymax + job.dj);

				totalcost += earliness;
				totalcost += tardiness;
				totalcost += ymaxplus * penymaxplus;
				totalcost += Cp;
			}
		}

		if (nomaint)
		{
			int earliness = max(0, job.dj - tt );
			int tardiness = max(0, tt - job.dj);
			int ymaxplus = max(0, tardiness - job.Ymax + job.dj);
			totalcost += earliness;
			totalcost += tardiness;
			totalcost += ymaxplus * penymaxplus;
		}

	}

	return totalcost;
}

double SetupModelTest(vector<vector<double>>& z2, int n, int* tau, double* prevdec, double* cordec, double* oper, double* gamma, double* sla, double* cor, double* prev, int opercost, int Cadd, int pensla, int Cc, int Cp, int Yp, int Yc, int L, int* SLA)
{

	double objresult = 0;
	int* cornum = new int[tt];
	int* prevnum = new int[tt];
	int* gammanum = new int[tt];
	int* slanum = new int[tt];
	int* opernum = new int[tt];

	for (int t = 0; t < tt; ++t)
	{
		cornum[t] = 0;
		prevnum[t] = 0;
		gammanum[t] = 0;
		slanum[t] = 0;
		opernum[t] = 0;
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
					cordec[t] += 1;
					objresult += Cc;
					nomaint = false;
					for (int cm = tau[j]; cm < min(tau[j] + Yc, tt); ++cm)
						cornum[cm] += 1;
				}
				else
				{
					prevdec[t] += 1;
					objresult += Cp;
					nomaint = false;
					for (int pm = t; pm < min(t + Yp, tt); ++pm)
						prevnum[pm] += 1;
				}
			}

		}

		if (nomaint && tt > tau[j])
		{
			cordec[tau[j]] += 1;
			objresult += Cc;
			for (int cm = tau[j]; cm < min(tau[j] + Yc, tt); ++cm)
				cornum[cm] += 1;
		}

	}


	for (int t = 0; t < tt; ++t)
	{
		int maintjobs = prevnum[t] + cornum[t];
		int available = n - maintjobs;
		slanum[t] = max(SLA[t] - available, 0);
		int exploit = min(available, SLA[t]);
		opernum[t] = exploit;
		gammanum[t] = max(maintjobs - L, 0);

		objresult += opercost * exploit + pensla * slanum[t] + Cadd * gammanum[t];
	}


	for (int t = 0; t < tt; ++t)
	{
		gamma[t] += gammanum[t];
		sla[t] += slanum[t];
		cor[t] += cornum[t];
		prev[t] += prevnum[t];
		oper[t] += opernum[t];
	}


	delete[] cornum;
	delete[] prevnum;
	delete[] gammanum;
	delete[] slanum;
	delete[] opernum;

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
	cplex.setParam(IloCplex::TiLim, 60 * 60*3); //10 mins new instead of 60
	//cplex.setParam(IloCplex::EpGap, 0.01);

	//cplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersFull);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect, 0); //depth first search

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
			uniform_int_distribution<int> ageunif(10, 20); //genc
			//uniform_int_distribution<int> ageunif(30, 40); //genc
			age = ageunif(rngmt);
		}

		else if (agetype == 1)
		{
			uniform_int_distribution<int> ageunif(50, 60); //ya�l�
			//uniform_int_distribution<int> ageunif(70, 80); //ya�l�
			//uniform_int_distribution<int> ageunif(120, 150); //ya�l�
			age = ageunif(rngmt);
		}

		else
		{
			//uniform_int_distribution<int> ageunif(10, 80); // kar���k
			//uniform_int_distribution<int> ageunif(30, 80); // kar���k
			uniform_int_distribution<int> ageunif(10, 60); // kar���k
			age = ageunif(rngmt);
		}

		//Jobs.emplace_back(j, age);
		//Jobs.emplace_back(j, Ymin, Ymax, dj, age);
		Jobs.emplace_back(j, Ymax, dj, age, 0, 0);
	}
}

void GenerateScenarios(vector<std::vector<int>>& scenK, vector<std::vector<std::vector<int>>>& tau, vector<Job>& Jobs, mt19937& rngmt, int type, int trial)
{
	double beta = 5; //shape //>1 olsun ki increasing failure rate olsun
	double alpha = 50; //scale - mtbf //24.7.25
	//double alpha = 70; //scale - mtbf //24.7.25
	//double alpha = 70; //scale - mtbf //24.7.25

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
	int seed = 0; //initial seed

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
	ofstream resultdetail;
	
	for (int sla = 0; sla < 2; ++sla)
		for (int age = 0; age < 3; ++age)  //0:stcoh; 1:det
		{
			Jobs.clear();

			int n = 57; // total number of jobs    //real
			if (sla == 1)
				n = 74;
			int L = ceil(n * 0.20); // track capacity //9.6.25

			mt19937 rngmt;

			rngmt.seed(seed); // to ensure the same solution when target changes 

			GenerateJobs(Jobs, n, rngmt, age);

			std::vector<std::vector<int>> scenK(n, std::vector<int>(K, 0));

			// 3D vector for tau[seedsayisi][n][kk]
			std::vector<std::vector<std::vector<int>>> tau(
				seedsayisi, std::vector<std::vector<int>>(n, std::vector<int>(kk, 0)));

			GenerateScenarios(scenK, tau, Jobs, rngmt, 0, -1);

			for (int typ = 1; typ < 2; ++typ)
			{

				if (typ == 0)
				{
					resultdetail.open("stochastic_" + to_string(sla) + "_" + to_string(age) + "_detail.txt");
					resultfile.open("stochastic_" + to_string(sla) + "_" + to_string(age) + "_expers.txt");
				}

				else
				{
					resultdetail.open("detexpers" + to_string(sla) + "_" + to_string(age) + "_detail.txt");
					//resultfile.open("detexpers" + to_string(sla) + "_" + to_string(age) + "_expers.txt");
					resultfile.open("heurexpers" + to_string(sla) + "_" + to_string(age) + ".txt");
				}

				double* tautest = new double[tt];

				for (int t = 0; t < tt; ++t)
					tautest[t] = 0.0;

				for (int i = 0; i < tt; ++i)
				{
					for (int j = 0; j < n; ++j)
						for (int k = 0; k < K; ++k)
							if (i == scenK[j][k])
								tautest[i] += 1;
				}

				//for (int t = 0; t < tt; ++t)
					//resultdetail << "tau[" << to_string(t) + "]=" << tautest[t] / K << "\n";

				delete[] tautest;

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

						double* testcost = new double[K];
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
						resultdetail << filename << "\n";

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

						std::vector<std::vector<double>> z(n, std::vector<double>(tt, 0.0));

						auto start = chrono::steady_clock::now();
						int neigh = 100;
						double currcost = 1e18;   // large initial value
						
						currcost = DispatchJobs(SLA, L, Jobs, tt, Yp, z,Cp, Cadd, opercost, pensla, penymaxplus, neigh);
											
						auto end = chrono::steady_clock::now();
						auto diff = end - start;
						double tottraintime = chrono::duration_cast<chrono::milliseconds>(diff).count() / 1000;


						////Get the best performing total cost seed on average
						double min = 10000000000000;
						int index = -1;

						//model2 Test
						auto start2 = chrono::steady_clock::now();

						double* sla_v = new double[tt];
						double* cornum = new double[tt];
						double* prevnum = new double[tt];
						double* gammanum = new double[tt];
						double* opernum = new double[tt];

						//25.8.25
						double* prevdec = new double[tt];
						double* cordec = new double[tt];

						for (int t = 0; t < tt; ++t)
						{
							sla_v[t] = 0.0;
							cornum[t] = 0.0;
							prevnum[t] = 0.0;
							gammanum[t] = 0.0;
							opernum[t] = 0.0;
							prevdec[t] = 0.0;
							cordec[t] = 0.0;
						}

						//Test the best random seed
						for (int k = 0; k < K; ++k)
						{
							int* testscenarios = new int[n];
							for (int j = 0; j < n; ++j)
								testscenarios[j] = scenK[j][k];

							//double objresult;

							testcost[k] = SetupModelTest(z, n, testscenarios, prevdec, cordec, opernum, gammanum, sla_v, cornum, prevnum, opercost, Cadd, pensla, Cc, Cp, Yp, Yc, L, SLA);

							delete[] testscenarios;
						}

						auto end2 = chrono::steady_clock::now();
						auto diff2 = end2 - start2;
						totaltestime = chrono::duration_cast<chrono::seconds>(diff2).count();

						////model1 test interval
						////model1 test interval
						double* CI_test = new double[2];
						CI_test = FindCI(K, testcost, 0);

						/////
						for (int t = 0; t < tt; ++t)
						{
							opernum[t] = opernum[t] / K;
							resultdetail << "opernum[" << to_string(t) + "]=" << opernum[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							sla_v[t] = sla_v[t] / K;
							resultdetail << "slavnum[" << to_string(t) + "]=" << sla_v[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							cornum[t] = cornum[t] / K;
							resultdetail << "cornum[" << to_string(t) + "]=" << cornum[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							prevnum[t] = prevnum[t] / K;
							resultdetail << "prevnum[" << to_string(t) + "]=" << prevnum[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							gammanum[t] = prevnum[t] / K;
							resultdetail << "gammanum[" << to_string(t) + "]=" << gammanum[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							prevdec[t] = prevdec[t] / K;
							resultdetail << "prevdecision[" << to_string(t) + "]=" << prevdec[t] << "\n";
						}

						for (int t = 0; t < tt; ++t)
						{
							cordec[t] = cordec[t] / K;
							resultdetail << "cordecision[" << to_string(t) + "]=" << cordec[t] << "\n";
						}

						////present averages
						double prev = 0.0;
						double cor = 0.0;
						double oper = 0.0;
						double slaviol = 0.0;
						double trackviol = 0.0;
						double avgtestcost = 0.0;
						double prevbeg = 0.0;
						double corbeg = 0.0;

						for (int t = 0; t < tt; ++t)
						{
							oper += opernum[t];
							cor += cornum[t];
							prev += prevnum[t];
							slaviol += sla_v[t];
							trackviol += gammanum[t];
							prevbeg += prevdec[t];
							corbeg += cordec[t];
						}

						oper /= tt;
						cor /= tt;
						prev /= tt;
						slaviol /= tt;
						trackviol /= tt;
						prevbeg /= tt;
						corbeg /= tt;

						for (int k = 0; k < K; ++k)
							avgtestcost += testcost[k];

						avgtestcost /= K;

						resultfile << "\n #prev\t" << "#cor \t" << "#oper \t" << "#sla_v \t" << "track_v\t" << "Testcost\t" << "Modelcost\t" << "Modeltime\t" << "Testime\t"
							<< "CI for UB \t" << "#prevdec \t " << "#cordec \t " << "seed \n";

						resultfile << prev << "\t" << cor << "\t" << oper << "\t" << slaviol << "\t" << trackviol << "\t" <<
							avgtestcost << "\t" << currcost << "\t" << tottraintime << "\t" << totaltestime << "\t"
							<< "[" << CI_test[0] << ";" << CI_test[1] << "]\t" << prevbeg << "\t" << corbeg << "\t" << seed << "\n\n";

						////14.4.25
						delete[] SLA;
						delete[] prevdec;
						delete[] cordec;
						delete[] sla_v;
						delete[] cornum;
						delete[] prevnum;
						delete[] opernum;
						delete[] gammanum;
						delete[] testcost;
						//delete[] CI_train;
						delete[] CI_test;
						//delete[] totraincost;
					}
				}

				resultfile.close();
				resultdetail.close();
			}

			++seed;
		}
}
