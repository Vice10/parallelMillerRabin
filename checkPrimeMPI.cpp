#include <stdio.h>
#include "mpi.h"
#include <iostream>
#include <time.h>
#include <chrono>

using namespace std;
typedef unsigned long long uu;
int procNum;
int procRank;
int isPrime;

uu primeList[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
		,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
		227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,
		367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
		521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,
		677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,
		857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019 };

uu GCD(uu u, uu v)
{
	while (v != 0) {
		uu r = u % v;
		u = v;
		v = r;
	}
	return u;
}

uu naiveFactor(const uu& n)
{
	for (const uu& prime : primeList)
		if (GCD(n, prime) != 1)
			return n / prime;
	return n;
}

uu binpower(uu base, uu e, const uu& mod)
{
	uu result = 1;
	base %= mod;
	while (e > 0) {
		if (e % 2 != 0)
			result = (result * base) % mod;
		base = (base * base) % mod;
		e /= 2;
	}
	return result;
}

bool isComposite(const uu& n, const uu& a, const uu& d, const uu& s)
{
	uu x = binpower(a, d, n);
	if (x == 1 || x == n - 1)
		return false;
	for (uu r = 1; r < s; r++) {
		x = (x * x) % n;
		if (x == n - 1)
			return false;
	}
	return true;
}

bool MillerRabin(const uu& n, const uu& iter) // true if n is prime
{
	if (n < 4)
		return n == 2 || n == 3;
	if (naiveFactor(n) != n)
		return false;
	uu s = 0;
	uu d = n - 1;
	while (d % 2 == 0) {
		d /= 2;
		s += 1;
	}
	for (uu i = 0; i < iter; i++) {
		uu a = 2 + rand() % (n - 3);
		if (isComposite(n, a, d, s))
			return false;
	}
	return true;
}

void processInit(uu& n, uu& iter, uu& localCount) {
	uu* pProcIterCount = new uu[procNum];
	if (procRank == 0) {
		cout << "Enter a number:\n";
		cin >> n;
		do {
			cout << "Enter iter count:\n";
			cin >> iter;
			if (iter < procNum)
				cout << "Iter count is >= procNum. Try again:\n";
		} while (iter < procNum);
		uu procIterNum = iter / procNum;
		uu iterLeft = iter - (procIterNum * procNum);
		for (int i = 0; i < procNum; i++)
			pProcIterCount[i] = procIterNum;
		for (uu i = 0; i < iterLeft; i++)
			pProcIterCount[i] += 1;
	}
	MPI_Bcast(&n, sizeof(n), MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Scatter(pProcIterCount, 1, MPI_UNSIGNED_LONG_LONG, &localCount, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	delete[] pProcIterCount;
}

void doCheck(uu n, uu localIter) {
	bool probPrime = MillerRabin(n, localIter);
	probPrime ? isPrime = 1 : isPrime = 0;
}

int main(int* argc, char** argv) {
	uu n;
	uu iter;
	uu localCount;

	MPI_Init(argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	processInit(n, iter, localCount);
	double t1 = MPI_Wtime();
	doCheck(n, localCount);
	int* resArr = nullptr;
	if (procRank == 0) {
		resArr = new(nothrow)int[procNum];
		if (resArr == nullptr) exit(854);
	}
	MPI_Gather(&isPrime, 1, MPI_INT, resArr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double t2 = MPI_Wtime();
	double elapsed = t2 - t1;
	if (procRank == 0) {
		int sum = 0;
		for (uu i = 0; i < procNum; i++)
			sum += resArr[i];
		if (sum == procNum)
			cout << n << " is probably prime\n";
		else
			cout << n << " is definitely composite\n";
	}
	cout << "Execution time: " << elapsed << " s";
	MPI_Finalize();
	return 0;
}