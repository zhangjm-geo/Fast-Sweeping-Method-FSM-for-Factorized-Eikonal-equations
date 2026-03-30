/*
	2D isotropic media, FSM for factored eikonal equation(Multiplicative factor)
	Revised from 3D Factored Fast Sweeping Method Written By Zhang Jianming, 2022.10.10 (Factored_FSM.cpp)
	Dunshi Wu, 2026-3-29, OK
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

#define HUGEV 1e20   // huge value for eikonal
#define NX 101       // number of rows for model
#define NY 101       // number of column for model
#define H  10.0      // grid spacing

template <typename T>
T factored_fsm2(T **tau, T **Tm, T v, T v0, T T0, T T0x, T T0y, int i, int j, int orix, int oriy);
void orient_judge(int &ind_beg, int &ind_end, int ori, int min_edge, int max_edge);
template <typename T>
T factored_fsm1(T **tau, T **Tm, T s, T T0, T T0x, T T0y, int i, int j, int orix, int oriy, int xy);
template <typename T>
bool causality1(T ta, T **Tm, T T0, int i, int j, int orix, int oriy, int xy);
template <typename T>
bool causality2(T ta, T **Tm, T T0, int i, int j, int orix, int oriy);
bool circle_judge(int ori, int i, int i_end);

int main(void)
{
	// position of source
	int sx, sy;
	sx = NX/2; sy = NY/2;

	// allocate memory
	float **vel = new float*[NX];
	float **s   = new float*[NX];
	float **Tm  = new float*[NX];
	float **T0  = new float*[NX];
	float **tau = new float*[NX];

	for (int i = 0; i < NX; i++) {
		vel[i] = new float[NY];
		s[i]   = new float[NY];
		Tm[i]  = new float[NY];
		T0[i]  = new float[NY];
		tau[i] = new float[NY];
	}

	// velocity model
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			vel[i][j] = 1000.0;
			s[i][j] = 1.0 / vel[i][j];
		}
	}
	float v0 = vel[sx][sy];
	float s0 = s[sx][sy];

	// initialization of Tm, T0, tau
	// Tm -- the eikonal ultimately need
	// T0 -- the theoretical eikonal factor
	// tau -- the factor capturing the singularity of the source
	// Tm = T0*tau
	for (int i = 0; i < NX; i++) {
		float dx = (i - sx) * H;
		for (int j = 0; j < NY; j++) {
			float dy = (j - sy) * H;
			Tm[i][j] = HUGEV;
			float d = sqrt(dx*dx + dy*dy);
			T0[i][j] = s0 * d;
			tau[i][j] = HUGEV;
		}
	}
	Tm[sx][sy]  = 0.0f;
	tau[sx][sy] = 1.0f;
	cout << "sx = " << sx << ", sy = " << sy << endl;
	
	// Fast Sweeping Method
	clock_t t1 = clock();

	int maxit = 2;
	for (int iter = 1; iter <= maxit; iter++) {
		
		// 4 sweeping directions controlled by orix and oriy
		for (int orix = 1; orix > -2; orix -= 2) {
			for (int oriy = 1; oriy > -2; oriy -= 2) {

				int i_beg, i_end, j_beg, j_end;
				orient_judge(i_beg, i_end, orix, 0, NX);
				orient_judge(j_beg, j_end, oriy, 0, NY);

				int i, j;
				bool i_flag, j_flag;

				i = i_beg;
				i_flag = true;
				while (i_flag) {
					j = j_beg;
					j_flag = true;
					while (j_flag) {

						float tmp_tau;
						int ix, jy;
						ix = i - orix;
						if (ix >= NX) ix = NX - 1; if (ix < 0) ix = 0;
						jy = j - oriy;
						if (jy >= NY) jy = NY - 1; if (jy < 0) jy = 0;

                        /*** Very Important!!! They are dependent on orix and oriy***/
						float T0x = (T0[i][j] - T0[ix][j]) / H;
						float T0y = (T0[i][j] - T0[i][jy]) / H;
						/************************************************************/

						tmp_tau = factored_fsm2(tau, Tm, vel[i][j], v0, T0[i][j], T0x, T0y, i, j, orix, oriy);
						tau[i][j] = min(tmp_tau, tau[i][j]);
						Tm[i][j] = tau[i][j] * T0[i][j];

						j += oriy;
						j_flag = circle_judge(oriy, j, j_end);
					}

					i += orix;
					i_flag = circle_judge(orix, i, i_end);
				}
			}
		}
	}

	clock_t t2 = clock();
	double time_used = (double)(t2 - t1) / CLOCKS_PER_SEC;
	cout << "Execution Time: " << time_used << " s" << endl;

	// ===================== save eikonal Tm ===============
	ofstream ofile("Tm.dat", ios::binary);
	for (int i = 0; i < NX; i++) {
		ofile.write((char*)Tm[i], NY * sizeof(float));
	}
	ofile.close();

	// ===================== free memory ===================
	for (int i = 0; i < NX; i++) {
		delete[] vel[i];
		delete[] s[i];
		delete[] Tm[i];
		delete[] T0[i];
		delete[] tau[i];
	}
	delete[] vel;
	delete[] s;
	delete[] Tm;
	delete[] T0;
	delete[] tau;

	return 0;
}

template <typename T>
T factored_fsm2(T **tau, T **Tm, T v, T v0, T T0, T T0x, T T0y, int i, int j, int orix, int oriy)
// solve the equation (17) to get tau[i][j] on a triangle (Luo ,2012, page 365),except that this
// function deals with the isotropic media.
{
	int ix = i - orix;
	if (ix >= NX) ix = NX - 1; if (ix < 0) ix = 0;
	int jy = j - oriy;
	if (jy >= NY) jy = NY - 1; if (jy < 0) jy = 0;

	T tau_x = tau[ix][j];
	T tau_y = tau[i][jy];

	// 1 - find root according to equation (19)
	T roota = factored_fsm1(tau, Tm, 1.0f/v, T0, T0x, T0y, i, j, orix, oriy, 1);
	T rootb = factored_fsm1(tau, Tm, 1.0f/v, T0, T0x, T0y, i, j, orix, oriy, 2);
	T root = min(roota, rootb);

	// 2 - find roots according to equation (17)
	/*** coefficients of quadratic equation ***/
	T s0 = 1.0f/v0;
	T k1 = v*v*(s0*s0 + (2*T0 / H)*(T0x + T0y) + 2 * (T0 / H)*(T0 / H));
	T k2 = -2.0f * v*v*(T0 / H)*(T0x*tau_x + T0y*tau_y + (T0 / H)*(tau_x + tau_y));
	T k3 = v*v*(T0 / H)*(T0 / H)*(tau_x*tau_x + tau_y*tau_y);

	T A = k1;
	T B = k2;
	T C = k3 - 1;
	T D = B*B - 4*A*C;
	/******************************************/
	
	if (D >= 0.0) {
		D = sqrt(D);
		T root1 = (-B + D) / (2.0*A);
		T root2 = (-B - D) / (2.0*A);
		bool ca1 = causality2(root1, Tm, T0, i, j, orix, oriy);
		bool ca2 = causality2(root2, Tm, T0, i, j, orix, oriy);

		if (ca1 && ca2) root = min(root, min(root1, root2));
		else if (ca1)   root = min(root, root1);
		else if (ca2)   root = min(root, root2);
	}
	return root;
}

void orient_judge(int &ind_beg, int &ind_end, int ori, int min_edge, int max_edge)
{
	if (ori == 1) {
		ind_beg = min_edge;
		ind_end = max_edge - 1;
	}
	else if (ori == -1) {
		ind_beg = max_edge - 1;
		ind_end = min_edge;
	}
}

template <typename T>
T factored_fsm1(T **tau, T **Tm, T s, T T0, T T0x, T T0y, int i, int j, int orix, int oriy, int xy)
// find root according to equation (19)
{
	int ix = i - orix;
	if (ix >= NX) ix = NX - 1; if (ix < 0) ix = 0;
	int jy = j - oriy;
	if (jy >= NY) jy = NY - 1; if (jy < 0) jy = 0;

	T tau_x = tau[ix][j];
	T tau_y = tau[i][jy];
	T tau_a, T0a;

	if (xy == 1) {
		tau_a = tau_x;
		T0a = T0x;
	}
	else {
		tau_a = tau_y;
		T0a = T0y;
	}

	T root = (T0*tau_a + H * s) / (T0 + T0a * H);
	bool ca = causality1(root, Tm, T0, i, j, orix, oriy, xy);

	if (ca) return root;
	else return HUGEV;
}

template <typename T>
bool causality1(T ta, T **Tm, T T0, int i, int j, int orix, int oriy, int xy)
// causality condition for equation (19)
{
	int ix = i - orix;
	if (ix >= NX) ix = NX - 1; if (ix < 0) ix = 0;
	int jy = j - oriy;
	if (jy >= NY) jy = NY - 1; if (jy < 0) jy = 0;

	T px = ta * T0 - Tm[ix][j];
	T py = ta * T0 - Tm[i][jy];
	T pa;

	if (xy == 1) 
		pa = px;
	else 
		pa = py;

	return (pa >= 0.0);
}

template <typename T>
bool causality2(T ta, T **Tm, T T0, int i, int j, int orix, int oriy)
// causality condition for roots of the quadratic equation
{
	int ix = i - orix;
	if (ix >= NX) ix = NX - 1; if (ix < 0) ix = 0;
	int jy = j - oriy;
	if (jy >= NY) jy = NY - 1; if (jy < 0) jy = 0;

	T px = ta * T0 - Tm[ix][j];
	T py = ta * T0 - Tm[i][jy];
	return (px >= 0.0 && py >= 0.0);
}

bool circle_judge(int ori, int i, int i_end)
{
	if (ori == 1)      return i <= i_end;
	else if (ori == -1) return i >= i_end;
	return true;
}
