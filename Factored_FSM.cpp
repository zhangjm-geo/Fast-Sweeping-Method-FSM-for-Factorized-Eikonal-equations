/*******************************************************************************/
/***************       3D Factored Fast Sweeping Method               **********/
/**************************Written By Zhang Jianming, 2022.10.10****************/
/***************************************CopyRight*******************************/

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>
#include<ctime>
using namespace std;
float pi=3.1415926;
int nx = 101;
int ny = 101;
int nz = 51;
float dx = 10;
float dy = 10;
float dz = 10;
float huge = 1e20;
float T_control(float a);
float min(float a,float b);
float max(float a,float b);
void BubbleSort(float* h, size_t len);
void Swap(float& a, float& b);
void FSM(float ***T,float***s,int i,int j,int k);
bool Circle_judge(int Ori,int i,int i_end);
void Orientation_judge(int &circle_begin,int &circle_end,int Ori,int MinEdge,int MaxEdge);
float Factored_FSM1(float ***tao,float ***T,float s,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz);
float Factored_FSM2(float ***tao,float ***T,float s0,float alp,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz);
float Factored_FSM3(float ***tao,float ***T,float s0,float alp,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ);
bool Causality1(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz);
bool Causality2(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz);
bool Causality3(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ);


int main()
{
	int i,j,k;
	//int x0 = nx/2;int y0=ny/2; int z0 = nz/2;
	int x0 = nx/2;int y0=ny/2; int z0 = 0;
	float ***s3=new float** [nx];
	float ***T_old3=new float** [nx];
	float ***T3=new float** [nx];
	float ***v3=new float** [nx];
	float ***tao=new float** [nx];
	float ***alp=new float** [nx];
	float ***T0=new float** [nx];
	for(int i=0;i<nx;i++)
	{
		s3[i]=new float *[ny];
		T_old3[i]=new float *[ny];
		T3[i]=new float *[ny];
		v3[i]=new float *[ny];
		tao[i]=new float *[ny];
		alp[i]=new float *[ny];
		T0[i]=new float *[ny];
		for(int j=0;j<ny;j++)
		{
			s3[i][j]=new float[nz];
			T_old3[i][j]=new float[nz];
			v3[i][j]=new float[nz];
			T3[i][j]=new float[nz];
			tao[i][j]=new float[nz];
			alp[i][j]=new float[nz];
			T0[i][j]=new float[nz];
		}
	}



	//initial for 3D
	for(i=0;i<nx;i++)   
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				v3[i][j][k]=1000;
				s3[i][j][k]=(1.0/v3[i][j][k]);
				T3[i][j][k]=huge;
				T_old3[i][j][k]=huge;
			}
	int sx=x0;int sy=y0;int sz=z0;
	T_old3[sx][sy][sz]=0.0;
	float s0=s3[sx][sy][sz];
	for(i=0;i<nx;i++)   
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				float d=sqrt(pow(i-sx,2)*dx*dx+pow(j-sy,2)*dy*dy+pow(k-sz,2)*dz*dz);
				alp[i][j][k]=s3[i][j][k]/s0;
				T0[i][j][k]=s0*d;
				tao[i][j][k]=huge;
			}
	tao[sx][sy][sz]=1.0;

	int loop=0,Max_loop=1;
	/***********FP FSM or Factored FSM for 3D*****************/
	float tt1=clock();
	for(loop=0;loop<Max_loop;loop++){
		int i_begin,i_end,j_begin,j_end,k_begin,k_end;
		bool i_flag,j_flag,k_flag;
		int OriX,OriY,OriZ;

		for(OriX=1;OriX>-2;OriX-=2){
			for(OriY=1;OriY>-2;OriY-=2){
				for(OriZ=1;OriZ>-2;OriZ-=2){
					//judge for FSM, need to change in this function
					Orientation_judge(i_begin , i_end , OriX , 0 , nx);
					Orientation_judge(j_begin , j_end , OriY , 0 , ny);
					Orientation_judge(k_begin , k_end , OriZ , 0 , nz);

					//judge for improved FSM, need to change in this function
					//Orientation_judge(i_begin , i_end , OriX , sx , nx);
					//Orientation_judge(j_begin , j_end , OriY , sy , ny);
					//Orientation_judge(k_begin , k_end , OriZ , sz , nz);

					i = i_begin;
					i_flag = true;
					while(i_flag){
						j = j_begin;
						j_flag = true;
						while(j_flag){
							k = k_begin;
							k_flag = true;
							while(k_flag){
								////FSM
								//FSM(T_old3,s3,i,j,k);

								//Factored FSM
								float tmp_tao;
								int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
								int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
								int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;
								float T0x=(T0[i][j][k]-T0[xx][j][k])/dx;
								float T0y=(T0[i][j][k]-T0[i][yy][k])/dy;
								float T0z=(T0[i][j][k]-T0[i][j][zz])/dz;

								tmp_tao=Factored_FSM3(tao,T_old3,s0,alp[i][j][k],T0[i][j][k],T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ);
								tao[i][j][k]=min(tao[i][j][k],tmp_tao);
								T_old3[i][j][k]=tao[i][j][k]*T0[i][j][k];

								k += OriZ;
								k_flag = Circle_judge(OriZ,k,k_end); 
							}
							j += OriY;
							j_flag = Circle_judge(OriY,j,j_end); 
						}
						i += OriX;
						i_flag = Circle_judge(OriX,i,i_end); 
					}
				}
			}
		}
	}


	float tt2=clock();
	cout<<"One shot modeling for 3DFSM is "<<(double)(tt2-tt1)/CLOCKS_PER_SEC<<" s."<<endl;
	float *Rdep=new float[nx];


	ofstream T3out("T3.dat");
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int k=0;k<nz;k++)
			{
				T_old3[i][j][k]-=sqrt(pow(i-x0,2)*dx*dx+pow(j-y0,2)*dy*dy+pow(k-z0,2)*dz*dz)*s0;
				T3out.write((char*)&T_old3[i][j][k],sizeof(float));
			}
		}
	}
	cout<<T_old3[nx-1][ny-1][nz-1]<<endl;;
	T3out.close();

	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			delete [] s3[i][j];
			delete [] T_old3[i][j];
			delete [] v3[i][j];
			delete [] T3[i][j];
			delete [] tao[i][j];
			delete [] alp[i][j];
			delete [] T0[i][j];
		}
		delete []s3[i];
		delete [] T_old3[i];
		delete []v3[i];
		delete [] T3[i];
		delete [] tao[i];
		delete [] alp[i];
		delete [] T0[i];
	}
	delete []s3;
	delete [] T_old3;
	delete []v3;
	delete [] T3;
	delete [] tao;
	delete [] alp;
	delete [] T0;
	return 0;
}//end of main




float T_control(float a)
{
	if(a>0) return  a;
	else return 0;

}
float min(float a,float b)
{
	if (a>b) return b;
	else return a;
}
float max(float a,float b)
{
	if (a>b) return a;
	else return b;
}
void BubbleSort(float* h, size_t len)
{
	if(h==NULL) return;
	if(len<=1) return;
	//i是次数，j是具体下标 
	for(int i=0;i<len-1;++i)
		for(int j=0;j<len-1-i;++j)
			if(h[j]>h[j+1])
				Swap(h[j],h[j+1]);               
	return;
} 
void Swap(float& a, float& b)
{
	float t=a;
	a=b;
	b=t;
	return;
}

float Factored_FSM1(float ***tao,float ***T,float s,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;

	float tao_xx=tao[xx][j][k];
	float tao_yy=tao[i][yy][k];
	float tao_zz=tao[i][j][zz];

	float da,tao_aa,T0a;
	if(xyz==1){da=dx;tao_aa=tao_xx;T0a=T0x;}
	if(xyz==2){da=dy;tao_aa=tao_yy;T0a=T0y;}
	if(xyz==3){da=dz;tao_aa=tao_zz;T0a=T0z;}
	float root=(s*da+T0*tao_aa)/(T0+T0a*da);
	if(Causality1(root,T0,T,i,j,k,OriX,OriY,OriZ,xyz)){return root;}
	else return huge;
}

float Factored_FSM2(float ***tao,float ***T,float s0,float alp,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;

	float tao_xx=tao[xx][j][k];
	float tao_yy=tao[i][yy][k];
	float tao_zz=tao[i][j][zz];
	float da,db,tao_aa,tao_bb,T0a,T0b;
	float Ta,Tb;
	if(xyz==12){da=dx;db=dy;tao_aa=tao_xx;tao_bb=tao_yy;T0a=T0x;T0b=T0y;Ta=T[xx][j][k];Tb=T[i][yy][k];}
	if(xyz==13){da=dx;db=dz;tao_aa=tao_xx;tao_bb=tao_zz;T0a=T0x;T0b=T0z;Ta=T[xx][j][k];Tb=T[i][j][zz];}
	if(xyz==23){da=dy;db=dz;tao_aa=tao_yy;tao_bb=tao_zz;T0a=T0y;T0b=T0z;Ta=T[i][yy][k];Tb=T[i][j][zz];}

	float root=huge,roota=huge,rootb=huge;
	//if(T[i][j][k]>Ta)roota=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,xyz/10);
	//if(T[i][j][k]>Ta)rootb=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,xyz%10);
	roota=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,xyz/10);
	rootb=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,xyz%10);
	root=min(roota,rootb);

	//if(T[i][j][k]>Ta&&T[i][j][k]>Tb){
	float A=T0*T0*db*db;
	float B=T0*T0*da*da;
	float A1=2.0*T0*T0a*da*db*db;
	float B1=2.0*T0*T0b*db*da*da;
	float A2=A+B+A1+B1+s0*s0*da*da*db*db;
	float B2=-((2.0*A+A1)*tao_aa+(2.0*B+B1)*tao_bb);
	float C2=A*tao_aa*tao_aa+B*tao_bb*tao_bb-s0*s0*alp*alp*da*da*db*db;
	float D2=sqrt(B2*B2-4.0*A2*C2);
	if(D2>=0)
	{
		float root1=(-B2+D2)/2.0/A2;
		float root2=(-B2-D2)/2.0/A2;
		if(Causality2(root1,T0,T,i,j,k,OriX,OriY,OriZ,xyz)&&Causality2(root2,T0,T,i,j,k,OriX,OriY,OriZ,xyz)){root=min(root,min(root1,root2));}
		else if(Causality2(root1,T0,T,i,j,k,OriX,OriY,OriZ,xyz)&&!Causality2(root2,T0,T,i,j,k,OriX,OriY,OriZ,xyz)) {root=min(root,root1);}
		else if(Causality2(root2,T0,T,i,j,k,OriX,OriY,OriZ,xyz)&&!Causality2(root1,T0,T,i,j,k,OriX,OriY,OriZ,xyz)) {root=min(root,root2);}
	}
	//}
	return root;
}


float Factored_FSM3(float***tao,float ***T,float s0,float alp,float T0,float T0x,float T0y,float T0z,int i,int j,int k,int OriX,int OriY,int OriZ){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;

	float tmp_tao=huge;
	float tao_xx=tao[xx][j][k];
	float tao_yy=tao[i][yy][k];
	float tao_zz=tao[i][j][zz];
	float root=huge;float roota=huge;float rootb=huge;float rootc=huge;
	//if(T[i][j][k]>T[xx][j][k])roota=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,1);
	//if(T[i][j][k]>T[i][yy][k])rootb=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,2);
	//if(T[i][j][k]>T[i][j][zz])rootc=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,3);
	roota=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,1);
	rootb=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,2);
	rootc=Factored_FSM1(tao,T,s0,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,3);
	root=min(roota,min(rootb,rootc));

	roota=huge;rootb=huge;rootc=huge;
	//if(T[i][j][k]>T[xx][j][k]&&T[i][j][k]>T[i][yy][k])roota=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,12);
	//if(T[i][j][k]>T[xx][j][k]&&T[i][j][k]>T[i][j][zz])rootb=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,13);
	//if(T[i][j][k]>T[i][yy][k]&&T[i][j][k]>T[i][j][zz])rootc=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,23);
	roota=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,12);
	rootb=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,13);
	rootc=Factored_FSM2(tao,T,s0,alp,T0,T0x,T0y,T0z,i,j,k,OriX,OriY,OriZ,23);
	root=min(roota,min(rootb,rootc));

	//if(T[i][j][k]>T[xx][j][k]&&T[i][j][k]>T[i][yy][k]&&T[i][j][k]>T[i][j][zz]){
	float A=T0*T0*dy*dy*dz*dz;
	float B=T0*T0*dx*dx*dz*dz;
	float C=T0*T0*dx*dx*dy*dy;
	float A1=2.0*T0*T0x*dx*dy*dy*dz*dz;
	float B1=2.0*T0*T0y*dy*dx*dx*dz*dz;
	float C1=2.0*T0*T0z*dz*dx*dx*dy*dy;
	float A2=A+B+C+A1+B1+C1+s0*s0*dx*dx*dy*dy*dz*dz;
	float B2=-((2.0*A+A1)*tao_xx+(2.0*B+B1)*tao_yy+(2.0*C+C1)*tao_zz);
	float C2=A*tao_xx*tao_xx+B*tao_yy*tao_yy+C*tao_zz*tao_zz-s0*s0*alp*alp*dx*dx*dy*dy*dz*dz;
	float D2=sqrt(B2*B2-4.0*A2*C2);
	if(D2>=0)
	{
		float root1=(-B2+D2)/2.0/A2;
		float root2=(-B2-D2)/2.0/A2;
		if(Causality3(root1,T0,T,i,j,k,OriX,OriY,OriZ)&&Causality3(root2,T0,T,i,j,k,OriX,OriY,OriZ)){root=min(root,min(root1,root2));}
		else if(Causality3(root1,T0,T,i,j,k,OriX,OriY,OriZ)&&!Causality3(root2,T0,T,i,j,k,OriX,OriY,OriZ)) {root=min(root,root1);}
		else if(Causality3(root2,T0,T,i,j,k,OriX,OriY,OriZ)&&!Causality3(root1,T0,T,i,j,k,OriX,OriY,OriZ)) {root=min(root,root2);}
	}
	//}
	return root;
}

bool Causality1(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;
	float px=tao*T0-T[xx][j][k];
	float py=tao*T0-T[i][yy][k];
	float pz=tao*T0-T[i][j][zz];
	float pa;
	if(xyz==1) pa=px;
	if(xyz==2) pa=py;
	if(xyz==3) pa=pz;

	if(pa>0){return true;}
	else{return false;}
}
bool Causality2(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ,int xyz){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;
	float px=tao*T0-T[xx][j][k];
	float py=tao*T0-T[i][yy][k];
	float pz=tao*T0-T[i][j][zz];
	float pa,pb;
	if(xyz==12) {pa=px;pb=py;}
	if(xyz==13) {pa=px;pb=pz;}
	if(xyz==23) {pa=py;pb=pz;}
	if(pa>0&&pb>0){return true;}
	else return false;
}
bool Causality3(float tao,float T0,float ***T,int i,int j,int k,int OriX,int OriY,int OriZ){
	int xx=i-OriX;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j-OriY;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k-OriZ;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;
	float px=tao*T0-T[xx][j][k];
	float py=tao*T0-T[i][yy][k];
	float pz=tao*T0-T[i][j][zz];
	if(px>0&&py>0&&pz>0){return true;}
	else return false;
}

void FSM(float ***T,float***s,int i,int j,int k){
	int xa=i+1;if(xa==nx) xa=nx-1;
	int xb=i-1;if(xb==-1) xb=0;
	int ya=j+1;if(ya==ny) ya=ny-1;
	int yb=j-1;if(yb==-1) yb=0;
	int za=k+1;if(za==nz) za=nz-1;
	int zb=k-1;if(zb==-1) zb=0;

	float T_xmin=33.0;float T_ymin=33.0; float T_zmin=33.0;float t_tmp;
	int sx,sy,sz;
	T_xmin=min(T[xa][j][k],T[xb][j][k]);if(T_xmin==T[xa][j][k])sx=1;else sx=-1;
	T_ymin=min(T[i][ya][k],T[i][yb][k]);if(T_ymin==T[i][ya][k])sy=1;else sy=-1;
	T_zmin=min(T[i][j][za],T[i][j][zb]);if(T_zmin==T[i][j][za])sz=1;else sz=-1;

	int xx=i+sx;if(xx==nx) xx=nx-1;if(xx==-1)xx=0;
	int yy=j+sy;if(yy==ny) yy=ny-1;if(yy==-1)yy=0;
	int zz=k+sz;if(zz==nz) zz=nz-1;if(zz==-1)zz=0;
	//dx!=dy!=dz
	float A=dx*dx*dy*dy+dy*dy*dz*dz+dz*dz*dx*dx;
	float B=-2.0*(dy*dy*dz*dz*T_xmin+dx*dx*dy*dy*T_zmin+dz*dz*dx*dx*T_ymin);
	float C=-pow(s[i][j][k],2)*dx*dx*dy*dy*dz*dz+T_xmin*T_xmin*dy*dy*dz*dz+T_ymin*T_ymin*dz*dz*dx*dx+T_zmin*T_zmin*dx*dx*dy*dy;
	float D=sqrt(B*B-4.0*A*C);
	float A1=dx*dx+dz*dz;
	float B1=-2.0*(dz*dz*T_xmin+dx*dx*T_zmin);
	float C1=-pow(s[i][j][k],2)*dx*dx*dz*dz+T_xmin*T_xmin*dz*dz+T_zmin*T_zmin*dx*dx;
	float D1=sqrt(B1*B1-4.0*A1*C1);
	float A2=dy*dy+dz*dz;
	float B2=-2.0*(dz*dz*T_ymin+dy*dy*T_zmin);
	float C2=-pow(s[i][j][k],2)*dy*dy*dz*dz+T_ymin*T_ymin*dz*dz+T_zmin*T_zmin*dy*dy;
	float D2=sqrt(B2*B2-4.0*A2*C2);
	float A3=dx*dx+dy*dy;
	float B3=-2.0*(dy*dy*T_xmin+dx*dx*T_ymin);
	float C3=-pow(s[i][j][k],2)*dx*dx*dy*dy+T_xmin*T_xmin*dy*dy+T_ymin*T_ymin*dx*dx;
	float D3=sqrt(B3*B3-4.0*A3*C3);
	T[i][j][k]=min(T[i][j][k],T_xmin+s[i][j][k]*dx);
	T[i][j][k]=min(T[i][j][k],T_ymin+s[i][j][k]*dy);
	T[i][j][k]=min(T[i][j][k],T_zmin+s[i][j][k]*dz);
	if(D1>=0&&T[i][j][k]>T_xmin&&T[i][j][k]>T_zmin&&T[i][j][k]<=T_ymin)
	{
		t_tmp=(-B1+D1)/2.0/A1;
		T[i][j][k]=min(T[i][j][k],t_tmp);
	}
	else if(D2>=0&&T[i][j][k]>T_ymin&&T[i][j][k]>T_zmin&&T[i][j][k]<=T_xmin)
	{
		t_tmp=(-B2+D2)/2.0/A2;
		T[i][j][k]=min(T[i][j][k],t_tmp);
	}
	else if(D3>=0&&T[i][j][k]>T_ymin&&T[i][j][k]>T_xmin&&T[i][j][k]<=T_zmin)
	{
		t_tmp=(-B3+D3)/2.0/A3;
		T[i][j][k]=min(T[i][j][k],t_tmp);
	}
	else if(D>=0&&T[i][j][k]>T_ymin&&T[i][j][k]>T_xmin&&T[i][j][k]>T_zmin)
	{
		t_tmp=(-B+D)/2.0/A;
		T[i][j][k]=min(T[i][j][k],t_tmp);
	}
	//SPAFSM
	T[i][j][k]=min(T[i][j][k],T[xx][j][zz]+sqrt(dx*dx+dz*dz)*s[xx][j][zz]);
	T[i][j][k]=min(T[i][j][k],T[i][yy][zz]+sqrt(dy*dy+dz*dz)*s[i][yy][zz]);
	T[i][j][k]=min(T[i][j][k],T[xx][yy][k]+sqrt(dy*dy+dx*dx)*s[xx][yy][k]);
	T[i][j][k]=min(T[i][j][k],T[xx][yy][zz]+sqrt(dy*dy+dx*dx+dz*dz)*s[xx][yy][zz]);

	//dx=dy=dz
	//float a[3]={T_xmin,T_ymin,T_zmin};
	//BubbleSort(a, 3);
	//float W=a[0];float V=a[1];float U=a[2];
	//T3[i][j][k]=W+s[i][j][k]*dx;
	//if(T3[i][j][k]>V) 
	//{
	//	T3[i][j][k]=((V+W+sqrt(2.0*pow(s[i][j][k]*dx,2)-pow(V-W,2))))/2.0;
	//	//T3[i][j][k]=((V+W+sqrt(-pow(V,2)-pow(W,2)+2.0*V*W))/(2.0*s[i][j][k]*dx*s[i][j][k]*dx))/2.0;
	//}
	//if(T3[i][j][k]>U)
	//{
	//	T3[i][j][k]=(V+W+U+sqrt(3.0*pow(s[i][j][k]*dx,2)-pow(V-W,2)-pow(V-U,2)-pow(U-W,2)))/3.0;
	//	float tssp=T3[xb][yb][zb]+s[xb][yb][zb]*dx*sqrt(3);
	//	T3[i][j][k]=min(T3[i][j][k],tssp);
	//	//float b=-2.0*(W+V+U);
	//	//T3[i][j][k]=(2.0*(W+V+U)+sqrt(b*b-12*(W*W+V*V+U*U-pow(s[i][j][k]*dx,2))))/6.0;
	//	//T3[i][j][k]=((2.0*pow(V+W+U,2)+sqrt(4.0*pow(V+W+U,2)))/(12.0*pow(W*W+V*V+U*U-s[i][j][k]*dx,2)))/2.0;
	//}
	//T[i][j][k]=min(T[i][j][k],T3[i][j][k]);
}
bool Circle_judge(int Ori,int i,int i_end){
	if(Ori == 1){
		if(i > i_end){
			return false;
		}
	} 
	else if(Ori == -1){
		if(i < i_end){
			return false;
		}
	}
	return true;
}
void Orientation_judge(int &circle_begin,int &circle_end,int Ori,int MinEdge,int MaxEdge){
	if(Ori == 1){
		//begin and end for FSM
		circle_begin = MinEdge;
		circle_end = MaxEdge - 1;
	}
	else{
		//begin and end for FSM
		circle_begin = MaxEdge - 1;
		circle_end = MinEdge;

		////begin and end for improved FSM
		//circle_begin = MinEdge;
		//if(MaxEdge==nz) circle_begin = nz-1;//to ensure sweeping from k>maxelev
		//circle_end = 0;
	}
}
