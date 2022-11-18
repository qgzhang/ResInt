#include "mex.h"
#include "math.h"
#include "time.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "float.h"
//#include "c_finddir.h"
//#define sign(x) (0 < x) - (x < 0)

/* Triangulation: estimates 3D coordinates of a scent point.
   Return:
        	X: 3D coordinates
        	mres: maximum reprojection error in L2 norm
        	nitr: number of internal iterations
   Input:
   			a1, a2, b1, b2, c, d: cameras setting
   			X0: initial coordinates
*/

int debug = 0;

double Inf = DBL_MAX;

double* getneggrad(double* A1, double* A2, double* B1, double* B2, double* C, double* D, \
	      int M, double* X, int nActive, int Active[nActive], double* grads[4] );

void getbound(double* C, double* D, int M, double* X, double* dir, double* lb0, double* ub0);

double stepsize(double* A1, double* A2, double* B1, double* B2, double* C, double* D, \
	            int M, double* X, double* dir, double lb0, double ub0);

void tri_pMat1(int m, int n, double A[m][n]);
void tri_pMat2(int m, int n, double* A);
void tri_pMatgrad(int m, int n, double* grads);
void tri_pVec(int m, double b[m]);
void tri_pVecInt(int m, int b[m]);

void meb(double* grads[4], int N, double* dir);


/***************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments. */    
    if (nrhs != 8)
    {
        mexErrMsgTxt("8 inputs required.");
    }
    else if (nlhs > 3)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
  
    if (debug) printf("just checked params\n"); fflush(stdout);

  	double *A1, *A2, *B1, *B2, *C, *D, *x0, *eps;
  	int M;
    /* Assign pointers to inputs.*/
    A1 = mxGetPr(prhs[0]);
    A2 = mxGetPr(prhs[1]);
    B1 = mxGetPr(prhs[2]);
    B2 = mxGetPr(prhs[3]);
     C = mxGetPr(prhs[4]);
     D = mxGetPr(prhs[5]);
    x0 = mxGetPr(prhs[6]);
    eps= mxGetPr(prhs[7]);
  
    M = mxGetM(prhs[0]);

	if (debug) printf("just got params\n"); fflush(stdout);


    /* Create matrices for the return argument. */    
    plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
	plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    /* Assign pointers to output.*/
    /*Ps = mxGetPr(plhs[0]);*/
    double* X_return    = mxGetPr(plhs[0]);
    double* mres_return = mxGetPr(plhs[1]);
    double* nitr_return = mxGetPr(plhs[2]);


	double X0[3] = {x0[0],x0[1],x0[2]};
	double X[3];
	int Active[4], nActive;
	
	double g1[3], g2[3], g3[3], g4[3];
	double* grads[4] = {g1,g2,g3,g4};

	double tmpforgrad;
	double r; // min dot product
	double norm;
	double lb0, ub0, *bound;
	double alpha;

	double mres0 = Inf;
	double mres = Inf/2;
	int nitr = 0;
	bool flag = 1;
	double tol = 0.00000001; //1e-8
	double a1[3], a2[3], b1, b2, c[3], d;
	double res[M], maxres;
	double dir[3];
	int i;
/*
    printf("\nA1 = ");  tri_pMat2(M,3,A1);
    printf("\nA2 = ");   tri_pMat2(M,3,A2);
    printf("\nC = ");   tri_pMat2(M,3,C);
    printf("\nB1 = ");    tri_pVec(M,B1);    
    printf("\nB2 = ");    tri_pVec(M,B2);    
    printf("\nD = ");    tri_pVec(M,D);    
    printf("\nX0 = ");    tri_pVec(3,X0); 
*/
	X[0] = X0[0]; X[1] = X0[1]; X[2] = X0[2]; 

	while ( fabs(mres0-mres)>tol && nitr<100 )
	{
		if (debug) printf("just enter while\n"); fflush(stdout);
		nitr++;
		// find the active constraints
		maxres = 0; 
		for (i = 0; i < M; ++i)
		{
			a1[0] = A1[i+0*M],  a1[1] = A1[i+1*M],  a1[2] = A1[i+2*M];
			a2[0] = A2[i+0*M],  a2[1] = A2[i+1*M],  a2[2] = A2[i+2*M];
			b1= B1[i];
	        b2= B2[i];
	        c[0] = C[i+0*M],  c[1] = C[i+1*M],  c[2] = C[i+2*M];
	        d= D[i];
	        res[i] = sqrt(   (a1[0]*X[0] + a1[1]*X[1] + a1[2]*X[2] + b1)   \
		                   * (a1[0]*X[0] + a1[1]*X[1] + a1[2]*X[2] + b1)   \
		                   + (a2[0]*X[0] + a2[1]*X[1] + a2[2]*X[2] + b2)   \
		                   * (a2[0]*X[0] + a2[1]*X[1] + a2[2]*X[2] + b2) ) \
		         	  /       (c[0]*X[0] +  c[1]*X[1] +  c[2]*X[2] + d);
		    if (maxres<res[i])
		    	maxres = res[i];
		}
        
//         printf("\n(%d) res = ", nitr);    tri_pVec(M,res); 
        
		mres0 = mres;
		mres = maxres;
        
//         printf("\n(%d) mres0 = %f, mres = %f", nitr, mres0, mres);
        
		if (mres>mres0)
		{
			flag = 0;
			break;
		}

		nActive = 0;
		for (i = 0; i < M; ++i)
 			if ( mres-res[i] < (*eps)*mres )
//            if ( mres-res[i] < 0.001 )
			{
				Active[nActive++] = i;
				if (nActive>=4)
				{
					break;
				}
			}

//        printf("\n(%d) nActive = %d, Active = ", nitr, nActive);    tri_pVecInt(nActive,Active); 
        if (debug) printf("before getgrad()\n");fflush(stdout);
		getneggrad( A1, A2, B1, B2, C, D, M, X, nActive, Active, grads );	
		if (debug) printf("after getgrad()\n");fflush(stdout);
		
//        printf("\n(%d) grads = ", nitr);   tri_pMatgrad(nActive,3,grads);

        if (debug) printf("before finddir()\n");fflush(stdout);
		//dir = c_finddir(grads, nActive, &r);


		meb(grads, nActive, dir);

		if (debug) printf("after finddir()\n");fflush(stdout);
		//free(grads);

       	norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
		if (norm<tol)
			break;

		
        
//         printf("\n(%d) norm = %f\n", nitr, norm);
        
		for (i=0;i<3;i++)
			dir[i] = dir[i] / norm;

//        printf("\n(%d) dir = ", nitr);    tri_pVec(3,dir); 
        if (debug) printf("before bound()\n"); fflush(stdout);
		//bound = getbound(C, D, M, X, dir);
		getbound(C, D, M, X, dir, &lb0, &ub0);
		if (debug) printf("after bound()\n");fflush(stdout);
		//lb0 = bound[0]; ub0 = bound[1];
		//free(bound);
        
//        printf("\n(%d) bound = [%f, %f]", nitr, lb0, ub0);
        if (debug) printf("before stepsize()\n");fflush(stdout);
		alpha = stepsize(A1, A2, B1, B2, C, D, M, X, dir, lb0, ub0);
        if (debug) printf("after stepsize()\n");fflush(stdout);
//        printf("\n(%d) alpha = %f\n", nitr, alpha);
	    if (alpha<tol)
			break;      

		for (i = 0; i < 3; ++i)
		{
			X0[i] = X[i];
			 X[i] = X[i] + alpha * dir[i];
		}
//         printf("\nflag(%d) = %d,  X = ", nitr, flag);  tri_pVec(3,X);
	}

	if (!flag)
	{
		X[0] = X0[0]; X[1] = X0[1]; X[2] = X0[2];
		mres = mres0;
	}

	// return X, mres, nitr;
	X_return[0] = X[0]; X_return[1] = X[1]; X_return[2] = X[2];   
 	*mres_return = mres;
 	*nitr_return = nitr;
}

/***************************************************************************/
double* getneggrad(double* A1, double* A2, double* B1, double* B2, double* C, double* D, \
	       int M, double* X,int nActive, int* Active, double* grads[4])
{
	//double *grads = (double *)malloc(nActive * 3 * sizeof(double));
	// grads = [g1 ; g2 ; ... ; gM];
	int i, j, actid;
	double a1[3];
	double a2[3];
	double b1;
    double b2;
    double c[3];
    double d;

    double a1xb1, a2xb2, cxd;
    double gx;
    double* gi;
    double norm;

	for (i = 0; i < nActive; ++i)
	{
		actid = Active[i];
		a1[0] = A1[actid+0*M],  a1[1] = A1[actid+1*M],  a1[2] = A1[actid+2*M];
		a2[0] = A2[actid+0*M],  a2[1] = A2[actid+1*M],  a2[2] = A2[actid+2*M];
         c[0] =  C[actid+0*M],   c[1] =  C[actid+1*M],   c[2] =  C[actid+2*M];
        b1= B1[actid];
        b2= B2[actid];
        d =  D[actid];

        a1xb1 = a1[0]*X[0] + a1[1]*X[1] + a1[2]*X[2] + b1;
        a2xb2 = a2[0]*X[0] + a2[1]*X[1] + a2[2]*X[2] + b2;
        cxd   = c[0]*X[0] + c[1]*X[1] + c[2]*X[2] + d;

		// GX = (A1*X+B1).^2 + (A2*X+B2).^2;
        gx = a1xb1*a1xb1 + a2xb2*a2xb2;

        gi = grads[i];
        for(j=0; j<3; j++)
        	//*(grads+i*3+j) = cxd/sqrt(gx) * (a1xb1*a1[j]+a2xb2*a2[j]) - sqrt(gx)*c[j];
        	gi[j] = cxd/sqrt(gx) * (a1xb1*a1[j]+a2xb2*a2[j]) - sqrt(gx)*c[j];

        //norm = sqrt( *(grads+i*3+0)**(grads+i*3+0)  + *(grads+i*3+1)**(grads+i*3+1) + \
        	         *(grads+i*3+2)**(grads+i*3+2) );
        norm = sqrt( gi[0]*gi[0] + gi[1]*gi[1] + gi[2]*gi[2] );

        for(j=0; j<3; j++)
        	//*(grads+i*3+j) = - *(grads+i*3+j) / norm;
        	gi[j] = - gi[j]/norm;
	}

	//return ;
}

/***************************************************************************/
void getbound(double* C, double* D, int M, double* X, double* dir, double* lb0, double* ub0)
{
    double lb = 0;
    double ub = Inf; /*/ just a large number;*/
    int i;
    double temp;
    double *result;
    double c[3];
    double d;

//     printf("C = "); tri_pMat2(M,3,C);
//     printf("D = "); tri_pVec(M,D);

    for(i=0;i<M;i++)
    {
        c[0] = C[i+0*M],  c[1] = C[i+1*M],  c[2] = C[i+2*M];
        d = D[i];
        temp = -(c[0]*X[0] + c[1]*X[1] + c[2]*X[2] + d) / (c[0]*dir[0] + c[1]*dir[1] + c[2]*dir[2]);

/*/         printf("temp = %f, cc.dir = %f \n", temp, (c1[0]*dir[0] + c1[1]*dir[1] + c1[2]*dir[2]));*/

        if( (c[0]*dir[0] + c[1]*dir[1] + c[2]*dir[2]) > 0 ) /*/ for lower bound*/
        {
            if(temp > lb)
                lb = temp;
//             printf("temp(%d) = %f, sign(%d) = %f, lb = %f\n", \
//             	i, temp, i, (c[0]*dir[0] + c[1]*dir[1] + c[2]*dir[2]), lb);
        }
        else if( (c[0]*dir[0] + c[1]*dir[1] + c[2]*dir[2]) < 0 )/*// for upper bound*/
        {
            if(temp < ub)
                ub = temp;
//             printf("temp(%d) = %f, sign(%d) = %f, ub = %f\n", \
//             	i, temp, i, (c[0]*dir[0] + c[1]*dir[1] + c[2]*dir[2]), ub);
        }

        
    }
    //result = (double*)malloc(2*sizeof(double));
    //result[0] = lb;
    //result[1] = ub;
    *lb0 = lb;
    *ub0 = ub;
/*//     printf("bound %f, %f \n", lb, ub);*/
    //return result;
}

/***************************************************************************/
double stepsize(double* A1, double* A2, double* B1, double* B2, double* C, double* D, \
	            int M, double* X, double* dir, double lb0, double ub0)
{
	if (lb0>=ub0)
		return 0;

	double a1[3];
	double a2[3];
	double b1;
    double b2;
    double c[3];
    double d;

	double tol1 = 0.000001; // 1e-6
	double tol2 = 0.00000001; // 1e-8
	double lb, ub, r;
	double Xn[3], Xp[3], Xc[3];
	double vni, vpi, vci;
	double vn, vp, vc;
	int i;

	ub = fmin(ub0, 10);
	lb = fmax(0, lb0);
	r = (ub+lb)/2;



	while (ub-lb>tol1)
	{
		r = (ub+lb)/2;
		for (i = 0; i < 3; ++i)
		{
			Xn[i] = X[i] + (r-tol2)*dir[i];
			Xp[i] = X[i] + (r+tol2)*dir[i];
			Xc[i] = X[i] + (r     )*dir[i];
		}

		vn = 0; vp = 0; vc = 0;
		for (i = 0; i < M; ++i)
		{
			a1[0] = A1[i+0*M],  a1[1] = A1[i+1*M],  a1[2] = A1[i+2*M];
			a2[0] = A2[i+0*M],  a2[1] = A2[i+1*M],  a2[2] = A2[i+2*M];
			b1= B1[i];
	        b2= B2[i];
	        c[0] = C[i+0*M],  c[1] = C[i+1*M],  c[2] = C[i+2*M];
	        d= D[i];
	        vni = sqrt(   (a1[0]*Xn[0] + a1[1]*Xn[1] + a1[2]*Xn[2] + b1)   \
	                    * (a1[0]*Xn[0] + a1[1]*Xn[1] + a1[2]*Xn[2] + b1)   \
	                    + (a2[0]*Xn[0] + a2[1]*Xn[1] + a2[2]*Xn[2] + b2)   \
	                    * (a2[0]*Xn[0] + a2[1]*Xn[1] + a2[2]*Xn[2] + b2) ) \
	         	  /        (c[0]*Xn[0] +  c[1]*Xn[1] +  c[2]*Xn[2] + d);
	        vpi = sqrt(   (a1[0]*Xp[0] + a1[1]*Xp[1] + a1[2]*Xp[2] + b1)   \
	                    * (a1[0]*Xp[0] + a1[1]*Xp[1] + a1[2]*Xp[2] + b1)   \
	                    + (a2[0]*Xp[0] + a2[1]*Xp[1] + a2[2]*Xp[2] + b2)   \
	                    * (a2[0]*Xp[0] + a2[1]*Xp[1] + a2[2]*Xp[2] + b2) ) \
	         	  /        (c[0]*Xp[0] +  c[1]*Xp[1] +  c[2]*Xp[2] + d);
	        vci = sqrt(   (a1[0]*Xc[0] + a1[1]*Xc[1] + a1[2]*Xc[2] + b1)   \
	                    * (a1[0]*Xc[0] + a1[1]*Xc[1] + a1[2]*Xc[2] + b1)   \
	                    + (a2[0]*Xc[0] + a2[1]*Xc[1] + a2[2]*Xc[2] + b2)   \
	                    * (a2[0]*Xc[0] + a2[1]*Xc[1] + a2[2]*Xc[2] + b2) ) \
	         	  /        (c[0]*Xc[0] +  c[1]*Xc[1] +  c[2]*Xc[2] + d);
	        if (vni>vn)
	            vn = vni;
	        if (vpi>vp)
	            vp = vpi;
	        if (vci>vc)
	            vc = vci;
		}

		if (vn>vc && vc>vp)
			lb = r;
		else if (vn<vc && vc<vp)
			ub = r;
		else
			break;
	}

	return r;
}

/***************************************************************************/
void tri_pMat1(int m, int n, double A[m][n])
{
	int i,j;
	printf("\n");
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printf("%.6f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/***************************************************************************/
void tri_pMat2(int m, int n, double* A)
{
	int i,j;
	printf("\n");
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printf("%.6f ", A[i+j*m]);
		}
		printf("\n");
	}
	printf("\n");
}

/***************************************************************************/
void tri_pMatgrad(int m, int n, double* grads)
{
	int i,j;
	printf("\n");
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printf("%.6f ", grads[i*n+j]);
		}
		printf("\n");
	}
	printf("\n");
}

/***************************************************************************/
void tri_pVec(int m, double b[m])
{
	int i;
	printf("\n");
	for (i = 0; i < m; ++i)
		printf("%.6f ", b[i]);
	printf("\n");
}

/***************************************************************************/
void tri_pVecInt(int m, int b[m])
{
	int i;
	printf("\n");
	for (i = 0; i < m; ++i)
		printf("%d ", b[i]);
	printf("\n");
}

/***************************************************************************/
/***************************************************************************/

static inline void cross(double* A, double* B, double* result)
{
	//(a2b3  -   a3b2,     a3b1   -   a1b3,     a1b2   -   a2b1)
	result[0] = A[1]*B[2] - A[2]*B[1];
	result[1] = A[2]*B[0] - A[0]*B[2];
	result[2] = A[0]*B[1] - A[1]*B[0];
}

static inline double dot(double* A, double* B)
{
	return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

static inline double norm3(double a1, double a2, double a3)
{
	return sqrt(a1*a1 + a2*a2 + a3*a3);
}

static inline double norm(double* A)
{
	return sqrt(A[0]*A[0] +A[1]*A[1] +A[2]*A[2]);
}

void circle(double* A, double* B, double* C, double* r, double* centre);

double meb3(double* A, double* B, double* C, double* centre);


void meb(double* grads[4], int N, double* dir)
{
	//printf("entre meb()\n");
	// N: number of gradients
	int d = 3; 	// dimension

	int i,j,f;
    
	if (N==1)
	{
		//printf("N==1\n");
		for (j = 0; j<d; j++)
			//dir[j] = *(grads+j);
			dir[j] = *(grads[0]+j);
        //*r = 0;
		return ;
	}
		
	if (N==2)
	{
		//printf("N==2\n");
		for (j = 0; j<d; j++)
			//dir[j] = (*(grads+0*d+j) + *(grads+1*d+j))/2;
			dir[j] = ( *(grads[0]+j) + *(grads[1]+j) ) / 2;
        //*r = norm3(dir[0]-*(grads+0*d+0), dir[1]-*(grads+0*d+1), dir[2]-*(grads+0*d+2));
		return ;
	}

	if (N==3)
	{
		//printf("N==3\n");
		double r;
		meb3(grads[0], grads[1], grads[2], dir);
		return ;
	}

	// N >= 4
	//printf("N>=4\n");
	f = 0;
	N = 4;
	int seq[4][4] = {{0,1,2,3}, {0,1,3,2}, {0,2,3,1}, {1,2,3,0}};
	//double pu[3];
	//double pv[3];
	//double pw[3];
	//double pr[3];
	double *pu, *pv, *pw, *pr;
	for (i = 0; i < N; ++i)
	{
		//for(j=0; j<d; j++)
		//{
		//	pu[j] = *(grads + seq[i][0]*d + j);
		//	pv[j] = *(grads + seq[i][1]*d + j);
		//	pw[j] = *(grads + seq[i][2]*d + j);
		//	pr[j] = *(grads + seq[i][3]*d + j);
		//}
		pu = grads[ seq[i][0] ];
		pv = grads[ seq[i][1] ];
		pw = grads[ seq[i][2] ];
		pr = grads[ seq[i][3] ];

		double r = meb3(pu,pv,pw, dir);

		if(norm3(dir[0]-pr[0],dir[1]-pr[1],dir[2]-pr[2]) <= r)
		{
			f=1;
			break;
		}
	}

	if(f==0)
	{
		//*r = 1;
		dir[0]=0; dir[1]=0; dir[2]=0; 
	}

	return ;
}

double meb3(double* A, double* B, double* C, double* centre)
{
	//printf("entre meb3()\n");
	int f=0;
	int seq[3][3] = { {0,1,2}, {0,2,1}, {1,2,0}};
	double* pu;
	double* pv;
	double* pr;
	double* all[3] = {A,B,C};
	int i,j;
	double r;
	for (i = 0; i < 3; ++i)
	{
		pu = all[ seq[i][0] ];
		pv = all[ seq[i][1] ];
		pr = all[ seq[i][2] ];
		for(j=0; j<3; j++)
			centre[j] = (pu[j]+pv[j])/2;

		r = norm3(centre[0]-pu[0],centre[1]-pu[1],centre[2]-pu[2]);
		if(norm3(centre[0]-pr[0],centre[1]-pr[1],centre[2]-pr[2]) <= r)
		{
			f=1;
			break;
		}
	}
	if (f==0)
		circle(pu,pv,pr, &r, centre);

	return r;
	//printf("exit meb3()\n");
}

void circle(double* A, double* B, double* C, double* r, double* centre)
{
	//printf("entre circle()\n");
	int i;
	double a[3] = {A[0]-C[0], A[1]-C[1], A[2]-C[2] };
	double b[3] = {B[0]-C[0], B[1]-C[1], B[2]-C[2] };
	double dota = dot(a,a);
	double na = sqrt(dota);
	double dotb = dot(b,b);
	double nb = sqrt(dotb);

	double aa[3] = {dota*b[0] - dotb*a[0], dota*b[1] - dotb*a[1], dota*b[2] - dotb*a[2]};

	double acrossb[3];
	cross(a,b,acrossb);

	double nab = norm( acrossb );

	double c_abc[3];
	cross(aa, acrossb, c_abc);

	*r = na*nb*norm3(a[0]-b[0], a[1]-b[1], a[2]-b[2]) / (2*nab);

	for(i=0; i<3; i++)
		centre[i] = c_abc[i]/(2*nab*nab) + C[i];
	//printf("exit circle()\n");
}














