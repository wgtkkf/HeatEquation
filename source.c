/* ********************************************* */
/* Coded by Takuro TOKUNAGA                      */
/* 2D heat conduction equation solved by FEM     */
/* Liner interpolation                           */
/* For final project of Adv. Conduction class    */
/* Last Modified: March 14, 2017                 */
/* About this code:                              */
/* required files:                               */
/* 1. rectangle1.msh: number of coord and nord   */
/* 2. rectangle2.msh: cord information           */
/* 3. rectangle3.msh: nord information           */
/* 4. parameters.txt                             */
/* 5. rectangle-t.bc                             */
/* 6. rectangle-num.bc                           */
/* 7. rectangle-n.bc                             */
/* ********************************************* */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define size0 3
#define size1 2
/* Need change here when change mesh */
#define nmax 287 /* 286 + 1 */
#define emax 501 /* 500 + 1 */

#define bcmax1 50
#define bcmax2 50
#define bcmax3 50
#define bcmax4 50
#define bcmax5 50
#define bcmax6 50
#define bcmax  50

void start();
void end();
void check();

/* main program */
int main(void)
{
	/* counters */
	int i, j, k, n;
	int counter=0;
	int stepnum=0;
	double time, tmax=0.001;
	
	/* matrixs */
	double area, area2, a2, area3, area12;
	int n1, n2, n3;
    int elem;
    double nop[emax][size0]={};
	double cord[nmax][size1]={};
	double b[size0]={}, c[size0]={}, dd[size0][size0]={};
	double eb[emax][size0]={}, ec[emax][size0]={};
	double x1, x2, x3, y1, y2, y3;
	double emm[emax][size0]={}, ehx[emax][size0][size0]={}, ehy[emax][size0][size0]={}, ess[emax][size0][size0]={};
	double essxbb[emax][size0][size0]={}, essxbc[emax][size0][size0]={}, essycc[emax][size0][size0]={}, essycb[emax][size0][size0]={};
	double lmm[nmax]={}, ilmm[nmax]={};
    double ass[emax][3][3]={};
	double ekx[emax][size0][size0][size0]={}, eky[emax][size0][size0][size0]={};
		
	double nbc1, nnbc1[bcmax1];
	double nbc2, nnbc2[bcmax2];
	double nbc3, nnbc3[bcmax3];
	double nbc4, nnbc4[bcmax4];
	double nbc5, nnbc5[bcmax5];
	double nbc6, nnbc6[bcmax6];
	double vbc1[bcmax1], vbc2[bcmax2], vbc3[bcmax3], vbc4[bcmax4], vbc5[bcmax5], vbc6[bcmax6];
    double tbc1[bcmax1];
		
	/* Parameters for heat conduction equation */
	FILE *fp0;
    double lambda; /* thermal conductiviy */
	double rho;    /* density             */
    double shc;    /* specific heat C     */
    double c1;     /* Coefficient         */
	
	/* element informations */
	FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;
	char s[20];
	int np;
	int ne;
	double dt;
    
    /* temperature */
    double t0[nmax]={};
    double tt0[nmax]={};
    double t1[nmax]={};
    double t00=0;
    double t11=0;
		
	/* cg method */
	double r[nmax];
	double p[nmax];
	double vectorb[nmax];
	double b2;
	double eps, eps1, eps2;
	double delta;
	int kend, icount;
	double ap[nmax];
	double rur0, rur1;
	double res2;
	double pap;
	double alpha;
	double beta;
	double v0, v1;
	double tmp;
	double vmean;
	double res;
	
	double nbcn;
	double nnbcn[bcmax][size1]={};/* check the size of table */
	double vbcn[bcmax];           /* check the size of table */
	
	int N, M, PG;
	double LL;
	
	double std = 10E-8;/* convergence */
	
	/* output file */
	char filename[256];
    char filename1[256];
	FILE *fp;
    FILE *fp_excel;
	div_t d;
	
	/* cpu time */
	clock_t stime, ftime;
	
	/* start of the main program */
	start();
	stime=clock();

	/* 1. read the parameters for newton flow */
	fp0=fopen("parameters.txt","r");
	if(fp0==NULL)
	{
		printf("failure\n");
        return -1;
	}
		
	fscanf(fp0,"%lf",&lambda);
	fscanf(fp0,"%lf",&rho);
	fscanf(fp0,"%lf",&shc);
	fclose(fp0);
		  
	printf("parameters for heat conduction\n");
    printf("thermal conductivity:%lf\n", lambda);
	printf("density:%lf\n", rho);
    printf("specific heat:%lf\n", shc);
	printf("\n");

	/* 2. read the element information */
	/* 2.1 read the mesh information */
	fp1=fopen("rectangle1.msh","r");
	if(fp1==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp2=fopen("rectangle2.msh","r");
	if(fp2==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp3=fopen("rectangle3.msh","r");
	if(fp3==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	
	/* *fp1: rectangle1.msh */
	while(fgets(s, 20, fp1)!=NULL){
		if(counter==1)
		{
			np = atoi(s);
		}else if(counter==2){
			ne = atoi(s);
		}else if(counter==3){
			dt = atof(s);
		}
		counter++;
	};
	counter=0;
	
	printf("number of the nodes:%d\n", np);
	printf("number of the elements:%d\n", ne);
	printf("dt:%lf\n", dt);
	
	/* show the cord information */
	/* read the coordinates of the mesh */
	for(i=1;i<=np;i++)
	{
		fscanf(fp2, "%d %lf %lf", &n, &cord[i][0], &cord[i][1]);
	}

    /* show the node information */
	for(i=1;i<=ne;i++)
	{
		fscanf(fp3, "%d %lf %lf %lf", &n, &nop[i][0], &nop[i][1], &nop[i][2]);
	}
    printf("test output\n");
    
    fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	/* 3. read the initial conditions */
	/* initialize */
	for(i=1;i<=np;i++)
	{
		t0[i] = 20.0;
		t1[i] = 20.0;
	}
			
	/* 4. read the boundary conditions */
	fp4=fopen("rectangle-num.bc","r");
	if(fp4==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp5=fopen("rectangle-t.bc","r");
	if(fp5==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp6=fopen("rectangle-n.bc","r");
	if(fp6==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	
	/* read the number of condition for temperature and N.B.*/
	fscanf(fp4, "%lf", &nbc1);
	fscanf(fp4, "%lf", &nbc2);
	fscanf(fp4, "%lf", &nbc3);
	fscanf(fp4, "%lf", &nbcn);
	fscanf(fp4, "%lf", &vmean);

	/* boundary condition for T */
	for(i=1;i<=(int)nbc1;i++)
	{
		fscanf(fp5, "%lf %lf", &nnbc1[i], &tbc1[i]);
	}
	/* boundary condition for natural B.C */
	for(i=1;i<=(int)nbcn;i++)
	{
		fscanf(fp6, "%d %d %d", &N, &M, &PG);
		if(M==1)
		{
			n1 = nop[n][1];
			n2 = nop[n][2];
		}else if(M==2)
		{
			n1 = nop[n][0];
			n2 = nop[n][2];
		}else if(M==3)
		{
			n1 = nop[n][0];
			n2 = nop[n][1];
		}
		nnbcn[i][0]=n1;
		nnbcn[i][1]=n2;
		
		x1 = cord[n1][0];
		y1 = cord[n1][1];
		x2 = cord[n2][0];
		y2 = cord[n2][1];
		LL = sqrt(pow((x1-x2),2) + pow((y1-y2),2));
		vbcn[i] = PG*LL;
	}
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);

	/* 5. calculation of the element matrix */
	for(elem=1;elem<=ne;elem++)
	{
		/* initialization */
		for(i=0;i<size0;i++)
		{
			emm[elem][i] = 0.0;
			for(j=0;j<size0;j++)
			{
				ehx[elem][i][j] = 0.0;
				ehy[elem][i][j] = 0.0;
				ess[elem][i][j] = 0.0;
                ass[elem][i][j] = 0.0;
				if(i==j)
				{
					dd[i][j] = 2.0;
				}else{
					dd[i][j] = 1.0;
				}
				for(k=0;k<size0;k++)
				{
					ekx[elem][i][j][k] = 0.0;
					eky[elem][i][j][k] = 0.0;
				}
			}
		}
		
		/* calculation area of element */
		/* n1 -- number of 1st node */
		/* n2 -- number of 1st node */
		/* n3 -- number of 1st node */
		/* area -- area of element */
		n1 = nop[elem][0];
		n2 = nop[elem][1];
		n3 = nop[elem][2];
		x1 = cord[n1][0];
		y1 = cord[n1][1];
		x2 = cord[n2][0];
		y2 = cord[n2][1];
		x3 = cord[n3][0];
		y3 = cord[n3][1];
		
		area2 = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
		area = 0.5*area2;		
		
		if(area<0)
		{
			printf("error: area is less than zero\n");
			exit(1);
		}
		
		a2 = 1/area2;
		area3 = area/3;
		area12 = area/12;
						
		/* differential of shape function */
		b[0] = (y2-y3)*a2;
		b[1] = (y3-y1)*a2;
		b[2] = (y1-y2)*a2;
		c[0] = (x3-x2)*a2;
		c[1] = (x1-x3)*a2;
		c[2] = (x2-x1)*a2;
        
		eb[elem][0] = b[0];
		eb[elem][1] = b[1];
		eb[elem][2] = b[2];
		ec[elem][0] = c[0];
		ec[elem][1] = c[1];
		ec[elem][2] = c[2];
		
		/* calculation of matrix */
		/* Hx, Hy */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				ehx[elem][i][j] = area3*b[j];
				ehy[elem][i][j] = area3*c[j];
			}
		}
		
		/* Kxx, Kyy */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				for(k=0;k<size0;k++)
				{
					ekx[elem][i][j][k] = area12*dd[i][j]*b[k];
					eky[elem][i][j][k] = area12*dd[i][j]*c[k];
				}
			}
		}
		
		/* S */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				ess[elem][i][j] = area*(b[i]*b[j]+c[i]*c[j]);
			}
		}
		
		/* Sxbb */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essxbb[elem][i][j] = 2*area*(b[i]*b[j]);
			}
		}

		/* Sxbc */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essxbc[elem][i][j] = area*(b[i]*c[j]);
			}
		}
				
		/* Sycc */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essycc[elem][i][j] = 2*area*(c[i]*c[j]);
			}
		}
		
		/* Sycb */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essycb[elem][i][j] = area*(c[i]*b[j]);
			}
		}
        
        /* A */
        for(i=0;i<size0;i++)
        {
            for(j=0;j<size0;j++)
            {
                ass[elem][i][j] = area*(b[i]*b[j]+c[i]*c[j]);
            }
        }
        
		/*  -                    */
		/* [M]lumped mass matrix */
		for(i=0;i<size0;i++)
		{
			emm[elem][i] = area3;
			n = nop[elem][i];
			lmm[n] = lmm[n] + area3;
		}
	}
	
	/* inverse matrix of the element matrix */
	for(i=1;i<=np;i++)
	{
		ilmm[i] = 1.0/lmm[i];
	}	
		
	/* start calculation */
	start();
    
    /* 4. start time loop */
    
    /* coefficient */
    c1 = (lambda*dt)/(rho*shc);
	while(time<=tmax){
		
		if(tmax<time)
		{
			break;
		}
		
        /* 5. calculate temperature */
        /* calculate r.h.s temperature (T) */

		for(i=1;i<=np;i++)
		{
            /* initialization */
			tt0[i] = 0.0;
		}
        
		for(n=1;n<=ne;n++)
		{
            for(i=0;i<size0;i++)
			{
				n1 = nop[n][i];
				for(j=0;j<size0;j++)
				{
					n2 = nop[n][j];
					tt0[n1] = tt0[n1] - ass[n][i][j]*t1[n2];
                }
            }
        }
        /* boundary condition */
        for(i=1;i<=nbc1;i++)
        {
            tt0[(int)nnbc1[i]] = 0.0;
        }
        
		
		/* boundary condition */
		for(i=1;i<=nbc1;i++)
		{
			t0[(int)nnbc1[i]] = tbc1[i];
		}
        
        /* Explicit Method */
        for(i=1;i<=np;i++)
        {
            t1[i] = t0[i] + c1*tt0[i]*ilmm[i];
        }
        
        /* 6. Update of Temperature */
        for(i=1;i<=np;i++)
        {
            t0[i] = t1[i];
        }
        
        
        if(time<tmax)
		{
            
			/* format of AVese */		
			d = div(stepnum, 500);
			
			/* write to files */
			if(stepnum==0)
			{
				
				/* 7. output of the result step0 */
				sprintf(filename,"./data/step%d.vtk",stepnum);
                /* output of the result step0 excel */
                sprintf(filename1,"./excel/excel.dat");
                
				/* 7.1 creation of files */
				if((fp=fopen(filename, "w"))==NULL){
					fprintf(stderr, "%s\n", filename);
				}
                /* creation of files for excel */
                if((fp_excel=fopen(filename1, "w"))==NULL){
                    fprintf(stderr, "%s\n", filename1);
                }
                
                /* 7.2 write words to files */
				/* format for PARAVIEW, temperature */
				fprintf(fp, "# vtk DataFile Version 2.0\n");
				fprintf(fp, "temperature\n");
				fprintf(fp, "ASCII\n");
				fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
				fprintf(fp, "POINTS %d float\n", np);
	
                /* cord information */
				for(i=1;i<=np;i++)
				{
					fprintf(fp, "%lf %lf 0.0\n", cord[i][0], cord[i][1]);
				}
				/* nord information */
				fprintf(fp, "CELLS %d %d\n", ne, 4*(emax-1));
				for(i=1;i<=ne;i++)
				{
					fprintf(fp, "3 %d %d %d\n", (int)nop[i][0]-1, (int)nop[i][1]-1, (int)nop[i][2]-1);
				}
				
				fprintf(fp, "CELL_TYPES %d\n", ne);
				for(i=1;i<=ne;i++)
				{
					fprintf(fp, "5\n");
				}
				
				/* temperature data */
				fprintf(fp, "POINT_DATA %d\n", np);
				fprintf(fp, "SCALARS point_scalars float\n");
                fprintf(fp, "LOOKUP_TABLE default\n");
    
				for(i=1;i<=np;i++)
				{
					fprintf(fp, "0.0\n");
				}
				fprintf(fp, "CELL_DATA %d\n", ne);
                
                /* excel data */
                fprintf(fp_excel, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", stepnum, time, t1[91], t1[92],t1[93], t1[94], t1[95], t1[96], t1[97], t1[98], t1[99],t1[100]);
            
			}else if(stepnum>0)
			{
								
				/* 7. output of the result */
				sprintf(filename,"./data/step%d.vtk",stepnum);
							
				/* 7.1 creation of files */
                if((fp=fopen(filename, "w"))==NULL){
                fprintf(stderr, "%s\n", filename);
                }
                
                /* creation of files for excel */
                if((fp_excel=fopen(filename1, "a"))==NULL){
                    fprintf(stderr, "%s\n", filename1);
                }
                
				/* 7.2 write words to files */
				/* format of PARAVIEW */
				 fprintf(fp, "# vtk DataFile Version 2.0\n");
				 fprintf(fp, "velosity\n");
				 fprintf(fp, "ASCII\n");
				 fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
				 fprintf(fp, "POINTS %d float\n", np);
                
				/* cord information */
				for(i=1;i<=np;i++)
				{
					fprintf(fp, "%lf %lf 0.0\n", cord[i][0], cord[i][1]);
				}
				/* nord information */
				fprintf(fp, "CELLS %d %d\n", ne, 4*(emax-1));
				for(i=1;i<=ne;i++)
				{
					fprintf(fp, "3 %d %d %d\n", (int)nop[i][0]-1, (int)nop[i][1]-1, (int)nop[i][2]-1);
				}
				
				fprintf(fp, "CELL_TYPES %d\n", ne);
				for(i=1;i<=ne;i++)
				{
					fprintf(fp, "5\n");
				}
				
                /* temperature data */
                fprintf(fp, "POINT_DATA %d\n", np);
                fprintf(fp, "SCALARS point_scalars float\n");
                fprintf(fp, "LOOKUP_TABLE default\n");
                
                for(i=1;i<=np;i++)
                {
                    fprintf(fp, "%lf\n", t1[i]);
                }
                fprintf(fp, "CELL_DATA %d\n", ne);
                
                /* excel data */
                fprintf(fp_excel, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", stepnum, time, t1[91], t1[92],t1[93], t1[94], t1[95], t1[96], t1[97], t1[98], t1[99],t1[100]);
			}
			
			fclose(fp);
            fclose(fp_excel);

			/* prpceed time */				
			time = time + dt;

			/* proceed the step */
            printf("step:%d\n", stepnum);
			stepnum++;
		}
    }

	/* 9. end of the maim program */
	ftime = clock();
	printf("%lf second\n", (double)(ftime-stime)/(double)CLOCKS_PER_SEC);
	end();
	return 0;
}

/* functions */
void start()
{
	printf("start of the calculation\n\n");	
}

void end()
{
	printf("end of the calculation\n");	
}

void check()
{
	printf("proceeding..\n");	
}
