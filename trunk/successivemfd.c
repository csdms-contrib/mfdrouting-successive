#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define niterations 100
#define fillincrement 0.01
#define oneoversqrt2 0.707106781186

float **topo,**toposave,*topovec,**depth,**slope,**flow,**flow1,**flow2,**flow3;
float **flow4,**flow5,**flow6,**flow7,**flow8,deltax;
int lattice_size_x,lattice_size_y,*iup,*idown,*jup,*jdown;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
      m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(n,arr,indx)
float arr[];
int indx[],n;
{
        unsigned long i,indxt,ir=n,itemp,j,k,l=1;
        int jstack=0,*istack;
        float a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void fillinpitsandflats(i,j)
int i,j;
{
    float min;

    min=topo[i][j];
    if (topo[iup[i]][j]<min) min=topo[iup[i]][j];
    if (topo[idown[i]][j]<min) min=topo[idown[i]][j];
    if (topo[i][jup[j]]<min) min=topo[i][jup[j]];
    if (topo[i][jdown[j]]<min) min=topo[i][jdown[j]];
    if (topo[iup[i]][jup[j]]<min) min=topo[iup[i]][jup[j]];
    if (topo[idown[i]][jup[j]]<min) min=topo[idown[i]][jup[j]];
    if (topo[idown[i]][jdown[j]]<min) min=topo[idown[i]][jdown[j]];
    if (topo[iup[i]][jdown[j]]<min) min=topo[iup[i]][jdown[j]];
    if ((topo[i][j]<=min)&&(i>1)&&(j>1)&&(i<lattice_size_x)&&(j<lattice_size_y))
     {topo[i][j]=min+fillincrement;
      fillinpitsandflats(i,j);
      fillinpitsandflats(iup[i],j);
      fillinpitsandflats(idown[i],j);
      fillinpitsandflats(i,jup[j]);
      fillinpitsandflats(i,jdown[j]);
      fillinpitsandflats(iup[i],jup[j]);
      fillinpitsandflats(idown[i],jup[j]);
      fillinpitsandflats(idown[i],jdown[j]);
      fillinpitsandflats(iup[i],jdown[j]);}
}

void mfdflowroute(i,j)
int i,j;
{    float tot;
 
     tot=0;
     if (topo[i][j]>topo[iup[i]][j]) 
      tot+=pow(topo[i][j]-topo[iup[i]][j],1.1);
     if (topo[i][j]>topo[idown[i]][j]) 
      tot+=pow(topo[i][j]-topo[idown[i]][j],1.1);
     if (topo[i][j]>topo[i][jup[j]]) 
      tot+=pow(topo[i][j]-topo[i][jup[j]],1.1);
     if (topo[i][j]>topo[i][jdown[j]]) 
      tot+=pow(topo[i][j]-topo[i][jdown[j]],1.1);
     if (topo[i][j]>topo[iup[i]][jup[j]]) 
      tot+=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[iup[i]][jdown[j]]) 
      tot+=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jup[j]]) 
      tot+=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jdown[j]]) 
      tot+=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[iup[i]][j]) 
      flow1[i][j]=pow(topo[i][j]-topo[iup[i]][j],1.1)/tot; 
       else flow1[i][j]=0;
     if (topo[i][j]>topo[idown[i]][j]) 
      flow2[i][j]=pow(topo[i][j]-topo[idown[i]][j],1.1)/tot; 
       else flow2[i][j]=0;
     if (topo[i][j]>topo[i][jup[j]]) 
      flow3[i][j]=pow(topo[i][j]-topo[i][jup[j]],1.1)/tot; 
       else flow3[i][j]=0;
     if (topo[i][j]>topo[i][jdown[j]]) 
      flow4[i][j]=pow(topo[i][j]-topo[i][jdown[j]],1.1)/tot; 
       else flow4[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jup[j]]) 
      flow5[i][j]=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1)/tot;
       else flow5[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jdown[j]]) 
      flow6[i][j]=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
       else flow6[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jup[j]]) 
      flow7[i][j]=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1)/tot;
       else flow7[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jdown[j]]) 
      flow8[i][j]=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
       else flow8[i][j]=0;
     flow[iup[i]][j]+=flow[i][j]*flow1[i][j];
     flow[idown[i]][j]+=flow[i][j]*flow2[i][j];
     flow[i][jup[j]]+=flow[i][j]*flow3[i][j];
     flow[i][jdown[j]]+=flow[i][j]*flow4[i][j];
     flow[iup[i]][jup[j]]+=flow[i][j]*flow5[i][j];
     flow[iup[i]][jdown[j]]+=flow[i][j]*flow6[i][j];
     flow[idown[i]][jup[j]]+=flow[i][j]*flow7[i][j];
     flow[idown[i]][jdown[j]]+=flow[i][j]*flow8[i][j];
}

void calculatealongchannelslope(i,j)
int i,j;
{    float down;

     down=0;
     if (topo[iup[i]][j]-topo[i][j]<down) down=topo[iup[i]][j]-topo[i][j];
     if (topo[idown[i]][j]-topo[i][j]<down) down=topo[idown[i]][j]-topo[i][j];
     if (topo[i][jup[j]]-topo[i][j]<down) down=topo[i][jup[j]]-topo[i][j];
     if (topo[i][jdown[j]]-topo[i][j]<down) down=topo[i][jdown[j]]-topo[i][j];
     if ((topo[iup[i]][jup[j]]-topo[i][j])*oneoversqrt2<down)
      down=(topo[iup[i]][jup[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[idown[i]][jup[j]]-topo[i][j])*oneoversqrt2<down) 
      down=(topo[idown[i]][jup[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[iup[i]][jdown[j]]-topo[i][j])*oneoversqrt2<down) 
      down=(topo[iup[i]][jdown[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[idown[i]][jdown[j]]-topo[i][j])*oneoversqrt2<down) 
      down=(topo[idown[i]][jdown[j]]-topo[i][j])*oneoversqrt2;
     slope[i][j]=fabs(down)/deltax;
}

main()
{    FILE *fp1,*fp2;
     float runoff,manningsn,*topovec;
     char inputfile[100],outputfile[100];
     int i,j,k,t,*topovecind;

     scanf("%s",&inputfile);
     scanf("%s",&outputfile);
     scanf("%d",&lattice_size_x);
     scanf("%d",&lattice_size_y);
     scanf("%f",&deltax);  /* m */
     scanf("%f",&runoff);  /* m */
     scanf("%f",&manningsn);
     fp1=fopen(inputfile,"r");
     fp2=fopen(outputfile,"w");
     setupgridneighbors();
     topo=matrix(1,lattice_size_x,1,lattice_size_y);
     toposave=matrix(1,lattice_size_x,1,lattice_size_y);
     topovec=vector(1,lattice_size_x*lattice_size_y);
     topovecind=ivector(1,lattice_size_x*lattice_size_y);
     depth=matrix(1,lattice_size_x,1,lattice_size_y);
     slope=matrix(1,lattice_size_x,1,lattice_size_y);
     flow=matrix(1,lattice_size_x,1,lattice_size_y);
     flow1=matrix(1,lattice_size_x,1,lattice_size_y);
     flow2=matrix(1,lattice_size_x,1,lattice_size_y);
     flow3=matrix(1,lattice_size_x,1,lattice_size_y);
     flow4=matrix(1,lattice_size_x,1,lattice_size_y);
     flow5=matrix(1,lattice_size_x,1,lattice_size_y);
     flow6=matrix(1,lattice_size_x,1,lattice_size_y);
     flow7=matrix(1,lattice_size_x,1,lattice_size_y);
     flow8=matrix(1,lattice_size_x,1,lattice_size_y);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fscanf(fp1,"%f",&topo[i][j]);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fillinpitsandflats(i,j); 
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {toposave[i][j]=topo[i][j];
        calculatealongchannelslope(i,j);}
     for (k=1;k<=niterations;k++)
      {printf("iteration = %d/%d\n",k,niterations);
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         {flow[i][j]=pow(runoff,1.66667)*sqrt(slope[i][j])*deltax/manningsn;
          topovec[(j-1)*lattice_size_x+i]=topo[i][j];}
          indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
          t=lattice_size_x*lattice_size_y+1;
          while (t>1)
           {t--;
            i=(topovecind[t])%lattice_size_x;
            if (i==0) i=lattice_size_x;
            j=(topovecind[t])/lattice_size_x+1;
            if (i==lattice_size_x) j--;
            mfdflowroute(i,j);
            if (slope[i][j]>0.001) depth[i][j]=pow(flow[i][j]*manningsn/
             (deltax*sqrt(slope[i][j])),0.6);
             else depth[i][j]=0;}
          for (j=1;j<=lattice_size_y;j++)
           for (i=1;i<=lattice_size_x;i++)
            if ((toposave[i][j]+depth[i][j])>topo[i][j])
		 topo[i][j]=toposave[i][j]+depth[i][j]/niterations;
          for (j=1;j<=lattice_size_y;j++)
           for (i=1;i<=lattice_size_x;i++)
            fillinpitsandflats(i,j);}
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fprintf(fp2,"%f\n",topo[i][j]-toposave[i][j]);
}
