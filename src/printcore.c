/*******************************************************************************
* \note TERMS OF USE:
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation. This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The user relies on the
* software, documentation and results solely at his own risk.
*
* \file     printcore.c
* \brief    Generates EPS image representing a reactor core with color
			coded per-assembly information such as power density. 
			For help use printcore -h
* \author   Radim Vocka
*
*******************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <limits.h>

typedef struct{
  int ncols;
  int nvals;
  int *precision;
} vf_header;

typedef struct{
  char img_type;
  char print_nan;
} imagecontrol;


void PrintProlog(FILE *psat,double pitch,int core_size,int symmetry,imagecontrol *IC,int *dbb,const char layout);
int* ReadConf(char *conffile,int *core_size,int *symmetry,double *pitch,char **outfile,
	      double *dxy,char *layout);
int GetLine(FILE *cist,char *buff,int buff_size,const char eol);
int Line2IntField(char *buff,int *field,int fieldsize);
void PrintCore1(FILE *psat,int core_size,int *Coremap,char *vfile,char *colorfile,double *dxy,double *printcolor,
		imagecontrol *IC,vf_header *vfhp,const char layout);
void PrintCore6(FILE *psat,int n_anneaux,int *Coremap,char *vfile,char *colorfile,double *dxy,double *printcolor,
		imagecontrol *IC,vf_header *vfhp);
void PrintSampleConffile();
void* ReadVfile(char *vfile,int *nc,int *nv,int **precision,vf_header *vfhp);
void* ReadColorFile(char *colorfile,int *nv);



int main(int argc,char *argv[])
{
  FILE *psat;
  int *colormap,i,j;
  int core_size,symmetry;
  char layout;
  double pitch,dxy[4],*printcolor=NULL;
  int dbb[4];
  char *conffile=NULL,*outfile,*outfile_cl=NULL,*vfile=NULL,*colorfile=NULL;
  imagecontrol IC;
  vf_header vfh,*vfhp=NULL;
  char exepath[PATH_MAX] = {0};
  char *c;
  for(i=0;i<4;++i){
    dbb[i]=0;
    dxy[i]=0;
  }
  IC.img_type='c';
  IC.print_nan='y';
  readlink("/proc/self/exe", exepath, sizeof(exepath));
  c=exepath+strlen(exepath);
  while(*c!='/') --c; 
  *c='\0';
  /* read the command line args */
  for(i=1;i<argc;++i){
    if(!strcmp(argv[i],"-h")){
      printf("prints the core in eps format\n");
      printf("Usage: printcore -f conffile [-h] [-hc] [-ae]\n");
      printf("       -v file ... file with values to print (ncols-1,nbvals,precisions[ncols-1])\n");
      printf("       -o file ... outputfile - overrides that of the conffile\n");
      printf("       -cf file ... file with colors for the map - color for each assembly\n");
      printf("                    structure: nbvals, (i,color)*nbvals, color=<0,12>\n");
      printf("       -h  ... prints this help\n");
      printf("       -hc ... prints out sample conffile\n");
      printf("       -hd ... input the header of the data file on cmd line (nbr columns, nbr lines, precision*nbr_columns)\n");
      printf("       -ae ... takes the conffile for the ete assembly (.../prog/printcore/conf/ete_assembly.txt)\n");
      printf("       -ae6 ... takes the conffile for the ete1/6 assembly (.../prog/printcore/conf/ete_assembly6.txt)\n");
      printf("       -c  min max column ... image with color scale\n");
      printf("       -dbb  dx0 dx1 dy0 dy1 .. delta bounding box - to finetune ...\n");
      printf("       -shift  dx dy .. shift the core image - to finetune ...\n");
      printf("       -a  pring an image of the assembly (instead of core)\n");
      printf("       -nn  no NaNs\n");
      exit(0);
    }
    else if(!strcmp(argv[i],"-hc")){
      PrintSampleConffile();
      exit(0);
    }
    else if(!strcmp(argv[i],"-f")){
      if(conffile) free(conffile);
      conffile=(char*)malloc(strlen(argv[++i])+1);
      strcpy(conffile,argv[i]);
    }
    else if(!strcmp(argv[i],"-o")){
      outfile_cl=(char*)malloc(strlen(argv[++i])+1);
      strcpy(outfile_cl,argv[i]);
    }
    else if(!strcmp(argv[i],"-ae")){
      if(conffile) free(conffile);
      conffile=(char*)malloc(strlen(exepath)+strlen("/../conf/conf_ete_assembly.txt")+1);
      sprintf(conffile,"%s%s",exepath,"/../conf/conf_ete_assembly.txt");
    }
    else if(!strcmp(argv[i],"-ae6")){
      if(conffile) free(conffile);
      conffile=(char*)malloc(strlen(exepath)+strlen("/../conf/conf_ete_assembly6.txt")+1);
      sprintf(conffile,"%s%s",exepath,"/../conf/conf_ete_assembly6.txt");
    }
    else if(!strcmp(argv[i],"-v")){
      vfile=(char*)malloc(strlen(argv[++i])+1);
      strcpy(vfile,argv[i]);
    }
    else if(!strcmp(argv[i],"-cf")){
      colorfile=(char*)malloc(strlen(argv[++i])+1);
      strcpy(colorfile,argv[i]);
    }
    else if(!strcmp(argv[i],"-c")){
      printcolor=malloc(3*sizeof(double));
      if(!printcolor){
	fprintf(stderr,"printcolor malloc failed. Aborting.\n");
	exit(1);
      }
      printcolor[0]=atof(argv[++i]);
      printcolor[1]=atof(argv[++i]);
      printcolor[2]=atof(argv[++i]);
    }
    else if(!strcmp(argv[i],"-dbb")){
      for(j=0;j<4;++j) dbb[j]=atof(argv[++i]);
    }
    else if(!strcmp(argv[i],"-shift")){
      for(j=2;j<4;++j) dxy[j]=atof(argv[++i]);
    }
    else if(!strcmp(argv[i],"-a")){
      IC.img_type='a';
    }
    else if(!strcmp(argv[i],"-nn")){
      IC.print_nan='n';
    }
    else if(!strcmp(argv[i],"-hd")){
      vfhp=&vfh;
      vfh.ncols=atoi(argv[++i]);
      vfh.nvals=atoi(argv[++i]);
      vfh.precision=malloc(sizeof(int)*vfh.ncols);
      for(j=0;j<vfh.ncols;++j)
	vfh.precision[j]=atoi(argv[++i]);
    }
    else{
      fprintf(stderr,"Unknown command line option: %s.Aborting.\n",argv[i]);
      exit(1);
    }
  }  
  if(!conffile){
    fprintf(stderr,"Conffile not specified (try -h for help). Aborting.\n");
    exit(1);
  }
  /* execute the program itself */
  colormap=ReadConf(conffile,&core_size,&symmetry,&pitch,&outfile,dxy,&layout);
  if(outfile_cl){
    free(outfile);
    outfile=outfile_cl;
  }
  psat=fopen(outfile,"w");
  PrintProlog(psat,pitch,core_size,symmetry,&IC,dbb,layout);
  if(symmetry==6) PrintCore6(psat,(core_size-1)/2,colormap,vfile,colorfile,dxy,printcolor,&IC,vfhp);
  else PrintCore1(psat,core_size,colormap,vfile,colorfile,dxy,printcolor,&IC,vfhp,layout);
  fclose(psat);
  return 0;
}


/*  returns first j position in hexagon with n_ann for line i  */
int minj(int i,int n_ann)
{
  if(i<n_ann) return 0;
  return i-n_ann;
}

/*  returns last j position in hexagon with n_ann for line i  */
int maxj(int i,int n_ann)
{
  if(i<n_ann) return n_ann+i;
  return 2*n_ann;
}

#ifdef _HPUX
int Min(int i,int j)
#else
inline int Min(int i,int j)
#endif
{
  if(i<j) return i;
  return j;
}

#ifdef _HPUX
int Max(int i,int j)
#else
inline int Max(int i,int j)
#endif
{
  if(i<j) return j;
  return i;
}


/* layout - (s)quare/(h)ex */
int* ReadConf(char *fname,int *coresize,int *sym,double *pit,char **outfile,double *dxy,char *layout)  
{
  int    *Colormap,alloc_size,core_size,NN,buff_size,i,j,n,pos,i_buff,repet;
  int    fieldsize,max_i,n_anneaux,symmetry;  
  char   *buff,*substr,c_flag,layout_loc[10];
  FILE   *cist;
  double pitch;
  
  /*  malloc the buffer */
  buff_size=100; buff=(char*)malloc(buff_size);
  if(buff==NULL){fprintf(stderr,"InitColormap: buff malloc failed. Aborting.\n");exit(1);}
  /*  read the geometry data */
  cist=fopen(fname,"r");
  if(!cist){
    fprintf(stderr,"ReadConf: file %s cannot be opened. Aborting.\n",fname);exit(1);
  }
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  *outfile=(char*)malloc(strlen(buff)+1);
  if(!*outfile){
    fprintf(stderr,"malloc of outfile failed. Aborting.\n");
    exit(1);
  }
  strcpy(*outfile,buff);
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  sscanf(buff,"%lf",&pitch);
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  sscanf(buff,"%s",layout_loc);
  *layout=layout_loc[0];
  if(*layout!='s' && *layout!='h'){
    fprintf(stderr,"unknown layout type (%c). Aborting.\n",*layout);
    exit(1);
  }
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  sscanf(buff,"%d",&n_anneaux);
  if(*layout=='h')
    alloc_size=2*n_anneaux+2;
  else
    alloc_size=n_anneaux+1;
  core_size=alloc_size-1;  
  NN=alloc_size*alloc_size;
  /*  core symmetry (1=full core, 6=1/6 core) */
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  sscanf(buff,"%d",&symmetry);
  if(symmetry !=1 && symmetry!=6){    
    fprintf(stderr,"unknown symmetry type (%d). Aborting.\n",symmetry);
    exit(1);
  }
  Colormap=(int*)malloc(NN*sizeof(int));
  if(!Colormap){
    fprintf(stderr,"Colormap malloc failed. Aborting.\n");
    exit(1);
  }  
  for(i=0;i<NN;++i) Colormap[i]=0;
  /*  read the colormap */
  if(*layout=='h'){
    if(symmetry==1) max_i=core_size-1;
    else max_i=(int)(core_size/2);
    for(i=max_i;i>=0;--i){
      while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
      j=minj(i,n_anneaux);
      if(symmetry==1) fieldsize=maxj(i,n_anneaux)-j+1;
      else fieldsize=i+1;      
      pos=i*alloc_size+j;      
      if(fieldsize!=(n=Line2IntField(buff,Colormap+pos,fieldsize))){
	fprintf(stderr,"InitColormap: err read only %d values out of ??. Aborting.\n",
		n);
	exit(1);
      }
    }
  }
  else{
    for(i=core_size-1;i>=0;--i){
      while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
      pos=i*alloc_size;      
      fieldsize=core_size;
      if(fieldsize!=(n=Line2IntField(buff,Colormap+pos,fieldsize))){
	fprintf(stderr,"InitColormap: err read only %d values out of ??. Aborting.\n",n);
	exit(1);
      }
    }
  }
  while(!GetLine(cist,buff,buff_size,'\n') || buff[0]=='#'){}
  sscanf(buff,"%lf%lf",dxy,dxy+1);
  *pit=pitch;
  *sym=symmetry;
  *coresize=core_size;
  fclose(cist);
  free(buff);
  return Colormap;
}


void PrintProlog(FILE *psat,double pitch,int core_size,int symmetry,imagecontrol *IC,int *dbb,const char layout)
{
  double unit;  
  int BBx,BBy;
  /* define the unit size (or pitch ...) */
  /* 28.346==1cm */
  unit=28.346*0.5*pitch;
  if(layout=='h'){
    if(symmetry==6){    
      BBx=unit*(core_size+2*1.2-1)+2+dbb[2];
      BBy=unit/sqrt(3)*((core_size-1)/2*3+4.2)+2+dbb[3];
    }
    else{
      BBx=2*unit*(core_size+0.2)+2+dbb[2];
      BBy=unit/sqrt(3)*((core_size-1)*3+4.2)+2+dbb[3];
    }  
  }
  else{
    BBx=2*unit*(core_size+0.1)+2+dbb[2];
    BBy=2*unit*(core_size+0.1)+2+dbb[3];
  }
  fprintf(psat,"%s\n","%!PS-Adobe-2.0 EPSF-2.0");
  fprintf(psat,"%s\n","%%Orientation: Portrait");
  fprintf(psat,"%s  %d %d ","%%BoundingBox:",dbb[0],dbb[1]);
  if(IC->img_type=='c')
    fprintf(psat,"%d %d\n",BBx,BBy);
  else
    fprintf(psat,"%d %d\n",BBy,BBx);
  fprintf(psat,"%s\n","%%Pages: 0");
  fprintf(psat,"%s\n","%%BeginSetup");
  fprintf(psat,"%s\n","%%EndSetup");
  fprintf(psat,"%s\n","%%EndComments");
  fprintf(psat,"%s\n","");
  fprintf(psat,"%s\n","/$F2psDict 300 dict def");
  fprintf(psat,"%s\n","$F2psDict begin");
  fprintf(psat,"%s\n","$F2psDict /mtrx matrix put");
  fprintf(psat,"%s\n","/col-1 {0 setgray} bind def");
  fprintf(psat,"%s\n","/col0 {0.000 0.000 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col1 {1.000 1.000 1.000 srgb} bind def");
  fprintf(psat,"%s\n","/col2 {0.000 0.000 0.950 srgb} bind def");
  fprintf(psat,"%s\n","/col3 {0.000 0.450 0.950 srgb} bind def");
  fprintf(psat,"%s\n","/col4 {0.000 0.750 0.800 srgb} bind def");
  fprintf(psat,"%s\n","/col5 {0.000 0.850 0.600 srgb} bind def");
  fprintf(psat,"%s\n","/col6 {0.000 0.950 0.400 srgb} bind def");
  fprintf(psat,"%s\n","/col7 {0.000 0.950 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col8 {0.700 0.950 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col9 {0.950 0.950 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col10 {0.990 0.750 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col11 {0.990 0.450 0.000 srgb} bind def");
  fprintf(psat,"%s\n","/col12 {0.950 0.000 0.000 srgb} bind def");
  fprintf(psat,"%s\n","");
  fprintf(psat,"%s\n","end");
  fprintf(psat,"%s\n","save");
  fprintf(psat,"%s\n","/srgb {setrgbcolor} bind def");
  fprintf(psat,"%s\n","/s {stroke} bind def");
  fprintf(psat,"%s\n","/f {fill} bind def");
  fprintf(psat,"%s\n","/n {newpath} bind def");
  fprintf(psat,"%s\n","/m {moveto} bind def");
  fprintf(psat,"%s\n","/rm {rmoveto} bind def");
  fprintf(psat,"%s\n","/rs {restore} bind def");
  fprintf(psat,"%s","/ut {");
  fprintf(psat,"%.3f",unit); 
  fprintf(psat,"%s\n"," mul} bind def");
  /* define the hex, pitch=2cm */
  fprintf(psat,"%s\n","/hex {1 ut -0.577 ut rmoveto 0 1.155 ut rlineto ");
  fprintf(psat,"%s\n","      -1 ut 0.577 ut rlineto -1 ut -0.577 ut rlineto ");
  fprintf(psat,"%s\n","      0 -1.155 ut rlineto 1 ut -0.577 ut rlineto"); 
  fprintf(psat,"%s\n","      closepath -1 ut 0.577 ut rmoveto } bind def");
  fprintf(psat,"%s\n","/sq {-1 ut -1 ut rmoveto 2 ut 0 ut rlineto ");
  fprintf(psat,"%s\n","      0 ut 2 ut rlineto -2 ut 0 ut rlineto ");
  fprintf(psat,"%s\n","      closepath 1 ut 1 ut rmoveto } bind def");
  if(IC->img_type=='c')
    fprintf(psat,"%s\n","/box {-1 ut -0.5 ut rmoveto 2 ut 0 rlineto 0 1 ut rlineto -2 ut 0 rlineto");
  else
    fprintf(psat,"%s\n","/box {-1 ut -0.5 ut rmoveto 1.4 ut 0 rlineto 0 1 ut rlineto -1.4 ut 0 rlineto");
  fprintf(psat,"%s\n","closepath 1 ut 0.5 ut rmoveto } bind def"); 
  fprintf(psat,"%s\n","");
  fprintf(psat,"%s\n","/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def");
  fprintf(psat,"%s\n","/$F2psEnd {$F2psEnteredState restore end} def");
  fprintf(psat,"%s\n","%%EndProlog");
}


/*  reads buff_size chars from a line terminated by eol */
/*  returns number chars read (0 if the line is empty) */
/*  returns -1 if 0 characters are read and EOF is reached */
int GetLine(FILE *cist,char *buff,int buff_size,const char eol)
{
  int i=0,flag=0;
  char c;
  while((c=fgetc(cist))!=EOF && c!=eol){
    if(i<buff_size-1){
      buff[i++]=c;
      if(c!=' ' && c!='\t') flag=1;
    }
  }
  buff[i]='\0';
  if(i==0 && c==EOF) return -1;
  return flag*i;
}


/*  ************************************************************** */
/*  split read line to integer values and fill them into the field */
/*  returns number of read values */
/*  ************************************************************** */
int Line2IntField(char *buff,int *field,int fieldsize)
{
  int pos,n;  
  char *substr;  
  n=0;while(buff[n]==' ' || buff[n]=='\t') ++n;
  substr=buff+n;
  for(pos=0;pos<fieldsize;++pos){
    if(*substr!='\0') field[pos]=atoi(substr);
    else return(pos);    
    while(*substr!=' ' && *substr!='\t' && *substr!='\0') ++substr;
    while(*substr==' ' || *substr=='\t') ++substr;
  }
  return(pos);  
}


void PrintCore6(FILE *psat,int n_anneaux,int *Coremap,char *vfile,char *colorfile,double *dxy,double *printcolor,
		imagecontrol *IC,vf_header *vfhp)
{
  double x,y,x1,y1,sq3,dx,dy,**res;
  char buffont[200],buff[20];  
  int core_size,alloc_size,ii,jj,i,j,count,count_color,k,*precision,colcolumn,color,*colors,ncolors;
  int ncols,nvals;
  if(printcolor) colcolumn=printcolor[2];
  sq3=sqrt(3);
  core_size=2*n_anneaux+1;
  alloc_size=core_size+1;
  if(vfile)
    res=(double**)ReadVfile(vfile,&ncols,&nvals,&precision,vfhp);
  if(colorfile)
    colors=(int*)ReadColorFile(colorfile,&ncolors);
  if(0 ){ /* debug */
    printf("haha %d %d\n",ncols,nvals);
    for(i=0;i<nvals;++i){
      printf("%d",i);      
      for(j=0;j<ncols;++j){
	printf(" %f",res[i][j]);
      }
      printf("\n");
    }
    exit(1);
  }
  fprintf(psat,"%s\n","$F2psBegin");
  fprintf(psat,"%s\n","0.07 ut setlinewidth");
  fprintf(psat,"%s\n","");
  fprintf(psat,"%s\n","/Times-Roman findfont");
  if(vfile){
    if(ncols==1)
      sprintf(buffont,"%f ut scalefont setfont",0.7);
    else if(ncols==3)
      sprintf(buffont,"%f ut scalefont setfont",1.6/(ncols));
    else
      sprintf(buffont,"%f ut scalefont setfont",1.2/(ncols));
    fprintf(psat,"%s\n",buffont);
  }
  else
    fprintf(psat,"%s\n","1 ut scalefont setfont");
  fprintf(psat,"%s\n","");
  count=0;
  count_color=0;
  /* print the color scale, if color */
  if(printcolor){
    if(IC->img_type=='c'){
      y=1.0/sqrt(3)*(n_anneaux*3+4.2)-0.25;
      x=1.0;
      dx=0;
      dy=-1;
    }
    else{
      y=2*(n_anneaux+1.2)-0.25;
      x=1.0/sqrt(3)*((n_anneaux)*3+4.2)+2-11.5*1.4;
      dx=1.4;
      dy=0;
    }
    sprintf(buff,"col%d",2);
    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
    fprintf(psat,"box\ngsave\n");
    fprintf(psat,"%s fill\ngrestore\ncol0\nstroke\n",buff);
    for(i=0;i<10;++i){
      if(IC->img_type=='c')
	y+=dy;
      else
	x+=dx;
      sprintf(buff,"col%d",3+i);
      fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
      fprintf(psat,"box\ngsave\n");
      fprintf(psat,"%s fill\ngrestore\ncol0\nstroke\n",buff);
    }
    if(IC->img_type=='c'){
      y=1.0/sqrt(3)*(n_anneaux*3+4.2)-0.25-0.5-0.2;
      x=2.5;
    }
    else{
      y=2*(n_anneaux+1.2)-0.25-1.2;
      x=1.0/sqrt(3)*((n_anneaux)*3+4.2)+2-11.5*1.4-0.7;
    }
    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
    fprintf(psat,"(%6.2f) show\n",printcolor[0]);
    for(i=1;i<10;++i){
      if(IC->img_type=='c'){
	y+=dy;
	dx=printcolor[0]+(printcolor[1]-printcolor[0])*i/9;
	fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	fprintf(psat,"(%6.2f) show\n",dx);
      }
      else{
	x+=dx;
	dy=printcolor[0]+(printcolor[1]-printcolor[0])*i/9;
	fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	fprintf(psat,"(%6.2f) show\n",dy);
      }
    }      
  }
  for(ii=n_anneaux;ii>=0;--ii){
    i=n_anneaux-ii;
    for(jj=n_anneaux-i;jj>=0;--jj){
      j=n_anneaux-i-jj;
      if(Coremap[ii*alloc_size+jj]){
	x=2*j+i+1.1;
	y=sq3*i+2.0/sq3+0.1;
	if(IC->img_type!='c'){
	  x-=1.1;
	  y-=(2.0/sq3+0.1);
	  x1=sq3/2*x-0.5*y;
	  y1=0.5*x+sq3/2*y;
	  x=x1+(2.0/sq3+0.1);
	  y=y1+1.1;
	}
	x+=dxy[2];
	y+=dxy[3];
	if(!printcolor){
	  if(IC->img_type=='c'){
	    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	    fprintf(psat,"hex\n");
	  }
	  else
	    fprintf(psat,"newpath %.2f ut %.2f ut 1.0 ut 0 360 arc closepath stroke\n",x,y);
	}
	else{
	  if(res[count][colcolumn]==9.9E99 || Coremap[i*alloc_size+j]<0) color=1;
	  else{
	    color=(res[count][colcolumn]-printcolor[0])/(printcolor[1]-printcolor[0])*9+3;
	    if(color<2) color=2;
	    if(color>12) color=12;
	  }
	  sprintf(buff,"col%d",color);
	  if(IC->img_type=='c'){
	    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	    fprintf(psat,"hex\ngsave\n");
	  }
	  else{
	    fprintf(psat,"newpath %.2f ut %.2f ut 1 ut 0 360 arc closepath\n",x,y);
	    fprintf(psat,"gsave\n");
	  }
	  fprintf(psat,"%s fill\ngrestore\n",buff);
	  if(!colorfile)
	    fprintf(psat,"col0\n");
	  else{
	    fprintf(psat,"col%d\n",colors[count_color]);
	    ++count_color;
	  }
	  fprintf(psat,"stroke\n");
	}
	if(vfile){
	  fprintf(psat,"col0\n");
	  dx=-0.85+dxy[0];
	  for(k=0;k<ncols;++k){
	    if(k==ncols-1)
	      fprintf(psat,"%s\n","/Times-Roman-Bold findfont");
	    else
	      fprintf(psat,"%s\n","/Times-Roman findfont");
	    fprintf(psat,"%s\n",buffont);
	    dy=0.15-1.2*k/ncols+dxy[1];
	    if(ncols==3)
	      dy=0.3-1.4*k/ncols+dxy[1];
	    fprintf(psat,"%.2f ut %.2f ut m\n",x+dx,y+dy);
	    if(res[count][k]==0/*  || k==0 */)
	      fprintf(psat,"(%6.0f) show\n",res[count][k]);
	    else if(res[count][k]==9.9E99) /* NaN  */
	      fprintf(psat,"(%6s) show\n","NaN");
	    else{
	      if(precision[k]==0)
		fprintf(psat,"(%6.0f) show\n",res[count][k]);
	      else if(precision[k]==1)
		fprintf(psat,"(%6.1f) show\n",res[count][k]);
	      else if(precision[k]==2)
		fprintf(psat,"(%6.2f) show\n",res[count][k]);
	      else 
		fprintf(psat,"(%6.3f) show\n",res[count][k]);
	    }
	  }	  
	}
	else{
	  dx=0.3;
	  dy=0.3;
	  if(count>10) dx=0.5;
	  fprintf(psat,"%.2f ut %.2f ut m\n",x-dx,y-dy);
	  fprintf(psat,"(%d) show\n",count);
	}
	++count;	
      }
    }
  }
  fprintf(psat,"%s\n","col0 stroke");
  fprintf(psat,"%s\n","$F2psEnd");
  fprintf(psat,"%s\n","rs");
}


void PrintCore1(FILE *psat,int core_size,int *Coremap,char *vfile,char *colorfile,double *dxy,double *printcolor,
		imagecontrol *IC,vf_header *vfhp,const char layout)
{
  double x,y,x1,y1,sq3,dx,dy,**res;
  int alloc_size,i,j,count,count_color,k,ncols,*precision,color,colcolumn,*colors,ncolors;
  char buffont[200];  
  char buff[20];
  int nvals;
  if(printcolor) colcolumn=printcolor[2]+0.5;
  sq3=sqrt(3);
  alloc_size=core_size+1;
  if(vfile)
    res=(double**)ReadVfile(vfile,&ncols,&nvals,&precision,vfhp);
  if(colorfile)
    colors=(int*)ReadColorFile(colorfile,&ncolors);
  fprintf(psat,"%s\n","$F2psBegin");
  fprintf(psat,"%s\n","0.05 ut setlinewidth");
  fprintf(psat,"%s\n","");
  fprintf(psat,"%s\n","/Times-Roman findfont");
  if(vfile){
    if(ncols==1)
      sprintf(buffont,"%f ut scalefont setfont",0.7);
    else if(ncols==3)
      sprintf(buffont,"%f ut scalefont setfont",1.6/(ncols));
    else
      sprintf(buffont,"%f ut scalefont setfont",1.2/(ncols));
    fprintf(psat,"%s\n",buffont);
  }
  else
    fprintf(psat,"%s\n","1 ut scalefont setfont");
  fprintf(psat,"%s\n","");
  count=0;
  count_color=0;
  /* print the color scale, if color */
  if(printcolor){
    if(IC->img_type=='c'){
      y=1.0/sqrt(3)*((core_size-1)*3+4.2)-0.25;
      x=1.0;
      dx=0;
      dy=-1;
    }
    else{
      y=2*(core_size+0.2)-0.25;
      x=1.0;
      dx=1.4;
      dy=0;
    }
    sprintf(buff,"col%d",2);
    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
    fprintf(psat,"box\ngsave\n");
    fprintf(psat,"%s fill\ngrestore\ncol0\nstroke\n",buff);
    for(i=0;i<10;++i){
      if(IC->img_type=='c')
	y+=dy;
      else
	x+=dx;
      sprintf(buff,"col%d",3+i);
      fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
      fprintf(psat,"box\ngsave\n");
      fprintf(psat,"%s fill\ngrestore\ncol0\nstroke\n",buff);
    }
    if(IC->img_type=='c'){
      y=1.0/sqrt(3)*((core_size-1)*3+4.2)-0.25-0.5-0.2;
      x=2.5;
    }
    else{
      y=2*(core_size+0.2)-0.25-1.2;
      x=0.7;
    }
    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
    fprintf(psat,"(%6.2f) show\n",printcolor[0]);
    for(i=1;i<10;++i){
      if(IC->img_type=='c'){
	y+=dy;
	dx=printcolor[0]+(printcolor[1]-printcolor[0])*i/9;
	fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	fprintf(psat,"(%6.2f) show\n",dx);
      }
      else{
	x+=dx;
	dy=printcolor[0]+(printcolor[1]-printcolor[0])*i/9;
	fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	fprintf(psat,"(%6.2f) show\n",dy);
      }
    }      
  }
  for(i=0;i<core_size;++i){
    for(j=0;j<core_size;++j){
      if(Coremap[i*alloc_size+j]){
	if(layout=='s'){
	  x=2*j+1.1;
	  y=2*i+1.1;
	}
	else{
	  x=2*j-i+(core_size-1)/2+1.1;
	  y=sq3*i+2.0/sq3+0.1;
	}
	if(IC->img_type!='c' && layout!='s'){
	  x-=(core_size+0.1);
	  y-=(sq3*(core_size-1)/2+2.0/sq3+0.1);
	  x1=sq3/2*x-0.5*y;
	  y1=0.5*x+sq3/2*y;
	  x=x1+(sq3*(core_size-1)/2+2.0/sq3+0.1);
	  y=y1+(core_size+0.1);
	}
	x+=dxy[2];
	y+=dxy[3];
	if(!printcolor){
	  if(IC->img_type=='c'){
	    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	    if(layout=='h')
	      fprintf(psat,"hex\n");
	    else
	      fprintf(psat,"sq\n");
	  }
	  else
	    fprintf(psat,"newpath %.2f ut %.2f ut 1.0 ut 0 360 arc closepath stroke\n",x,y);
	}
	else{
	  if(res[count][colcolumn]==9.9E99 || Coremap[i*alloc_size+j]<0) color=1;
	  else{
	    color=(res[count][colcolumn]-printcolor[0])/(printcolor[1]-printcolor[0])*9+3;
	    if(color<2) color=2;
	    if(color>12) color=12;
	  }
	  sprintf(buff,"col%d",color);
	  if(IC->img_type=='c'){
	    fprintf(psat,"%.2f ut %.2f ut m\n",x,y);
	    if(layout=='h')
	      fprintf(psat,"hex\ngsave\n");
	    else
	      fprintf(psat,"sq\ngsave\n");
	  }
	  else{
	    fprintf(psat,"newpath %.2f ut %.2f ut 1 ut 0 360 arc closepath\n",x,y);
	    fprintf(psat,"gsave\n");
	  }
	  fprintf(psat,"%s fill\ngrestore\n",buff);
	  if(!colorfile)
	    fprintf(psat,"col0\n");
	  else{
	    fprintf(psat,"col%d\n",colors[count_color]);
	    ++count_color;
	  }
	  fprintf(psat,"stroke\n");
	}
	if(vfile && Coremap[i*alloc_size+j]>0){
	  fprintf(psat,"col0\n");
	  dx=-0.85+dxy[0];
	  for(k=0;k<ncols;++k){
	    if(k==ncols-1)
	      fprintf(psat,"%s\n","/Times-Roman-Bold findfont");
	    else
	      fprintf(psat,"%s\n","/Times-Roman findfont");
	    fprintf(psat,"%s\n",buffont);
	    dy=0.15-1.2*k/ncols+dxy[1];
	    if(ncols==3)
	      dy=0.3-1.4*k/ncols+dxy[1];
	    fprintf(psat,"%.2f ut %.2f ut m\n",x+dx,y+dy);
	    if(res[count][k]==0/*  || k==0 */)
	      fprintf(psat,"(%6.0f) show\n",res[count][k]);
	    else if(res[count][k]==9.9E99) /* NaN  */
	      if(IC->print_nan=='y')
		fprintf(psat,"(%6s) show\n","NaN");
	      else
		fprintf(psat,"(%6s) show\n","   ");
	    else{
	      if(precision[k]==0)
		fprintf(psat,"(%6.0f) show\n",res[count][k]);
	      else if(precision[k]==1)
		fprintf(psat,"(%6.1f) show\n",res[count][k]);
	      else if(precision[k]==2)
		fprintf(psat,"(%6.2f) show\n",res[count][k]);
	      else 
		fprintf(psat,"(%6.3f) show\n",res[count][k]);
	    }
	  }	  
	}
	else if(0){	  
	  dx=0.3;
	  dy=0.3;
	  if(count>10) dx=0.5;
	  fprintf(psat,"%.2f ut %.2f ut m\n",x-dx,y-dy);
	  fprintf(psat,"(%d) show\n",count);
	}	
	if(Coremap[i*alloc_size+j]>0) ++count;	
      }
    }
  }
  fprintf(psat,"%s\n","col0 stroke");
  fprintf(psat,"%s\n","$F2psEnd");
  fprintf(psat,"%s\n","rs");
}








void PrintSampleConffile()
{
  printf("%s\n","#Sample conffile:");
  printf("%s\n","#output file");
  printf("%s\n","mycore.eps");
  printf("%s\n","#pitch");
  printf("%s\n","2");
  printf("%s\n","#layout (h)ex/(s)quare");
  printf("%s\n","h");
  printf("%s\n","#nanneaux(hex)/size(square)");
  printf("%s\n","7");
  printf("%s\n","#symmetry(1,6, for square only 1!)");
  printf("%s\n","6");
  printf("%s\n","#couleurs");
  printf("%s\n","0 3 1 1 2 1 1 1");
  printf("%s\n"," 3 2 2 1 1 2 0");
  printf("%s\n","  3 1 1 1 1 0");
  printf("%s\n","   3 2 1 1 0");
  printf("%s\n","    3 1 2 0");
  printf("%s\n","     3 2 0");
  printf("%s\n","      3 0");
  printf("%s\n","       0");
  printf("%s\n","");
  printf("%s\n","#dxy (text dx and dy shift - fine tune");
  printf("%s\n","0.0 0.0");
  printf("%s\n","");
}


void* ReadVfile(char *vfile,int *nc,int *nv,int **precision,vf_header *vfhp)
{
  FILE *cist;
  double **res;
  int ncols,nvals,i,j;
  
  cist=fopen(vfile,"r");
  if(!cist){
    fprintf(stderr,"ReadVfile: cannot open file %s. Aborting.\n",vfile);
    exit(1);
  }
  if(!vfhp){
    fscanf(cist,"%d %d",&ncols,&nvals);
    *nc=ncols;
    *nv=nvals;
    *precision=(int*)malloc(sizeof(int)**nc);
    for(i=0;i<*nc;++i)   fscanf(cist,"%d",(*precision)+i);
  }
  else{
    *nc=vfhp->ncols;
    *nv=vfhp->nvals;
    ncols=*nc;
    nvals=*nv;
    *precision=vfhp->precision;
  }
  res=(double**)malloc(nvals*sizeof(double*));
  for(i=0;i<nvals;++i)
    res[i]=(double*)malloc(ncols*sizeof(double));
  for(i=0;i<nvals;++i){
    fscanf(cist,"%d",&j);
    for(j=0;j<ncols;++j)
      fscanf(cist,"%lf",&(res[i][j]));
  }
  fclose(cist);
  return res;
}

void* ReadColorFile(char *colorfile,int *nv)
{
  FILE *cist;
  int *res;
  int nvals,i,j;
  
  cist=fopen(colorfile,"r");
  if(!cist){
    fprintf(stderr,"ReadColorFile: cannot open file %s. Aborting.\n",colorfile);
    exit(1);
  }
  fscanf(cist,"%d",&nvals);
  *nv=nvals;
  res=(int*)malloc(nvals*sizeof(int));
  for(i=0;i<nvals;++i){
    fscanf(cist,"%d%d",&j,&(res[i]));
  }
  fclose(cist);
  return res;
}


