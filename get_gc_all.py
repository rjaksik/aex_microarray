#!/usr/bin/python2

# get_gc_ctrl  v1.0  02.04.13

import os, re, math, sys, getopt 
import rpy2.robjects as robjects
import numpy
from numpy import *
import scipy.stats
import time
import datetime


def usage():
   print """
   get_gc_ctrl.py directory trans_gc trans_gc_colNr cdf_name tmreport_file output_file 

   DESCRIPTION
     calculates difference in probeset signal changes between two samples and their
     average slope difference + mean intensity FC, the statistics are averaged between probesest of similar GC

   INPUT:     
     directory - directory containing individual experiment subdirectories 
     cdf_name - name of the CDF file to be used ex: hgu133ahsrefseq for Dai custom RefSeq based CDF
     tmreport_file - file with tm plot regression parameters

   OUTPUT:
     output_file - file with results in a format: Sample1<tab>Sample2<tab>stats....

   FLAGS:
     --rma - analyze entire experiment normalized with RMA (default)
     --mas - analyze entire experiment normalized with MAS5
     --gcrma - analyze entire experiment normalized with GCRMA
     --plier - analyze entire experiment normalized with PLIER
     --farms - analyze entire experiment normalized with FARMS
     --mbei - analyze entire experiment normalized with MBEI
     --dat=N - analyze only a specific dataset number N
     -u - continue analysis based on results file
     -m - use median array as reference

    """

def main(argv):
  try:
    #processing options and parameters
    opts, args = getopt.gnu_getopt(argv, "u", ['rma','mas','gcrma','plier','farms','mbei','dat=','med'])
  except:
    #incorrect parameters
    print "Error: Incorrect parameters"
    usage()
    sys.exit()

  if (len(args)<6) or (len(args)>6) :
    #application requires 3 parameters
    usage()
    sys.exit()
  else:

    #----starting the script----
    start = time.clock()

    expdir=args[0];
    gc_file_name=args[1];
    gc_file_col=int(args[2])-1;
    cdfname=args[3];
    tmrep_file_name=args[4];
    outfile=args[5];

    if expdir[len(expdir)-1]=='/':
       expdir=expdir[:-1];


    # GC% intervals
    gc_inter=[0,35,40,45,50,55,60,65,100];

    # process arguments
    Nset=0;
    procMode=0;
    Nmet='RMA';
    update=False;
    medref=False;
    for opt,oarg in opts:
       opt=opt.replace('--','-')
       if opt in ('-dat'):
          Nset=int(oarg);
       if opt in ('-rma'):
          procMode=0;
          Nmet='RMA';
       if opt in ('-gcrma'):
          procMode=1;
          Nmet='GCRMA';
       if opt in ('-plier'):
          procMode=2;
          Nmet='PLIER';
       if opt in ('-farms'):
          procMode=3;
          Nmet='FARMS';
       if opt in ('-mas'):
          procMode=4;
          Nmet='MAS5';
       if opt in ('-mbei'):
          procMode=5;
          Nmet='MBEI';
       if opt in ('-u'):
          update=True;
       if opt in ('-med'):
          medref=True;
    print Nmet;

    #assign an analysis summary file
    global report_file
    report_file=expdir+'_'+Nmet+'.report';


    #read the tmreport file
    tmrep_slope_db={};
    tmrep_mean_db={};
    tmrep_mean_exp_db={};
    tmrep_mean_exp_db_ct={};
    tmrep_file=open(tmrep_file_name,'r');
    for line in tmrep_file.readlines():
       line=line.strip();
       pline=line.split('\t');
       if not pline[0][0]=='#':
          tmrep_slope_db[pline[0]+':'+pline[1]]=float(pline[4]);
          tmrep_mean_db[pline[0]+':'+pline[1]]=float(pline[8]);
          if pline[0] in tmrep_mean_exp_db:
             tmrep_mean_exp_db[pline[0]]+=float(pline[8]);
             tmrep_mean_exp_db_ct[pline[0]]+=1;
          else:
             tmrep_mean_exp_db[pline[0]]=float(pline[8]);
             tmrep_mean_exp_db_ct[pline[0]]=1;

    tmrep_file.close()
    
    #average the mean_exp_db
    for exp in tmrep_mean_exp_db:
       tmrep_mean_exp_db[exp]=tmrep_mean_exp_db[exp]/tmrep_mean_exp_db_ct[exp];

    #read the gc file
    gc_db={};
    gc_file=open(gc_file_name,'r');
    for line in gc_file.readlines():
       line=line.strip();
       pline=line.split('\t');
       if not pline[0][0]=='#':
          ID=pline[0].strip('_at');
          gc_db[ID]=float(pline[gc_file_col]);
    gc_file.close()
    
    #list experiment that where already analyzed for the update mode
    analyzed_dat=[]
    if update:
       out_file=open(outfile,'r');
       for line in out_file.readlines():
          line=line.strip();
          pline=line.split();
          exper=pline[1]; #sec column only for lines that contain the summary
          if not (exper in analyzed_dat):
             analyzed_dat.append(exper);

    #list all experiments
    os.system('ls %s >explist_%s.tmp' %(expdir,expdir))

    #read the experiment list
    TotalSets=0;
    list_file=open('explist_'+expdir+'.tmp','r')
    exptab=[];
    for exper in list_file.readlines():
       exper=exper.strip();
       if os.path.isdir(expdir+'/'+exper):
          if exper<>'back':
             exptab.append(exper);
             TotalSets+=1;
    list_file.close()
    #os.system('rm explist_'+expdir+'.tmp')

    #analyze each experiment
    Nstart=0;
    Nstop=TotalSets;


    # change the analysis range if only one experiment will be processed
    if Nset>0:
       if Nset-1>=TotalSets: #id above the allowed range
          Nstart=0;
          Nstop=0;
          print 'invalid dataset id: '+str(Nset)
       else:                 #correct id
          Nstart=Nset-1;
          Nstop=Nset;

    #create header if neccesery
    if not update:
       tmpstr='#Experiment\tMethod\tSample1\tSample2\tSlopeDiff\tGCvsLFC_rho\tGCvsLFC_p-val';
       if Nset<=1:
          os.system('echo "'+tmpstr+'" > '+outfile);

    #analyze datasets
    analyzed_ct=0;
    failed_ct=0;
    first=True;
    list_file=open('explist_'+expdir+'.tmp','r')
    for k in range(Nstart,Nstop):
       exper=exptab[k]; # analyze this experiment
       analyze=True;
       if update:
          if exper in analyzed_dat:
             analyze=False;
             print '  skipping dataset: %s nr: %s' %(exper,k+1);
       
       if analyze:
           if os.path.isdir(expdir+'/'+exper):

              print '=== Analyzing: %s nr: %s ===' %(exper,k+1);

              #extract data from cel files using R's parsing system
              code='library(Biobase,verbose=F,warn.conflicts=F);'
              code+='library(BiocGenerics,verbose=F,warn.conflicts=F);'
              code+='library(AnnotationDbi,verbose=F,warn.conflicts=F,quietly=T);'
              code+='library(affy,verbose=F,warn.conflicts=F);'                                
              code+='LFS=list.celfiles(path="'+expdir+'/'+exper+'", full.names=TRUE);'
              code+='N=length(LFS);'
              result=robjects.r(code)
                            
              Nsamp=robjects.r['N']
              Nsamp=Nsamp[0]                                        

              ############################# multi - sample #############################
              if procMode>=0:                    
                    success=False;
                    try:
                       code='raw<-ReadAffy(filenames=LFS,cdfname="'+cdfname+'");' 
                       if procMode==0: #RMA
                          code+='normdata <- rma(raw,verbose=F,);';
                       if procMode==1: #GC-RMA
                          code+='library(gcrma,verbose=F,warn.conflicts=F);';
                          code+='normdata <- gcrma(raw,verbose=F,);'
                       if procMode==2: #PLIER
                          code+='library(plier,verbose=F,warn.conflicts=F);';
                          code+='normdata <- justPlier(raw);';
                       if procMode==3: #FARMS
                          code+='library(farms,verbose=F,warn.conflicts=F);';
                          code+='normdata <- q.farms(raw,verbose=F,);'
                       if procMode==4: #FARMS
                          code+='library(preprocessCore,verbose=F,warn.conflicts=F);';
                          code+='normdata <- mas5(raw,verbose=F);'
                       if procMode==5: #MBEI
                          code+='normdata <- expresso(raw, normalize.method="invariantset",bg.correct=FALSE,pmcorrect.method="pmonly",summary.method="liwong",verbose = F);'
                       result=robjects.r(code)
                       code='write.table(exprs(normdata), file="'+expdir+'/'+exper+'.setdata", sep="\t",  quote=F);'                 
                       result=robjects.r(code)

                       code='G=dim(exprs(normdata))[1];'
                       code+='N=length(LFS);'
                       result=robjects.r(code)
                       Nset=robjects.r['G'];
                       Nset=Nset[0];
                       Nsamp=robjects.r['N']
                       Nsamp=Nsamp[0];

                       code+='rm(list=ls());'
                       result=robjects.r(code)
                       success=True;
                    except:                      
                       print "Failed to analyze: "+exper
                       update_report('Failed (R)',exper,'');
                       failed_ct+=1

                    if success==True:       
                       #read the data into dictionary
                       fline=True
                       Ninter=len(gc_inter)-1;
                       expdata=zeros([Nset,Nsamp+1],float);
                       data_file=open(expdir+'/'+exper+'.setdata')
                       z=0;
                       for line in data_file.readlines():
                          line=line.strip()
                          pline=line.split('\t')
                          if fline:
                             headers=pline
                          else:      
                             ID=pline[0].strip('_at');
                             if ID in gc_db and ID.find('AFFX')<0:
                                expdata[z,0]=gc_db[ID];
                                expdata[z,1:Nsamp+1]=pline[1:len(pline)];
                                z+=1;

                          fline=False                 
                       data_file.close()                       
                       os.system('rm '+expdir+'/'+exper+'.setdata')

                       final_corr=zeros([Nsamp*Nsamp,2],float)
                       final_corr_ct=0;
                       results='';
                       success=True;

                       maxmdiff=0;
                       maxslope=0;
                       
                       #generating report - vs ref array
                       if medref:
                         refarray=numpy.median(expdata[:,2:],1)
                         for i in range (0,Nsamp):
                            expnum=zeros(Ninter,int);
                            tempexp=zeros(Ninter,float);

                            fc=expdata[:,i+1]-refarray;

                            idtm1=exper+':'+headers[i];
                            if idtm1 in tmrep_slope_db:
                               mdiff=log2(tmrep_mean_db[idtm1]/tmrep_mean_exp_db[exper]);

                               if maxmdiff<mdiff:
                                  maxmdiff=mdiff;
                            
                               for nj in range(0,z):
                                  GCval=expdata[nj,0];
                                  for ni in range(0,Ninter):
                                     if GCval>gc_inter[ni] and GCval<=gc_inter[ni+1]:
                                        expnum[ni]+=1;
                                        tempexp[ni]+=fc[nj];
   
                               tmpstr='';
                               for ni in range(0,Ninter):
                                  tmpstr+='\t%.3f' %(tempexp[ni]/expnum[ni]);

                               os.system('echo "%s\t%s\t%s\t%s\t%s%s" >> %s' %(expdir,exper,Nmet,headers[i],mdiff,tmpstr,outfile));
                            else:
                               success=False;
                               print "Missing Tm data: "+idtm1;


                       #generating report - pairs of samples
                       if not medref:
                         for j in range (0,Nsamp):
                           for i in range (0,Nsamp):
                             if i>j:
                                expnum=zeros(Ninter,int);
                                tempexp=zeros(Ninter,float);

                                idtm1=exper+':'+headers[i];
                                idtm2=exper+':'+headers[j];
                                if idtm1 in tmrep_slope_db:
                                   if idtm2 in tmrep_slope_db:
                                      fc=expdata[:,i+1]-expdata[:,j+1];                                
                                      sdiff=tmrep_slope_db[idtm1] - tmrep_slope_db[idtm2]
                                      mdiff=tmrep_mean_db[idtm1]/tmrep_mean_db[idtm2]
                                      #if mdiff<1:
                                      #   fc=fc * -1;
                                      #   fc=expdata[:,j+1]-expdata[:,i+1];                                
                                      #   sdiff=tmrep_slope_db[idtm2] - tmrep_slope_db[idtm1]
                                      #   mdiff=tmrep_mean_db[idtm2]/tmrep_mean_db[idtm1];
                                      
                                      if maxmdiff<mdiff:
                                         maxmdiff=mdiff;

                                      if maxslope<abs(sdiff):
                                         maxslope=abs(sdiff);
                                         
                                      for nj in range(0,z):
                                         GCval=expdata[nj,0];
                                         for ni in range(0,Ninter):
                                            if GCval>gc_inter[ni] and GCval<=gc_inter[ni+1]:
                                               expnum[ni]+=1;
                                               tempexp[ni]+=fc[nj];

                                      tmpstr='';
                                      for ni in range(0,Ninter):
                                         tmpstr+='\t%.3f' %(tempexp[ni]/expnum[ni]);
                                      
                                      os.system('echo "%s\t%s\t%s\t%s\t%s\t%s\t%s%s" >> %s' %(expdir,exper,Nmet,headers[i],headers[j],sdiff,mdiff,tmpstr,outfile));
                                   else:
                                      success=False;
                                      print "Missing Tm data: "+idtm2;                                      
                                else:
                                   success=False;
                                   print "Missing Tm data: "+idtm1;

                       os.system('echo "MAXDIFF\t%s\t%s" >> %s' %(maxmdiff,maxslope,outfile));

                       if success:
                          #generate analysis log entry
                          update_report('Success',exper+' ['+str(Nsamp)+']','');
                          analyzed_ct+=1
                       else:
                          update_report('Failed',exper+' ['+str(Nsamp)+']','');
                          failed_ct+=1
              first=False;

    list_file.close()
    
    #analysis time
    end = time.clock()
    calc_time = end - start;
    if calc_time>60:
       print 'Code time: %.d min %.d sec' % (round(calc_time/60),calc_time % 60)
    else:
       print 'Code time: %.d sec' % (calc_time)

    #generate summary report if more than one dataset was analyzed
    if Nset==0:
       rep_str='Batch success - analyzed: %s   failed: %s   total time: %s min' %(analyzed_ct,failed_ct,round(calc_time/60))
       update_report(rep_str,'','');

       

def update_report(status,exper,sampName):
   #generate data for the report file
   now = datetime.datetime.now();
   os.system('echo "%s\tCTRL\t%s\t%s\t%s" >> %s' %(now.strftime("%Y-%m-%d %H:%M"),status,exper,sampName,report_file))


if __name__ == "__main__":
    main(sys.argv[1:])


