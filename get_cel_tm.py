#!/usr/bin/python2

# get_cel_tm  v1.6.1  31.08.12 

import os, re, math, sys, getopt 
import rpy2.robjects as robjects
import numpy
from numpy import *
import time
import datetime
import subprocess

def usage():
   print """
   get_cel_tm.py directory pgc_file cdf_name output_file 

   DESCRIPTION:
     calculates Tm plot values for each sample in directory based on specified pgc_file

   INPUT:     
     directory - directory containing individual experiment subdirectories 
     pgc2_file - file containing GC nucleotide counts for each probe - format: ProbesetID<tab>ProbeID<tab>ProbeSequence<tab>GC<tab>ProbeLoc<tab>AnStats
     cdf_name - name of the CDF file to be used ex: hgu133ahsrefseq for Dai custom RefSeq based CDF

   OUTPUT:
     output_file - file where first column includes the sample name and remaining 26 expression level % for each Tm

   FLAGS:
     basic:
     -u - continue analysis based on results file
     --scale - turn on probeset expression scaling
     --exg4 - exclude probes with (G)4 motif
     --dat=N - analyze only a specific dataset number N
     --An=N - calculate (A)n motif statistics for motif of length N (default 24)

     sumarization options:
     --median (default)
     --mean - calculate mean value for each tm instead of median

     sample processing options:
     --single - analyze one sample at a time (default)
     --multi-raw - analyze antire experiment
     --multi-mas - analyze entire experiment processed with MAS5 (MAS-bg corr, MAS-pm corr)
     --multi-rma - analyze entire experiment processed with RMA (RMA-bg corr, quantile norm)
     --multi-gcrma - analyze entire experiment processed with GCRMA (GCRMA-bg corr, quantile norm)
     --multi-sgcrma - analyze entire experiment processed with slow version of GCRMA that models each experiment (GCRMA-bg corr, quantile norm)
     --multi-plier - analyze entire experiment processed with PLIER (quantile norm)
     --multi-farms - analyze entire experiment processed with FARMS (quantile norm)
     --multi-mbei - analyze entire experiment processed with MBEI (invariant-set/dChip/Li-Wong normalization)
    """
def main(argv):
  try:
    #processing options and parameters
    opts, args = getopt.gnu_getopt(argv, "u", ['mean','median','linreg','no-scale','single','multi-raw','multi-rma','multi-gcrma','multi-farms','multi-plier','multi-mas','multi-mbei','exg4','dat='])
  except:
    #incorrect parameters
    print "Error: Incorrect parameters"
    usage()
    sys.exit()

  if (len(args)<4) or (len(args)>4) :
    #application requires 3 parameters
    usage()
    sys.exit()
  else:

    #----starting the script----
    start = time.clock()

    #process all inputs
    expdir=args[0];
    pgc_file=args[1];
    cdfname=args[2];
    outfile=args[3];

    #processing all flags
    procMode=0;   # single-raw
    noscale=True; # dont scale values
    summode=0;    # use median
    exg4=False;   # exclude probes with G4 motif
    Nset=0;       # analyze all sets
    update=False; # dont update the previous results
    global AnLen
    AnLen=24;
    Nmet='RawS'

    for opt,oarg in opts:
       opt=opt.replace('--','-')
       if opt in ('-single'): 
          Nmet='RawS'
          procMode=0;
       if opt in ('-multi-raw'):
          Nmet='RAW'
          procMode=1;
       if opt in ('-multi-rma'):
          Nmet='RMA'
          procMode=2;
       if opt in ('-multi-gcrma'):
          Nmet='GCRMA'
          procMode=3;
       if opt in ('-multi-plier'):
          Nmet='PLIER'
          procMode=4;
       if opt in ('-multi-farms'):
          Nmet='FARMS'
          procMode=5;
       if opt in ('-multi-mas'):
          Nmet='MAS5'
          procMode=6;
       if opt in ('-multi-sgcrma'):
          Nmet='slowGCRMA'
          procMode=7;
       if opt in ('-multi-mbei'):
          Nmet='MBEI'
          procMode=8;
       if opt in ('-scale'):
          noscale=False;
       if opt in ('-mean'):
          summode=1;       
       if opt in ('-exg4'):
          exg4=True;
       if opt in ('-dat'):
          Nset=int(oarg);
       if opt in ('-u'):
          update=True;
       if opt in ('-An'):
          AnLen=int(oarg);

    #assign an analysis summary file
    global report_file
    report_file=expdir+'_'+Nmet+'.report';


    #list all experiments
    os.system('ls %s >explist_%s.tmp' %(expdir,expdir))


    #list experiment that where already analyzed for the update mode
    analyzed_dat=[]
    if update:
       out_file=open(outfile,'r');
       for line in out_file.readlines():
          line=line.strip();
          pline=line.split();
          exper=pline[0];
          if not (exper in analyzed_dat):
             analyzed_dat.append(exper);
    
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
    abort_pgc=False;
    if Nset>0: 
       if Nset-1>=TotalSets: #id above the allowed range
          Nstart=0;
          Nstop=0;
          print 'invalid dataset id: '+str(Nset)
       else:                 #correct id
          Nstart=Nset-1;
          Nstop=Nset;
          
          #abort procedure for PGC read
          exper=exptab[Nstart]; # analyze this experiment       
          if exper in analyzed_dat:
             abort_pgc=True;
             

    #create header if neccesery 
    if not update:
       tmpstr='#Experiment\tSample\tPlatform\tMethod\tTmSlope\tTmIntersec\tMedianExp\tMeanExp\tRawMeanExp\tMeanWithG4\tMeanWithoutG4\tAn_motif_minus\tAn_motif_plus\tAn_motif_minus_sc\tAn_motif_plus_sc';
       for i in range(0,26):
          tmpstr+='\t'+str(64.9+41*(i-16.4)/25);
       if Nset<=1:
          os.system('echo "'+tmpstr+'" > '+outfile);       

    #read the pgc file if neccesery
    if not abort_pgc:
       pgc_gc = {};
       pgc_g4 = {};
       pgc_an = {};
       pgc_id = {};
       Nz=0;
       pgc_file=open(pgc_file,'r');
       for line in pgc_file.readlines():
          line=line.strip();
          pline=line.split('\t');
          if not pline[0][0]=='#':
             prb=pline[0]+':'+pline[1]; 
             pgc_gc[prb]=int(pline[3]);
             pgc_an[prb]=int(pline[5]);          

             #szukamy gggg w sekwencji sondy
             seq=pline[2].lower();
             if seq.find('gggg')+1>0:  
                has_g=1;             
             else:
                has_g=0;
             pgc_g4[prb]=has_g;

             pgc_id[Nz]=prb;         
             Nz+=1;
       pgc_file.close()


    #analyze datasets 
    analyzed_ct=0;
    failed_ct=0;
    for k in range(Nstart,Nstop):
       exper=exptab[k].strip(); # analyze this experiment       
       analyze=True;
       if update:
          if exper in analyzed_dat:          
             analyze=False;
             print '  skipping dataset: %s nr: %s' %(exper,k+1);

       if analyze:          
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
              
              ############### raw data - single sample ######################
              if procMode==0:
               for samp in range(0,Nsamp):
                success=False;
                sampName=exper+' (unknown_sample)'
                try:
                    code='raw<-ReadAffy(filenames=LFS['+str(samp+1)+'],cdfname="'+cdfname+'");'
                    code+='probe <- pm(raw);'
                    code+='sampName=LFS['+str(samp+1)+'];'
                    code+='write.table(probe, file="'+expdir+'/'+exper+'.probedata", sep="\t",  quote=F);'                 
                    result=robjects.r(code);
                    sampName=robjects.r['sampName'];
                    sampName=sampName[0];

                    code='firstVal=probe[1];'
                    code+='avg_rawexp=as.vector(colMeans(pm(raw)));'
                    result=robjects.r(code);
                    firstVal=robjects.r['firstVal'][0];
                    avg_rawexp=robjects.r['avg_rawexp'];                                                      
                    if is_number(firstVal):                                              
                       success=True; 

                except Exception,e:
                    print e
                    print "Failed to analyze: "+sampName
                    update_report('Failed (R)',exper,sampName);
                    failed_ct+=1
               
                if success==True:    
                   #read the data into dictionary
                   fline=True
                   expdata={}
                   data_file=open(expdir+'/'+exper+'.probedata')
                   for line in data_file.readlines():
                      line=line.strip();
                      pline=line.split('\t');
                      if fline:
                         header=pline[0];
                      else:                    
                         expdata[int(pline[0].strip())]=pline[1];
                      fline=False;                 
                   data_file.close()                 
                   os.system('rm '+expdir+'/'+exper+'.probedata')
                   #subprocess.Popen('rm '+expdir+'/'+exper+'.probedata',shell=True);

                   #obliczanie statystyk Tm dla kazdej probki
                   finalRes=calculateTmStats(expdata,1,pgc_id,pgc_gc,pgc_g4,pgc_an,noscale,summode,exg4,avg_rawexp);

                   #generating report
                   tempstr='';
                   for i in range(0,37):
                      tempstr+='\t%.3f' %(finalRes[i,0]);
                   
                   os.system('echo "'+exper+'\t'+header+'\t'+expdir+'\t'+Nmet+tempstr+'" >> ' +outfile);                   
                   #subprocess.Popen('echo "'+exper+'\t'+header+'\t'+expdir+'\t'+Nmet+tempstr+'" >> ' +outfile,shell=True);

                   update_report('Success',exper,sampName);                   
                   analyzed_ct+=1;
                else:
                   print "Failed to analyze: "+sampName+" - damaged CEL file"
                   update_report('Failed (Py)',exper,sampName);                   
                   failed_ct+=1

              ############################# multi - sample #############################
              if procMode>0:
                    #PART I: analysis in R
                    if procMode==3 or procMode==7:
                       code='library(gcrma,verbose=F,warn.conflicts=F);';
                    if procMode==4:
                       code='library(plier,verbose=F,warn.conflicts=F);';
                    if procMode==5:
                       code='library(farms,verbose=F,warn.conflicts=F);';                    
                    code+='library(preprocessCore,verbose=F,warn.conflicts=F);';
                    result=robjects.r(code)
                    
                    success=False;
                    try:
                       if procMode==1:
                          code='raw<-ReadAffy(filenames=LFS, cdfname="'+cdfname+'");'
                          code+='Tab <- pm(raw);'
                          #code='for(i in 1:length(LFS)){'
                          #code+='   raw<-ReadAffy(filenames=LFS[i], cdfname="'+cdfname+'");'
                          #code+='   probe <- pm(raw);'
                          #code+='   if(i==1){Tab=probe;}'
                          #code+='   else {Tab=cbind(Tab,probe);} }'
                          result=robjects.r(code)
                          code='write.table(Tab, file="'+expdir+'/'+exper+'.probedata", sep="\t",  quote=F);'                 
                          result=robjects.r(code)
                          code='firstVal=Tab[1,1];'
                          result=robjects.r(code);
                          
                       else:                       
                          code='raw<-ReadAffy(filenames=LFS,cdfname="'+cdfname+'");' 
                          if procMode==2: #RMA
                             code+='probebg <- apply(pm(raw), 2, bg.adjust);'
                             code+='normdata <- normalize.quantiles(probebg);';
                          if procMode==3: #GC-RMA
                             code+='probebg <- bg.adjust.gcrma(raw);'
                             code+='normdata <- normalize.quantiles(pm(probebg));';
                          if procMode==4: #PLIER 
                             #code+='probebg <- bg.adjust.gcrma(raw);'
                             code+='normdata <- normalize.quantiles(pm(raw));';
                          if procMode==5: #FARMS 
                             code+='normdata <- normalize.quantiles(pm(raw));';
                          if procMode==6: #MAS 
                             code+='probebg <- bg.correct.mas(raw);'
                             code+='normdata <- pmcorrect.mas(probebg)';
                          if procMode==7: #slow GC-RMA
                             code+='probebg <- bg.adjust.gcrma(raw,fast=F);'
                             code+='normdata <- normalize.quantiles(pm(probebg));';
                          if procMode==8: #MBEI
                             code+='normdata <- pm(normalize.AffyBatch.invariantset(raw));';
                          result=robjects.r(code)

                          code='firstVal=normdata[1,1];'
                          result=robjects.r(code);
                          
                          code='rownames(normdata)=rownames(pm(raw));';
                          code+='colnames(normdata)=colnames(pm(raw));';
                          result=robjects.r(code)
                          code='write.table(round(normdata,digits=3), file="'+expdir+'/'+exper+'.probedata", sep="\t",  quote=F);'                 
                          result=robjects.r(code)

                       code='avg_rawexp=as.vector(colMeans(pm(raw)));'
                       result=robjects.r(code);
                       avg_rawexp=robjects.r['avg_rawexp'];                                                                             
                       firstVal=robjects.r['firstVal'][0];                                                      

                       code+='rm(list=ls())'
                       result=robjects.r(code)
                          
                       #check if data format is correct   
                       if is_number(firstVal):
                          success=True;
                       else:
                          print "Error: Data contains NAs"

                    except Exception,e:
                       print e
                       print "Failed to analyze: "+exper
                       update_report('Failed (R)',exper,'');
                       failed_ct+=1;
                    
                    # PART II: data analysis in Py
                    if success==True:  
                       try:                       
                          #read the data into dictionary
                          fline=True
                          expdata={}
                          data_file=open(expdir+'/'+exper+'.probedata','r')
                          for line in data_file.readlines():
                             line=line.strip()
                             pline=line.split('\t')
                             if fline:
                                headers=pline
                             else:                    
                                expdata[int(pline[0].strip())-1]=pline[1:len(pline)]
                             fline=False                 
                          data_file.close()
                          Nsamp=len(headers)
                          os.system('rm '+expdir+'/'+exper+'.probedata')
                          #subprocess.Popen('rm '+expdir+'/'+exper+'.probedata',shell=True);

                          finalRes=calculateTmStats(expdata,Nsamp,pgc_id,pgc_gc,pgc_g4,pgc_an,noscale,summode,exg4,avg_rawexp);
                          expdata={};

                          #generating report
                          for j in range (0,Nsamp):
                             tempstr='';
                             for i in range(0,37):
                                tempstr+='\t%.3f' %(finalRes[i,j]);
                             os.system('echo "'+exper+'\t'+headers[j]+'\t'+expdir+'\t'+Nmet+tempstr+'" >> ' +outfile);      
                             #subprocess.Popen('echo "'+exper+'\t'+headers[j]+'\t'+expdir+'\t'+Nmet+tempstr+'" >> ' +outfile,shell=True);

                          update_report('Success',exper+' ['+str(Nsamp)+']','');
                          analyzed_ct+=1
                    
                       except Exception,e:
                          print e
                          update_report('Failed (Py)',exper+' ['+str(Nsamp)+']','');
                          failed_ct+=1
                    else:
                       print "BioC analysis failed - omitting Python processing"

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

def calculateTmStats(expdata,Nsamp,pgc_id,pgc_gc,pgc_g4,pgc_an,noscale,summode,exg4,avg_rawexp):

   #declering global variable in order to modify it in the function
   max_iter=len(pgc_id);              
              
   ctidx=zeros(26,int);
   tempvec={};
   tempitems=0;

   expRes=zeros([max_iter,Nsamp],float)
   gcRes=zeros(max_iter,int)
   expResCt=0;

   tempexp=zeros(Nsamp,float);
   tmaxexp=zeros(Nsamp,float);
   tminexp=zeros(Nsamp,float);
   fline=True;
   curr_iter=0;
   prev_prbset='';
   used_probes={};

   #parameters for the G4 motif stats
   sum_g4 = zeros(Nsamp,float);
   sum_nog4 = zeros(Nsamp,float);
   sum_an_plus = zeros(Nsamp,float);
   sum_an_minus = zeros(Nsamp,float);
   sum_an_plus_sc = zeros(Nsamp,float);
   sum_an_minus_sc = zeros(Nsamp,float);
   num_g4 = 0;
   num_nog4 = 0;
   num_an_plus = 0;
   num_an_minus = 0;

   for z in range(0,max_iter):
      key=pgc_id[z];
      temp=key.split(':');
      prbset=temp[0];
      prbid=temp[1];
      curr_iter+=1;
      
      if prbset==prev_prbset or fline:                    
         if not(exg4 and pgc_g4[key]>0):                #if the flag exg4 is set and probe containg G4 motif than it will be omitted                          
            tempvec[tempitems]=int(prbid);
            tempitems+=1;
            used_probes[prbid]=1;                       
      if fline:
         prev_prbset=prbset;
         fline=False;

      if prbset!=prev_prbset or curr_iter==max_iter:      
         ## create table with all expr values
         tempTab = zeros([tempitems,Nsamp],float);
         tempexp = zeros(Nsamp+1,float);  #+1 so it will always be a vector not scalar

         ## copy values to temp table and summarize G4 stats
         for i in range(0,tempitems):        #for each probeset probe
            probe = tempvec[i];
            if curr_iter==max_iter:
               probeset = prbset;
            else:
               probeset = prev_prbset;

            #convert and assign values 
            if Nsamp==1:
               tempexp[0]=float(expdata[probe]);
               tempTab[i,0]=tempexp[0];
            else:
               for j in range(0,Nsamp):  #for each sample
                  try:
                     tempexp[j]=float(expdata[probe][j]);
                  except:
                     print "Error: invalid data in sample nr "+str(j+1)+": "+expdata[probe][j];
                     tempexp[j]=100; #substitute
                  tempTab[i,j] = tempexp[j];                                                    

            #summarize G4 probe signals
            ID=probeset+':'+str(probe);
            if pgc_g4[ID]>0:
               for j in range(0,Nsamp):  #for each sample
                  sum_g4[j]+=tempexp[j];
               num_g4+=1;
            else:
               for j in range(0,Nsamp):  #for each sample
                  sum_nog4[j]+=tempexp[j];
               num_nog4+=1;
          
         ## get min and max values of a set
         if tempitems>1:
            for i in range(0,Nsamp):
               tmaxexp[i]=max(tempTab[:,i]);
               tminexp[i]=min(tempTab[:,i]);

         ## summarize (A)n probes
         if tempitems>2:
          for i in range(0,tempitems):    #for each probeset probe
            ID=probeset+':'+str(tempvec[i])
            if pgc_an[ID]>=AnLen:
               #print 'P\t'+ID+'\t'+str((tempTab[i,0]-tminexp[0])/(tmaxexp[0]-tminexp[0]))
               for j in range(0,Nsamp):  #for each sample
                  sum_an_plus_sc[j]+=(tempTab[i,j]-tminexp[j])/(tmaxexp[j]-tminexp[j]);
                  sum_an_plus[j]+=tempTab[i,j];
               num_an_plus+=1;
            elif pgc_an[ID]<=-AnLen:
               #print 'M\t'+ID+'\t'+str((tempTab[i,0]-tminexp[0])/(tmaxexp[0]-tminexp[0]))
               for j in range(0,Nsamp):  #for each sample
                  sum_an_minus_sc[j]+=(tempTab[i,j]-tminexp[j])/(tmaxexp[j]-tminexp[j]);
                  sum_an_minus[j]+=tempTab[i,j];
               num_an_minus+=1;
         
         ## scalling
         if not noscale:
            for j in range(0,Nsamp):  #for each sample
               maxexp=max(tempTab[:,j]);
               for i in range(0,tempitems): #for each probeset probe
                  tempTab[i,j]=100.0*tempTab[i,j]/maxexp;
         
         ## assigning to experiment summary table
         for i in range(0,tempitems): #for each probeset probe
            probe=tempvec[i];
            if not (noscale and (probe in used_probes)):  #probes cannot repeat if scaling is turned off               
               ID=probeset+':'+str(probe);
               gcidx=pgc_gc[ID];
               gcRes[expResCt]=gcidx;
               for j in range(0,Nsamp):  #for each sample
                  expRes[expResCt,j]=tempTab[i,j];            
               ctidx[gcidx]+=1;
               expResCt+=1;

         ## assign value for the next probeset        
         if (exg4 and pgc_g4[key]>0):
            tempitems=0;
         else:
            tempitems=1;
            tempvec[0]=int(prbid);
            used_probes[prbid]=1;

      prev_prbset=prbset;
   
   #calculating median/mean for each Tm
   expRes=expRes[0:expResCt,:]
   gcRes=gcRes[0:expResCt]

   expRes=numpy.array(expRes);
   gcRes=numpy.array(gcRes);
   finalRes=zeros([37,Nsamp],float)
   for i in range(0,26):
      ind=numpy.where(gcRes==i);
      temptab=expRes[ind,:][0];    
      for j in range(0,Nsamp):  #for each sample
         if ctidx[i]>0:
            if summode==0:
               finalRes[i+11,j]=numpy.median(temptab[:,j]);             
            elif summode==1:
               finalRes[i+11,j]=numpy.mean(temptab[:,j]);                 
         else:
            finalRes[i+11,j]=0;
   
   #correlation analysis
   for k in range(0,Nsamp):  #for each sample         
      yvec=expRes[:,k]
      (finalRes[0,k],finalRes[1,k])=polyfit(gcRes,yvec,1);
      finalRes[2,k]=numpy.median(yvec); #median expression
      finalRes[3,k]=numpy.mean(yvec);   #average expression
      finalRes[4,k]=avg_rawexp[k];      #average raw expression 
      finalRes[5,k]=sum_g4[k]/num_g4;     #average expression with G4
      finalRes[6,k]=sum_nog4[k]/num_nog4; #average expression without G4
      finalRes[7,k]=sum_an_minus[k]/num_an_minus;  #average expression before An
      finalRes[8,k]=sum_an_plus[k]/num_an_plus; #average expression after An
      finalRes[9,k]=sum_an_minus_sc[k]/num_an_minus;  #average expression before An (scaled)
      finalRes[10,k]=sum_an_plus_sc[k]/num_an_plus; #average expression after An (scaled)

   return finalRes;


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def update_report(status,exper,sampName):
   #generate data for the report file
   now = datetime.datetime.now();
   text = '%s\tTM\t%s\t%s\t%s ' %(now.strftime("%Y-%m-%d %H:%M"),status,exper,sampName)
   print text
   #subprocess.Popen('echo '+text+'>>'+report_file,shell=True);
   os.system('echo "%s" >> %s' %(text,report_file))


if __name__ == "__main__":
    main(sys.argv[1:])
