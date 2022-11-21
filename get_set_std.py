#!/usr/bin/python2

import os, re, math, sys, getopt 
import rpy2.robjects as robjects
import numpy
from numpy import *

def usage():
   print """
   get_set_std.py directory pgc2_file cdf_file 

   DESCRIPTION
     calculates inter-proeset variance for each sample based on specified pgc_file

   INPUT:     
     directory - directory containing individual experiment subdirectories 
     pgc2_file - file containing GC nucleotide counts for each probe - format: ProbeID<tab>Probeset<tab>GC

   OUTPUT:
     individual file for each experiment including STDs for each probe-set (rows) and each sample (columns)
     file naming: experiment_name.setstd
    """
def main(argv):
  breakapp=False
  try:
    #processing options and parameters    
    opts, args = getopt.gnu_getopt(argv,'m')
    #verifying parameters number
    if (len(args)<3) or (len(args)>3):
    	usage()
    	breakapp=True
  except:
    #incorrect parameters
    usage()
    breakapp=True
    print "Error: Incorrect parameters"
  
  if breakapp:
     sys.exit()
    
  
  #----starting the script----
  expdir=args[0];
  pgc_file=args[1];
  cdf_file=args[2]

  nmean=False;
  prefx='';

  for opt,oarg in opts:
    if opt in ("-m"):
       nmean=True;
       prefx='_mean';

  #read the pgc2 file
  pgc_set={};
  pgc_id=zeros(1000000,int);
  Nz=0;
  pgc_file=open(pgc_file,'r');
  for line in pgc_file.readlines():
     line=line.strip();
     pline=line.split('\t');
     if not pline[0][0]=='#':
        prb=int(pline[0].strip())+1;
        pgc_set[Nz]=pline[1];
        pgc_id[Nz]=prb;
        Nz+=1;
  pgc_file.close()

  #list all experiments
  os.system('ls '+expdir+'>explist.tmp')

  #for each experiment from the list
  list_file=open('explist.tmp','r')
  for exper in list_file.readlines():
     exper=exper.strip();
     
     if  os.path.isdir(expdir+'/'+exper):
       if not os.path.exists(expdir+'/'+exper+'.setvar'+prefx):
	    code='library(affy);'                  
            code+='LFS=list.celfiles(path="'+expdir+'/'+exper+'", full.names=TRUE);'
            code+='N=length(LFS);'
            result=robjects.r(code)
	    #code='library(makecdfenv);'
            #cdfName(rawdata)
	    #code+='"'+expdir+'"=make.cdf.env("'+cdf_file+'");'

            result=robjects.r(code)
            Nsamp=robjects.r['N']
            Nsamp=Nsamp[0]
              
            for samp in range(0,Nsamp):
              success=False;
              try:
                code='raw<-ReadAffy(filenames=LFS['+str(samp+1)+']);' #cdfname="hgu133plus2hsrefseq"
                code+='probe <- pm(raw);'
                code+='sampName=LFS['+str(samp+1)+'];'
                code+='write.table(probe, file="'+expdir+'/'+exper+'.probedata", sep="\t",  quote=F);'                 
                result=robjects.r(code);
                success=True ;
                sampName=robjects.r['sampName']
              except:
                sampName=robjects.r['sampName']
                print "Failed to analyze: "+sampName
                
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

                #calculating STD in each probeset              
                max_iter=len(pgc_set);                                            
                tmp_file=open('tmp.setvar','w')
                tmp_file.writelines('%s\t%s\n' %('Probeset',header))
                ctidx=zeros(26,int);
                tempvec=zeros(1000,float);
                tempitems=0;
                fline=True;
                curr_iter=0;
                prev_prbset='';              
                for z in range(0,Nz):
                  key=pgc_id[z];
                  prbset=pgc_set[z];
                  curr_iter+=1;
                  #print z
                  #print prbset

                  if fline:
                     prev_prbset=prbset;
                     fline=False;
                     tempvec[tempitems]=int(key);
                     tempitems+=1;
                  elif prbset==prev_prbset:                    
                     tempvec[tempitems]=int(key);
                     tempitems+=1;
                  elif prbset!=prev_prbset or curr_iter==max_iter:      
                     #create table with all expr values
                     tempTab = zeros([tempitems],float);
                     for i in range(0,tempitems):        #for each probeset probe
                        probe = tempvec[i];
                        tempTab[i] = expdata[probe];                                                    
                        #print '%s\t%.3f\t%s' %(probe,tempTab[i])

                     #calculating STD
                     setstd=numpy.std(tempTab,ddof=1)**2;
                     setmean=numpy.mean(tempTab);
                     
                     if nmean:
                        setstd=setstd/setmean;
                     
                     #saving to temp file
                     tmp_file.writelines('%s\t%.4f\n' %(prev_prbset,setstd))

                     #initializing variables
                     tempitems=1;
                     tempvec[0]=key;

                  prev_prbset=prbset;

                  
                tmp_file.close()

                outfName=expdir+'/'+exper+'.setvar'+prefx;
                if os.path.exists(outfName):
                   os.system('cut -f2 tmp.setvar>tmp')
                   os.system('mv tmp tmp.setvar')
                   os.system('paste '+outfName+' tmp.setvar>tmp')
                   os.system('mv tmp '+outfName)
                else:
                   os.system('mv tmp.setvar '+outfName)
                                      
  list_file.close()


if __name__ == "__main__":
     main(sys.argv[1:])
