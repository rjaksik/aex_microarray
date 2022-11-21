#!/usr/bin/python2

#Convert files to xda
#ls |grep "E-" |sed 's/E-/apt-cel-convert -i -f xda E-/g' |sed ':a;N;$!ba;s/\n/\/*\n/g' >TMP

#clean data out of files other than CEL
#grep DatHeader -a Rat230_2/*  -r > Rat230_2.LISTA
#grep "Rat230_2.1sq" Rat230_2.LISTA -v | grep LISTA -v |cut -f1 -d":" | sed 's/Rat230_2/rm Rat230_2/g'  >Rat230_2.RM


import os, re, math, sys, getopt

def usage():
   print """
   get_aex_data.py platform_ID platform_name output_directory

   DESCRIPTION
     downloads the entire dataset from the ArrayExpress database for a specific platform    

   INPUT:
     platform_ID - as used by the ArrayExpress, ex.: A-AFFY-34
     platform_name - as used by affymetrix, ex.: HG-U133A
     output_directory

   FLAGS:
     --update=year - update the dataset with experiments conducted on specific year only

   OUTPUT:
     folders divided by experiment ID, ex.: E-GEOD-28072 containing only the CEL files for specified platform
    """

def main(argv):
  try:
    #processing options and parameters
    opts, args = getopt.gnu_getopt(argv, "", ['update='])
  except:
    #incorrect parameters
    print "Error: Incorrect parameters"
    usage()
    sys.exit()

  if (len(args)<3) or (len(args)>3) :
    #application requires 3 parameters
    usage()
    sys.exit()
  else:

    #----starting the script----

    platform_ID=args[0]
    platform_name=args[1]
    outdir=args[2]

    update=0;
    for opt,oarg in opts:
       opt=opt.replace('--','-')
       if opt in ('-update'):
          update=oarg;


    if not os.path.exists(outdir):
       os.system('mkdir '+outdir)

    log_file=outdir+'/aex_log.'+platform_name
    os.system('echo>'+log_file)

    #download experiments list
    os.system('wget "http://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=&species=&array=%s&sortby=releasedate&sortorder=descending&expandefo=on" -O %s.aex_list' %(platform_ID,outdir+'/'+platform_ID))
    
    list_file = outdir+'/'+platform_ID+'.aex_list';

    #passing experiments conducted only on specific day
    if update<>0:
       print 'Downloading experiments conducted in year: '+update;
       os.system('grep '+update+'- '+list_file+'>'+list_file+'.tmp')
       os.system('mv '+list_file+'.tmp '+list_file)
       

    list_file = open(list_file,'r')

    #for each experiment from the list
    for line in list_file.readlines():
           pline=line.split('\t')
           experiment=pline[0].strip()
           if experiment!='Accession' and experiment!='':
              #check if raw data is available
              if pline[6]=='Data is not available':
                 os.system('echo %s\t0\t0\t"%s">>%s' %(experiment,pline[1],log_file))
              #do not download experiments exceeding 1000 samples
              elif int(pline[4])>1000: 
                 os.system('echo %s\t%s\t0\t"%s">>%s' %(experiment,pline[2],pline[1],log_file))
              #check if we already have this experiment
              elif not os.path.exists(outdir+'/'+experiment):
                 os.system('mkdir %s/%s' %(outdir,experiment))
                
                 #for each package file in the experiment
                 for i in range(1,25):
                    #download package file (*.zip) and save it to
                    zipfile='%s.raw.%d.zip' %(experiment,i);
                    os.system('wget http://www.ebi.ac.uk/arrayexpress/files/%s/%s -P %s/%s' %(experiment,zipfile,outdir,experiment))
                    #extract each package file and remove the zip file
                    if os.path.exists(outdir+'/'+experiment+'/'+zipfile):
                       os.system('unzip %s/%s/%s -d %s/%s' %(outdir,experiment,zipfile,outdir,experiment))
                       os.system('rm %s/%s/%s' %(outdir,experiment,zipfile))
                       os.system('apt-cel-convert -f xda -i %s/%s/*.CEL' %(outdir,experiment))
                       os.system('apt-cel-convert -f xda -i %s/%s/*.cel' %(outdir,experiment))
  
                 #find and remove CEL files for other platform   
                 os.system('grep DatHeader %s/%s/* -a | cut -f3- -d"/" |grep "%s.1sq" -v -i>%s.cellist' %(outdir,experiment,platform_name,experiment))
                 os.system('cut -f1 %s.cellist -d":" >%s.delete' %(experiment,experiment))
                 cpn_file = open('%s.delete' %(experiment))
                 count_lines = len(cpn_file.readlines());
                 cpn_file.close()
                 cpn_file = open('%s.delete' %(experiment))
                 for line2 in cpn_file.readlines():
                    line2=line2.strip();
                    os.system('rm %s/%s/%s' %(outdir,experiment,line2))                    
                 cpn_file.close()
                 os.system('echo %s\t%s\t%d\t"%s">>%s' %(experiment,pline[4],int(pline[4])-count_lines,pline[1],log_file))
                 os.system('rm %s.delete' %(experiment))
                 os.system('rm %s.cellist' %(experiment))

                 #find and remove files other than CEL
                 os.system('find %s/%s ! -iname "*.cel"  -type f -exec rm {} \;' %(outdir,experiment));

                    
    list_file.close()


if __name__ == "__main__":
    main(sys.argv[1:])


