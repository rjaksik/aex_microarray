#!/usr/bin/python2

import os, re, math



def usage():
   print """
   get_cel_dates.py directory output_file

   DESCRIPTION
     ectracts sample processing dates (yy/mm/dd) out of each cell file in specified directory

   INPUT:
     directory - directory containing the cell files in possible subdirectories

   OUTPUT:
     file_name date
    """

def main(*args,**opts):
 if (len(args)<2) or (len(args)>2):
    usage()
 else:
    useddir=args[0]
    out_file_name=args[1]

    fname=out_file_name.split('/');
    fname=fname[len(fname)-1];

    print fname

    #get all headers
    if not os.path.exists('lista.tmp.'+fname):
       os.system('grep DatHeader %s/* -a -r -i>lista.tmp.%s' % (useddir,fname))

    #analize each header
    list_file,out_file = open('lista.tmp.'+fname,'r'),open(out_file_name,'w')
    for line in list_file.readlines():
           pline=line.split(':')
           file_name=pline[0].strip()
           if len(pline)>2:
              #converting MM/DD/YY to YY/MM/DD
              temp_str=pline[2].strip();
              idx=temp_str.find('/');
              date=temp_str[idx+4:idx+6]+'/'+temp_str[idx-2:idx+3];
              dateCtrl=date.count('/');
              if dateCtrl<>2:
                 date='NA'
              if len(pline)>=4:
                 time=temp_str[idx+7:idx+9]+':'+pline[3]+':'+pline[4][0:2];
              else:
                 time='NA';
              out_file.writelines(file_name+'\t'+date+'\t'+time+'\n');
           else:
              out_file.writelines(file_name+'\tNA\n')
    list_file.close()
    out_file.close()    


def getOptsArgs2(opts='',long_opts=[]):
    import getopt,re,sys
    inp = sys.argv[1:]

    ## pre-process input
    for i in range(len(inp)):
        # convert - to {}
        if inp[i]=='-':
            inp[i]='{}'
        # convert -longarg to --longarg
        if inp[i].lstrip('-') in [x.rstrip('=') for x in long_opts] and \
           len(inp[i]) >= 2 and inp[i][0]=='-' and inp[i][1]!='-':
            # check for ambiguity
            import warnings 
            warnings.simplefilter('ignore', DeprecationWarning)
            import sets
            warnings.resetwarnings()
            longchars = sets.Set(inp[i].lstrip('-'))
            optchars =sets.Set(opts.replace(':',''))
            if longchars.issubset(optchars):
                raise Exception('Ambiguous parameter %s: Use -%s or change the order of your flags\n' % (inp[i],inp[i]))
            else:
                inp[i]='-%s' % inp[i]
    ## parse input
    o, a = getopt.gnu_getopt(inp,opts,long_opts)

    o = dict([(x[0].strip('-'),x[1]) for x in o])
    return o, a


if __name__=='__main__':
    import sys
    opts=''      # example: 'at:v'
    long_opts=['up','down','gene','cdsup','cdsdown'] # example: ['test=','blah','boo=']
    o, a = getOptsArgs2(opts,long_opts)
    main(*a,**dict(o))
