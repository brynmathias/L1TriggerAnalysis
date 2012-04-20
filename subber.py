#!/usr/bin/env python
import ROOT as r
import os, commands, sys
import argparse
import errno
def ensure_dir(path):
    try:
      os.makedirs(path)
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
      else: raise




def main():
  """docstring for main"""
  parser = argparse.ArgumentParser(description="Run the L1 analysis on the batch queue.")
  parser.add_argument("--input",type = str,dest='inputFile')
  parser.add_argument("--output",type=str,default="./",dest='outputDir')
  parser.add_argument("--trigger",type=str,default="HLT_IsoMu*",dest='Trigger')
  parser.add_argument("--queue",type=str,default = "hepshort.q",dest='queue')
  parser.add_argument("--analysis",type=str,default="L1JetAnalysis",dest='analysis')
  args = parser.parse_args()

  flist = file(args.inputFile).readlines()
  print flist
  print "Will submit", len(flist) , "jobs"
  ensure_dir("./%s/"%args.outputDir)
  text = "%s subScript.sh %s %s %s %s"%(args.queue,args.analysis,args.inputFile,args.outputDir,args.Trigger)
  log = open("./%s/Log.txt"%args.outputDir,'w')
  log.write(text)
  for i in range(0,len(flist)):
    os.system("qsub -q %s subScript.sh %s %s %s %s %s"%(args.queue,args.analysis,args.inputFile,args.outputDir,args.Trigger,i))



if __name__ == '__main__':
  main()

