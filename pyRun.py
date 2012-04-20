import ROOT as r
import sys


inpu = sys.argv[2]
Analysis = sys.argv[1]
outFile = sys.argv[3]
outDir = sys.argv[3]
Trigger = sys.argv[4]
r.gROOT.ProcessLine(".x ~/cms04/CMSSW_5_2_3_patch3/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C")
flist = file(inpu).readlines()
print Analysis, inpu, outFile, Trigger
idx = sys.argv[5]
#print int(idx) , "INDEX!!!!"
f = flist[int(idx)]
fi = f.strip()
print "file" , [f][0:-1] , "outdir" , str(outDir)+"/"+str(fi).rpartition("/")[2] 

print Analysis , inpu, outFile,Trigger, idx
r.gROOT.ProcessLine(".L "+str(Analysis+"+"))
if Analysis == "L1EnergySumAnalysis": m = r.L1EnergySumAnalysis()
if Analysis == "L1JetAnalysis.C" : m = r.L1JetAnalysis()
m.Open(f[0:-1])
m.run(-1,str(outDir)+"/"+str(fi).rpartition("/")[2], Trigger)
#m.run(910000,str(outFile),bool(useUnCor))
#m.Delete()
