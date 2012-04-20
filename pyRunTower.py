import ROOT as r
import sys


r.gROOT.ProcessLine(".x ~/cms03/CMSSW_4_2_5/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C")

r.gROOT.ProcessLine(".L L1TowerAnalysis.C+")
m = r.L1TowerAnalysis()
n = r.L1TowerAnalysis()
o = r.L1TowerAnalysis()
m.OpenWithList("Tower0.txt")
n.OpenWithList("Tower5.txt")
o.OpenWithList("Tower10.txt")
m.run(-1,"./Tower0.root","HLT_ZeroBias")
n.run(-1,"./Tower5.root","HLT_ZeroBias")
o.run(-1,"./Tower10.root","HLT_ZeroBias")