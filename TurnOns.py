#!/usr/bin/env python

from plottingUtils import *


turnOns = {
 "RefJet":["JetEt10","JetEt20","JetEt30","JetEt40","JetEt50","JetEt60","JetEt70","JetEt80","JetEt90","JetEt100","JetEt110","JetEt120","JetEt130","JetEt140","JetEt150"],
  "NoRatio":["EnCorrelation","L1HT",],
  "RecoHT":["RecoHTL1100","RecoHTL1150","RecoHTL175","RecoHTL150","RecoHTL1175","RecoHTL1125",]
}


fname = "SingleMu_L1Trigger.root"

def main():
  """docstring for main"""
  c1 = Print("out.pdf")
  c1.open()
  for key,triggerL in turnOns.iteritems():

    c1.canvas.SetLogy(False)
    f = r.TFile.Open(fname)
    if "NoRatio" not in key:
      denominator = f.Get(key)
      # if "HT"  in key:
      denominator.Rebin(10)
      for trig in triggerL:

        numerator = f.Get(trig)
        # if "HT" in key:
        numerator.Rebin(10)
        TurnOn = r.TGraphAsymmErrors()
        TurnOn.SetTitle(trig)
        TurnOn.Divide(numerator,denominator)
        TurnOn.Draw("ap")
        c1.Print()
    else:
      for hist in triggerL:
        hist = f.Get(hist)
        if isinstance(hist,r.TH2F):
          hist.GetXaxis().SetRangeUser(0.,300.)
          hist.GetYaxis().SetRangeUser(0.,300.)
          hist.Draw("COLZ")
          c1.Print()

        if isinstance(hist,r.TH1F):
          c1.canvas.SetLogy()
          hist.Draw("hist")
          hist.GetXaxis().SetRangeUser(10.,300.)
          c1.Print()
          
  c1.close()
  pass
  
  
  
if __name__ == "__main__":
  main()