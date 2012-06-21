#!/usr/bin/env python

from plottingUtils import *


turnOns = {



#"Denomiantor":[List of numerators]  
 "RefJet":["Jet16","Jet36","Jet52","Jet92",],
 "RecoHT":["RecoHTL150","RecoHTL175","RecoHTL1100","RecoHTL1150",],
# Special set so that we dont need to have a different loop for drawing 2d histos as well as other 1d distros
  "NoRatio":["EnCorrelation","L1HT",],
  "METReference":["MET30Pass","MET50Pass","MET70Pass"]

}


fname = "./rootFiles/test.root"

def main():
  """docstring for main"""
  c1 = Print("out.pdf")
  c1.DoPageNum = False
  c1.open()
  for key,triggerL in turnOns.iteritems():
    turnOnList = []
    c1.canvas.SetLogy(False)
    f = r.TFile.Open(fname)
    if "NoRatio" not in key:
      denominator = f.Get(key)
      # if "HT"  in key:
      denominator.Rebin(10)
      c = 1
      for trig in triggerL:
        if c == 5: c+=1
        numerator = f.Get(trig)
        # if "HT" in key:
        numerator.Rebin(10)
        TurnOn = r.TGraphAsymmErrors()
        # TurnOn.SetTitle(trig)
        TurnOn.Divide(numerator,denominator)
        TurnOn.SetMarkerColor(c)
        TurnOn.SetLineColor(c)
        c+=1
        TurnOn.Draw("ap")
        turnOnList.append(TurnOn)
        c1.Print()
      c1.Clear()
      mg = r.TMultiGraph()
      for g in turnOnList:
        mg.Add(g)
      leg = Legend(x1= 0.5, y1=0.25, x2 = 0.8, y2 = 0.8)
      for g,name in zip(turnOnList,triggerL):
        leg.AddEntry(g,name,"pl")
      mg.Draw("ap")
      mg.GetXaxis().SetRangeUser(0.,300.)
      mg.GetXaxis().SetTitle("E_{T} (GeV)")
      mg.GetYaxis().SetTitle("Efficiency")
      leg.Draw()
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