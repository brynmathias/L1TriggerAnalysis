#!/usr/bin/env python

from plottingUtils import *

def errorFun(x, par):
  return 0.5*par[0]*(1. + r.TMath.Erf( (x[0] - par[1]) / (math.sqrt(2.)*par[2]) ))



turnOns = {
#"Denomiantor":[List of numerators]  
 "RefJet":["Jet16","Jet36","Jet52",]#"Jet92",],
 "RecoHT":["RecoHTL150","RecoHTL175","RecoHTL1100","RecoHTL1150",],
  "METReference":["MET30Pass","MET50Pass","MET70Pass"],
  "MHTReference":["MHT30Pass","MHT50Pass"],
  "SumEtReference":["SumEt60","SumEt100"],
# Special set so that we dont need to have a different loop for drawing 2d histos as well as other 1d distros
  "NoRatio":["EnCorrelation","L1HT","ResolutionAsFnOfeta","ResolutionAsFnOfpT"],

}


fname = "./L1MuHPF_NoFastJet.root"

def main():
  c1 = Print("L1MuHPF_NoFastJet.pdf")
  c1.DoPageNum = False
  c1.open()
  for key,triggerL in turnOns.iteritems():
    turnOnList = []
    names = []
    xMax = 100.
    xAxisTitle = ""
    if key != "NoRatio": yAxisTitle = "Efficiency"
    if "MET" in key: 
      xAxisTitle = "#slash{E}_{T} (GeV)"
      names = ["L1MET 30","L1MET 50","L1MET 70"]
    if "SumEtReference" in key: 
      xAxisTitle = "#sum E_{T} (GeV)"
      names = ["L1SumET 60","L1SumET 100"]
      xMax = 1000.
    if "RecoHT" in key: 
      xAxisTitle = "H_{T} (GeV)"
      names = ["L1HTT 50","L1HTT 75","L1HTT 100", "L1HTT 150"]
      xMax = 650.
    if "Jet" in key: 
      xAxisTitle = "E_{T} (GeV)"
      names = ["L1 SingleJet 16","L1 SingleJet 36","L1 SingleJet 52","L1 SingleJet 92"]
    if "MHTReference" in key: 
      xAxisTitle = "#slash{H}_{T} (GeV)"
      names = ["L1HTM 30","L1HTM 50","L1HTM 70",]
    
    c1.canvas.SetLogy(False)

    f = r.TFile.Open(fname)
    if "NoRatio" not in key:
      fitText = ""
      denominator = f.Get(key)
      # if "HT"  in key:
      denominator.Rebin(4)
      c = 1
      for trig in triggerL:
        if c == 5: c+=1
        numerator = f.Get(trig)
        fermiFunction = r.TF1("fermiFunction",errorFun,10.,100.,3)
        fermiFunction.SetParameters(1.00,10.,1.)
        fermiFunction.SetParNames("#epsilon","#mu","#sigma")

        
        # if "HT" in key:
        numerator.Rebin(4)
        TurnOn = r.TGraphAsymmErrors()
        TurnOn.Divide(numerator,denominator)
        TurnOn.Fit(fermiFunction,"%f"%(Cor),"%f"%(Cor),10.,100.)
        fitText += "%s #sigma = %f, #mu = %f \n" %(trig,fermiFunction2.GetParameter(2),fermiFunction2.GetParameter(1))
        TurnOn.SetMarkerColor(c)
        TurnOn.SetLineColor(c)
        TurnOn.GetXaxis().SetTitle(xAxisTitle)
        TurnOn.GetYaxis().SetTitle(yAxisTitle)
        c+=1
        TurnOn.Draw("ap")
        turnOnList.append(TurnOn)
        c1.Print()
        c1.canvas.Print("%s.png"%(trig))
      c1.Clear()
      print fitText
      mg = r.TMultiGraph()
      for g in turnOnList:
        mg.Add(g)
      leg = Legend(x1= 0.5, y1=0.25, x2 = 0.8, y2 = 0.8)
      for g,name in zip(turnOnList,names):
        leg.AddEntry(g,name,"pl")
      mg.Draw("ap")
      mg.GetXaxis().SetRangeUser(0.,xMax)
      mg.GetXaxis().SetTitle(xAxisTitle)
      mg.GetYaxis().SetTitle("Efficiency")
      leg.Draw()
      c1.canvas.Print("multiGrap_%s.png"%key)
      c1.Print()
    else:
      for hist in triggerL:
        c1.canvas.SetLogy(False)
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
