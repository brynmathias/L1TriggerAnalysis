#!/usr/bin/env python

from plottingUtils import *
r.gStyle.SetOptFit(0)
def errorFun(x, par):
  return 0.5*par[0]*(1. + r.TMath.Erf( (x[0] - par[1]) / (math.sqrt(2.)*par[2]) ))



# "EmulatedJet16","EmulatedJet20","EmulatedJet36","EmulatedJet52","EmulatedJet68","EmulatedJet92",

turnOns = {
#"Denomiantor":[List of numerators]  
 "RefJet":["Jet16","Jet36","Jet52","Jet92"],#,"EmulatedJet20","EmulatedJet36","EmulatedJet52","EmulatedJet68","EmulatedJet92",],
 # "RecoHT":["RecoHTL150","RecoHTL175","RecoHTL1100","RecoHTL1150",],
 # "RefJetEmu":["EmulatedJet16","EmulatedJet36","EmulatedJet52"],
  # "METReference":["MET30Pass","MET50Pass","MET70Pass"],
  # "MHTReference":["MHT30Pass","MHT50Pass"],
  # "SumEtReference":["SumEt60","SumEt100"],
# Special set so that we dont need to have a different loop for drawing 2d histos as well as other 1d distros
  # "NoRatio":["L1CorVsUnCor","RecoJetvsL1Jet","RecoJetvsL1JetEmu"],

}


files = ["./IsoMu24_L1Fast",]
#fname = "./SingleMu.root"
def main():
  for fname in files:
    
    c1 = Print("%s.pdf"%fname)
    c1.DoPageNum = False
    fitText = ""
    everything = []
    for key,triggerL in turnOns.iteritems():
    
      turnOnList = []
      names = []
      xMax = 150.
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
      if "DeltaR" in key:
        xAxisTitle = "#delta R"
      c1.canvas.SetLogy(False)

      f = r.TFile.Open(fname+".root")
      if "NoRatio" not in key:
        denominator = f.Get(key)
        # if "HT"  in key:
        denominator.Rebin(4)
        everything.append(denominator)
        c = 1
        for trig in triggerL:
          fitLow = 0.
          fitHigh = 500.
          fitMid = 10.
          if "Jet16" in trig: 
              fitLow = 15.
              fitMid = 15.
              fitHigh = 70.
          if "Jet36" in trig: 
              fitLow = 15.
              fitMid = 60.
              fitHigh = 100.
          if  "Jet52" in trig: 
              fitLow = 20.
              fitHigh = 140.
          if  "Jet92" in trig: 
              fitLow = 20.
              fitHigh = 200.
              fitMid = 100.
          if trig == "RecoHTL150": 
              fitMid = 75.
              fitHigh = 400.
          if trig == "RecoHTL175": 
              fitMid = 100.
              fitHigh = 400.              
          if trig == "RecoHTL1100": 
              fitMid = 175.
              fitHigh = 400.
          if trig == "RecoHTL1150": 
              fitMid = 250.
              fitHigh = 400.
          if trig == "MET30Pass":
            fitHigh = 200.
            fitMid = 70.
          if trig == "MET50Pass":
            fitHigh = 200.
            fitMid = 80.
          if trig == "MET70Pass":
            fitHigh = 200.
            fitMid = 100.
          if c == 5: c+=1
          
          numerator = f.Get(trig)
          everything.append(numerator)
          fermiFunction = r.TF1("fermiFunction",errorFun,fitLow,fitHigh,3)
          fermiFunction.SetParameters(1.00,fitMid,1.)
          fermiFunction.SetParNames("#epsilon","#mu","#sigma")
          everything.append(fermiFunction)
        
          # if "HT" in key:
          numerator.Rebin(4)
          TurnOn = r.TGraphAsymmErrors()
          # TurnOn.SetTitle(trig)
          TurnOn.Divide(numerator,denominator)
          TurnOn.Fit(fermiFunction,"0.","1.",fitLow,fitHigh)
          fitText += "%s #sigma = %f, #mu = %f #epslion = %f, Chi2 = %f\n" %(trig,fermiFunction.GetParameter(2),fermiFunction.GetParameter(1),fermiFunction.GetParameter(0),fermiFunction.GetChisquare()/fermiFunction.GetNDF())
          TurnOn.SetMarkerColor(c)
          TurnOn.SetLineColor(c)
          TurnOn.GetXaxis().SetTitle(xAxisTitle)
          TurnOn.GetYaxis().SetTitle(yAxisTitle)
          TurnOn.GetYaxis().SetTitleSize(TurnOn.GetXaxis().GetTitleSize())
          TurnOn.GetYaxis().SetTitleOffset(0.8)
          c+=1
          TurnOn.Draw("ap")
          TurnOn.GetXaxis().SetRangeUser(0.,200.)
          
          turnOnList.append(TurnOn)
          fermiFunction.Draw("same")
          c1.Print()
          # c1.canvas.Print("%s.png"%(trig))
        c1.Clear()
        mg = r.TMultiGraph()
        for g in turnOnList:
          mg.Add(g)
        leg = Legend(x1= 0.5, y1=0.25, x2 = 0.8, y2 = 0.8)
        for g,name in zip(turnOnList,names):
          leg.AddEntry(g,name,"pl")
        mg.Draw("ap")
        mg.GetXaxis().SetRangeUser(0.,220.)
        mg.GetXaxis().SetTitle(xAxisTitle)
        mg.GetYaxis().SetTitleSize(mg.GetXaxis().GetTitleSize())
        mg.GetYaxis().SetTitleOffset(0.8)
        mg.GetYaxis().SetTitle("Efficiency")
        leg.Draw()
        # c1.canvas.Print("multiGrap_%s.png"%key)
        c1.Print()
      else:
        for hist in triggerL:
          xAxisTitle=None
          yAxisTitle=None
          
          if "RecoJetvsL1Jet" in hist:
              xAxisTitle = "Offline Jet E_{T} (GeV)"
              yAxisTitle = "L1 Jet E_{T} (GeV)"
              
          c1.canvas.SetLogy(False)
          hist = f.Get(hist)
          hist.Scale(100./hist.Integral())
          # hist.SetTitle("")
          if xAxisTitle is None:xAxisTitle=hist.GetXaxis().GetTitle()
          if yAxisTitle is None:yAxisTitle=hist.GetYaxis().GetTitle()
          if isinstance(hist,r.TH2F):
            hist.GetXaxis().SetRangeUser(0.,254.)
            hist.GetYaxis().SetRangeUser(0.,254.)
            hist.GetXaxis().SetTitle(xAxisTitle)
            hist.GetYaxis().SetTitle(yAxisTitle)
            hist.GetZaxis().SetTitle("A.U.")
            hist.Draw("COLZ")
            c1.Print()

          if isinstance(hist,r.TH1F):
            c1.canvas.SetLogy()
            hist.Draw("hist")
            hist.GetXaxis().SetRangeUser(10.,300.)
            hist.GetXaxis().SetTitle(xAxisTitle)
            hist.GetYaxis().SetTitle("Events / 0.053")
            c1.Print()

    
    print fitText
    c1.Print()    
    c1.close()
  
  
  
  
if __name__ == "__main__":
  main()
