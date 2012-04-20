
import ROOT as r
import sys, os

def ensure_dir(path):
    try:
      os.makedirs(path)
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
      else: raise


print " Use Is python subber.py Analysis InPutFileList OutPutDir Trigger" 
# inpu = file("./dec22List.txt","r")
# inpu = file("./incompList.txt","r")
# inpu = file("./HTDataSet.txt","r")
# inpu = file("./Pf.txt","r")
# inpu = file("./newList.txt","r")
# inpu = file("./muList.txt","r")
# inpu = file("./Un.txt","r")
# inpu = file("./Dn.txt","r")
inpu = sys.argv[2]
flist = file(inpu).readlines()
print flist
print "Will submit", len(flist) , "jobs"

for i in range(0,len(flist)):
  os.system("qsub -q hepshort.q subScript.sh " + sys.argv[1] +" "+sys.argv[2] +" " + sys.argv[3] + " " + " "+ sys.argv[4] + " " + str(i))


#Now make an output text file containing the parameters that the job was submitted with
ensure_dir(sys.argv[3])

log = open(sys.argv[3]+"/Log.txt",'w')
text =  sys.argv[1] +" "+sys.argv[2] +" " + sys.argv[3] + " " + " "+ sys.argv[4] + "\n"
log.write(text)

