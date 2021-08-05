#converts file to format to be used by bdt

from ROOT import *
import numpy as np
from array import array

nbins1x = 3
nbins2x = 12
nbins3x = 12
nbins1y = 96
nbins2y = 12
nbins3y = 6
lvl1 = nbins1x * nbins1y
lvl2 = nbins2x * nbins2y
lvl3 = nbins3x * nbins3y

def get_z(myindex):
    if myindex == 504:
        return 0
    elif myindex == 505:
        return 1
    elif myindex == 506:
        return 2
    elif (myindex >= lvl1 + lvl2):
        return 2
    elif (myindex >= lvl1):
        return 1
    else:
        return 0
    pass

def get_y(myindex,zbin):
    if (myindex >=504):
        return -1
    elif (zbin==0):
        return myindex % nbins1y
    elif (zbin==1):
        return myindex % nbins2y
    else:
        return myindex % nbins3y
    pass

def get_x(myindex,ybin,zbin):
    if (myindex >=504):
        return -1
    elif (zbin==0):
        return (myindex - ybin)//nbins1y
    elif (zbin==1):
        return (myindex - lvl1 - ybin)//nbins2y
    else:
        return (myindex - lvl1 - lvl2 - ybin)//nbins3y
    pass


#def sendhelp(inname = "axion1_100GeV_20k.root", output_prefix = "new"):
inname = "gamma_100GeV_20k.root"
output_prefix = "new"
treename = "fancy_tree"
outname = output_prefix+"_"+inname
fin = TFile(inname)
tin = fin.Get(treename)
fout = TFile(outname,"recreate");
tout= TTree(treename,treename)
    
#assignments so python doesn't complain when setting up the branches.
sampling1 = TH2F("","",3,-240.,240.,480//5,-240.,240.)
sampling2 = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
sampling3 = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)
    
#y = "energy"
#exec("%s = %s" % (y,"tin.cell_0"))
#print("initial energy is ",energy)

#cells = np.empty(507)

total_e = array('f', [ 0. ])
firstlayer_e = array('f', [ 0. ])
secondlayer_e = array('f', [ 0. ])
thirdlayer_e = array('f', [ 0. ])
lateral_depth = array('f', [ 0. ])
lateral_depth2 = array('f', [ 0. ])
firstlayer_x = array('f', [ 0. ])
firstlayer_x2 = array('f', [ 0. ])
secondlayer_x = array('f', [ 0. ])
secondlayer_x2 = array('f', [ 0. ])
front_energy = array('f', [ 0. ])
middle_energy = array('f', [ 0. ])
frac_first = array('f', [ 0. ])
frac_second = array('f', [ 0. ])
frac_third = array('f', [ 0. ])
shower_depth_gamma = array('f', [ 0. ])
second_lateral_width_gamma = array('f', [ 0. ])
first_lateral_width_gamma = array('f', [ 0. ])
shower_depth_width_gamma = array('f', [ 0. ])

#branch definitions
tout.Branch("total_e", total_e,"total_e/F");
tout.Branch("firstlayer_e",firstlayer_e,"firstlayer_e/F");
tout.Branch("secondlayer_e",secondlayer_e,"secondlayer_e/F");
tout.Branch("thirdlayer_e",thirdlayer_e,"thirdlayer_e/F");
tout.Branch("lateral_depth",lateral_depth,"lateral_depth/F");
tout.Branch("lateral_depth2",lateral_depth2,"lateral_depth2/F");
tout.Branch("firstlayer_x",firstlayer_x,"firstlayer_x/F");
tout.Branch("firstlayer_x2",firstlayer_x2,"firstlayer_x2/F");
tout.Branch("secondlayer_x",secondlayer_x,"secondlayer_x/F");
tout.Branch("secondlayer_x2",secondlayer_x2,"secondlayer_x2/F");
tout.Branch("frac_first",frac_first,"frac_first/F");
tout.Branch("frac_second",frac_second,"frac_second/F");
tout.Branch("frac_third",frac_third,"frac_third/F");
tout.Branch("shower_depth_gamma",shower_depth_gamma,"shower_depth_gamma/F");
tout.Branch("second_lateral_width_gamma",second_lateral_width_gamma,"second_lateral_width_gamma/F");
tout.Branch("first_lateral_width_gamma",first_lateral_width_gamma,"first_lateral_width_gamma/F");
tout.Branch("shower_depth_width_gamma",shower_depth_width_gamma,"shower_depth_width_gamma/F");

for i in range(tin.GetEntries()):
    tin.GetEntry(i)
        
    y = "energy"
    exec("%s = %s" % (y,"tin.cell_0"))
    #print("cell initial energy is ",energy)
    
    total_e = 0.
    firstlayer_e = 0.
    secondlayer_e = 0.
    thirdlayer_e = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    firstlayer_x = 0.
    firstlayer_x2 = 0.
    secondlayer_x = 0.
    secondlayer_x2 = 0.

    if i%100 == 0:
        print(i,tin.GetEntries())
        
    for j in range(504):
        exec("%s = %s" % (y,"tin.cell_"+str(j)))
        #print("filled energy is ",energy)
        
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        
        lateral_depth+= zbin*energy
        lateral_depth2+= zbin*zbin*energy
        yvalue = 0.;
        xvalue = 0.;
        
        total_e += energy
        
        if (zbin==0):
            firstlayer_e += energy
            sampling1.Fill(sampling1.GetXaxis().GetBinCenter(xbin+1),sampling1.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = sampling1.GetXaxis().GetBinCenter(xbin+1)
            yvalue = sampling1.GetYaxis().GetBinCenter(ybin+1)
            firstlayer_x += xvalue * energy
            firstlayer_x2 += xvalue * xvalue * energy
            
        elif (zbin==1):
            secondlayer_e += energy
            sampling2.Fill(sampling2.GetXaxis().GetBinCenter(xbin+1),sampling2.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = sampling2.GetXaxis().GetBinCenter(xbin+1)
            yvalue = sampling2.GetYaxis().GetBinCenter(ybin+1)
            secondlayer_x += xvalue * energy
            secondlayer_x2 += xvalue * xvalue * energy
            
        elif (zbin == 2):
            thirdlayer_e += energy
            sampling3.Fill(sampling3.GetXaxis().GetBinCenter(xbin+1),sampling3.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = sampling3.GetXaxis().GetBinCenter(xbin+1)
            yvalue = sampling3.GetYaxis().GetBinCenter(ybin+1)
        pass

    frac_first=firstlayer_e/total_e
    frac_second=secondlayer_e/total_e
    frac_third=thirdlayer_e/total_e
    shower_depth_gamma=lateral_depth/total_e
    second_lateral_width_gamma = ((secondlayer_x2/secondlayer_e) - (secondlayer_x/secondlayer_e)**2)**0.5
    first_lateral_width_gamma = (firstlayer_x2/firstlayer_e - (firstlayer_x/firstlayer_e)**2)**0.5
    shower_depth_width_gamma = (lateral_depth2/total_e - (lateral_depth/total_e)**2)**0.5
    tout.Fill()
    pass

tout.Write()
fout.Close()
