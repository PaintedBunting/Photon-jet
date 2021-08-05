from ROOT import *
import numpy as np

nbins1x = 3;
nbins2x = 12;
nbins3x = 12;
nbins1y = 96;
nbins2y = 12;
nbins3y = 6;
lvl1 = nbins1x * nbins1y;
lvl2 = nbins2x * nbins2y;
lvl3 = nbins3x * nbins3y;

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


#opening the files. Feel free to remove or add them. 
gfile = TFile("gamma_100GeV_20k.root")
efile = TFile("electron_100GeV_19k.root")
pfile = TFile("pi0_100GeV_20k.root")
a1file = TFile("axion1_100GeV_20k.root")
a2file = TFile("axion2_100GeV_20k.root")

gtree = gfile.Get("fancy_tree")
etree = efile.Get("fancy_tree")
ptree = pfile.Get("fancy_tree")
a1tree = a1file.Get("fancy_tree")
a2tree = a2file.Get("fancy_tree")


#setting up the canvas, pads and corresponding histogram stacks
c = TCanvas("a","a",520,3720)

c.Divide(1, 6, 0.02, 0.01)
one = c.cd(1)
two = c.cd(2)
three = c.cd(3)
four = c.cd(4)
five = c.cd(5)
six = c.cd(6)

Fraction_in_thirdlayer = THStack("","Fraction in the third layer")
Fraction_not_in = THStack("","Fraction not in the third layer")
Middle_lateral_width = THStack("","Middle lateral width")
Front_lateral_width = THStack("","Front lateral width")
Shower_Depth = THStack("","Shower depth")
Shower_Depth_width = THStack("","Shower depth width")


#code for gamma file

gzsegmentation = TH1F("","",3,np.array([-240.,-150.,197.,240.]))
gsampling1_eta = TH2F("","",3,-240.,240.,480//5,-240.,240.)
gsampling2_eta = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
gsampling3_eta = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)

gFraction_in_thirdlayer = TH1F("gamma","",100,0,0.01)
Fraction_in_thirdlayer.Add(gFraction_in_thirdlayer)

gFraction_not_in = TH1F("gamma","",100,0,0.01)
Fraction_not_in.Add(gFraction_not_in)

gMiddle_lateral_width = TH1F("gamma","",100,0,100)
Middle_lateral_width.Add(gMiddle_lateral_width)

gFront_lateral_width = TH1F("gamma","",100,0,100)
Front_lateral_width.Add(gFront_lateral_width)

gShower_Depth = TH1F("gamma","",100,0,1.5)
Shower_Depth.Add(gShower_Depth)

gShower_Depth_width = TH1F("gamma","",100,0,1.0)
Shower_Depth_width.Add(gShower_Depth_width)


for i in range(min(1000,gtree.GetEntries())):
    gtree.GetEntry(i)
    if (i%100==0):
        print(i,gtree.GetEntries())
        pass
    y = "energy"
    exec("%s = %s" % (y,"gtree.cell_0"))
    total_energy = 0.
    third_layer = 0.
    not_in = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    second_layer_X = 0.
    second_layer_X2 = 0.
    first_layer_X = 0.
    first_layer_X2 = 0.
    front_energy = 0.
    middle_energy = 0.
    for j in range(507):
        exec("%s = %s" % (y,"gtree.cell_"+str(j)))
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        gzsegmentation.Fill(gzsegmentation.GetXaxis().GetBinCenter(zbin+1),energy)
        zvalue = gzsegmentation.GetXaxis().GetBinCenter(zbin+1)
        yvalue = 0.;
        xvalue = 0.;
        total_energy+=energy
        lateral_depth+=zbin*energy
        lateral_depth2+=zbin*zbin*energy
        if (zbin==2):
            third_layer+=energy
            pass
        if (xbin < 0 or ybin < 0):
            not_in+=energy
            pass
        if (zbin==0):
            gsampling1_eta.Fill(gsampling1_eta.GetXaxis().GetBinCenter(xbin+1),gsampling1_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = gsampling1_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = gsampling1_eta.GetYaxis().GetBinCenter(ybin+1)
        elif (zbin==1):
            gsampling2_eta.Fill(gsampling2_eta.GetXaxis().GetBinCenter(xbin+1),gsampling2_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = gsampling2_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = gsampling2_eta.GetYaxis().GetBinCenter(ybin+1)
        else:
            gsampling3_eta.Fill(gsampling3_eta.GetXaxis().GetBinCenter(xbin+1),gsampling3_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = gsampling3_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = gsampling3_eta.GetYaxis().GetBinCenter(ybin+1)
            pass
        if (zbin==0):
            first_layer_X += xvalue*energy
            first_layer_X2 += xvalue*xvalue*energy
            front_energy+=energy
        elif (zbin==1):
            second_layer_X += xvalue*energy
            second_layer_X2 += xvalue*xvalue*energy
            middle_energy+=energy
            pass
        pass
    gFraction_in_thirdlayer.Fill(third_layer/total_energy)
    gFraction_not_in.Fill(not_in/total_energy)
    gShower_Depth.Fill(lateral_depth/total_energy)
    if (middle_energy > 0):
        gMiddle_lateral_width.Fill(((second_layer_X2/middle_energy) - (second_layer_X/middle_energy)**2)**0.5)
        pass
    if (front_energy > 0):
        gFront_lateral_width.Fill((first_layer_X2/front_energy - (first_layer_X/front_energy)**2)**0.5)
        pass
    gShower_Depth_width.Fill((lateral_depth2/total_energy - (lateral_depth/total_energy)**2)**0.5)
    pass

gFraction_in_thirdlayer.SetLineColor(kPink)
gFraction_not_in.SetLineColor(kPink)
gMiddle_lateral_width.SetLineColor(kPink)
gFront_lateral_width.SetLineColor(kPink)
gShower_Depth.SetLineColor(kPink)
gShower_Depth_width.SetLineColor(kPink)


#code for electron file

ezsegmentation = TH1F("","",3,np.array([-240.,-150.,197.,240.]))
esampling1_eta = TH2F("","",3,-240.,240.,480//5,-240.,240.)
esampling2_eta = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
esampling3_eta = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)

eFraction_in_thirdlayer = TH1F("electron","",100,0,0.01)
Fraction_in_thirdlayer.Add(eFraction_in_thirdlayer)

eFraction_not_in = TH1F("electron","",100,0,0.01)
Fraction_not_in.Add(eFraction_not_in)

eMiddle_lateral_width = TH1F("electron","",100,0,100)
Middle_lateral_width.Add(eMiddle_lateral_width)

eFront_lateral_width = TH1F("electron","",100,0,100)
Front_lateral_width.Add(eFront_lateral_width)

eShower_Depth = TH1F("electron","",100,0,1.5)
Shower_Depth.Add(eShower_Depth)

eShower_Depth_width = TH1F("electron","",100,0,1.0)
Shower_Depth_width.Add(eShower_Depth_width)


for i in range(min(1000,etree.GetEntries())):
    etree.GetEntry(i)
    if (i%100==0):
        print(i,etree.GetEntries())
        pass
    y = "energy"
    exec("%s = %s" % (y,"etree.cell_0"))
    total_energy = 0.
    third_layer = 0.
    not_in = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    second_layer_X = 0.
    second_layer_X2 = 0.
    first_layer_X = 0.
    first_layer_X2 = 0.
    front_energy = 0.
    middle_energy = 0.
    for j in range(507):
        exec("%s = %s" % (y,"etree.cell_"+str(j)))
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        ezsegmentation.Fill(ezsegmentation.GetXaxis().GetBinCenter(zbin+1),energy)
        zvalue = ezsegmentation.GetXaxis().GetBinCenter(zbin+1)
        yvalue = 0.;
        xvalue = 0.;
        total_energy+=energy
        lateral_depth+=zbin*energy
        lateral_depth2+=zbin*zbin*energy
        if (zbin==2):
            third_layer+=energy
            pass
        if (xbin < 0 or ybin < 0):
            not_in+=energy
            pass
        if (zbin==0):
            esampling1_eta.Fill(esampling1_eta.GetXaxis().GetBinCenter(xbin+1),esampling1_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = esampling1_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = esampling1_eta.GetYaxis().GetBinCenter(ybin+1)
        elif (zbin==1):
            esampling2_eta.Fill(esampling2_eta.GetXaxis().GetBinCenter(xbin+1),esampling2_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = esampling2_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = esampling2_eta.GetYaxis().GetBinCenter(ybin+1)
        else:
            esampling3_eta.Fill(esampling3_eta.GetXaxis().GetBinCenter(xbin+1),esampling3_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = esampling3_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = esampling3_eta.GetYaxis().GetBinCenter(ybin+1)
            pass
        if (zbin==0):
            first_layer_X += xvalue*energy
            first_layer_X2 += xvalue*xvalue*energy
            front_energy+=energy
        elif (zbin==1):
            second_layer_X += xvalue*energy
            second_layer_X2 += xvalue*xvalue*energy
            middle_energy+=energy
            pass
        pass
    eFraction_in_thirdlayer.Fill(third_layer/total_energy)
    eFraction_not_in.Fill(not_in/total_energy)
    eShower_Depth.Fill(lateral_depth/total_energy)
    if (middle_energy > 0):
        eMiddle_lateral_width.Fill(((second_layer_X2/middle_energy) - (second_layer_X/middle_energy)**2)**0.5)
        pass
    if (front_energy > 0):
        eFront_lateral_width.Fill((first_layer_X2/front_energy - (first_layer_X/front_energy)**2)**0.5)
        pass
    eShower_Depth_width.Fill((lateral_depth2/total_energy - (lateral_depth/total_energy)**2)**0.5)
    pass

eFraction_in_thirdlayer.SetLineColor(kGreen)
eFraction_not_in.SetLineColor(kGreen)
eMiddle_lateral_width.SetLineColor(kGreen)
eFront_lateral_width.SetLineColor(kGreen)
eShower_Depth.SetLineColor(kGreen)
eShower_Depth_width.SetLineColor(kGreen)


#code for pi0 file

pzsegmentation = TH1F("","",3,np.array([-240.,-150.,197.,240.]))
psampling1_eta = TH2F("","",3,-240.,240.,480//5,-240.,240.)
psampling2_eta = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
psampling3_eta = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)

pFraction_in_thirdlayer = TH1F("pi0","",100,0,0.01)
Fraction_in_thirdlayer.Add(pFraction_in_thirdlayer)

pFraction_not_in = TH1F("pi0","",100,0,0.01)
Fraction_not_in.Add(pFraction_not_in)

pMiddle_lateral_width = TH1F("pi0","",100,0,100)
Middle_lateral_width.Add(pMiddle_lateral_width)

pFront_lateral_width = TH1F("pi0","",100,0,100)
Front_lateral_width.Add(pFront_lateral_width)

pShower_Depth = TH1F("pi0","",100,0,1.5)
Shower_Depth.Add(pShower_Depth)

pShower_Depth_width = TH1F("pi0","",100,0,1.0)
Shower_Depth_width.Add(pShower_Depth_width)


for i in range(min(1000,ptree.GetEntries())):
    ptree.GetEntry(i)
    if (i%100==0):
        print(i,ptree.GetEntries())
        pass
    y = "energy"
    exec("%s = %s" % (y,"ptree.cell_0"))
    total_energy = 0.
    third_layer = 0.
    not_in = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    second_layer_X = 0.
    second_layer_X2 = 0.
    first_layer_X = 0.
    first_layer_X2 = 0.
    front_energy = 0.
    middle_energy = 0.
    for j in range(507):
        exec("%s = %s" % (y,"ptree.cell_"+str(j)))
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        pzsegmentation.Fill(pzsegmentation.GetXaxis().GetBinCenter(zbin+1),energy)
        zvalue = pzsegmentation.GetXaxis().GetBinCenter(zbin+1)
        yvalue = 0.;
        xvalue = 0.;
        total_energy+=energy
        lateral_depth+=zbin*energy
        lateral_depth2+=zbin*zbin*energy
        if (zbin==2):
            third_layer+=energy
            pass
        if (xbin < 0 or ybin < 0):
            not_in+=energy
            pass
        if (zbin==0):
            psampling1_eta.Fill(psampling1_eta.GetXaxis().GetBinCenter(xbin+1),psampling1_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = psampling1_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = psampling1_eta.GetYaxis().GetBinCenter(ybin+1)
        elif (zbin==1):
            psampling2_eta.Fill(psampling2_eta.GetXaxis().GetBinCenter(xbin+1),psampling2_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = psampling2_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = psampling2_eta.GetYaxis().GetBinCenter(ybin+1)
        else:
            psampling3_eta.Fill(psampling3_eta.GetXaxis().GetBinCenter(xbin+1),psampling3_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = psampling3_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = psampling3_eta.GetYaxis().GetBinCenter(ybin+1)
            pass
        if (zbin==0):
            first_layer_X += xvalue*energy
            first_layer_X2 += xvalue*xvalue*energy
            front_energy+=energy
        elif (zbin==1):
            second_layer_X += xvalue*energy
            second_layer_X2 += xvalue*xvalue*energy
            middle_energy+=energy
            pass
        pass
    pFraction_in_thirdlayer.Fill(third_layer/total_energy)
    pFraction_not_in.Fill(not_in/total_energy)
    pShower_Depth.Fill(lateral_depth/total_energy)
    if (middle_energy > 0):
        pMiddle_lateral_width.Fill(((second_layer_X2/middle_energy) - (second_layer_X/middle_energy)**2)**0.5)
        pass
    if (front_energy > 0):
        pFront_lateral_width.Fill((first_layer_X2/front_energy - (first_layer_X/front_energy)**2)**0.5)
        pass
    pShower_Depth_width.Fill((lateral_depth2/total_energy - (lateral_depth/total_energy)**2)**0.5)
    pass

pFraction_in_thirdlayer.SetLineColor(kCyan)
pFraction_not_in.SetLineColor(kCyan)
pMiddle_lateral_width.SetLineColor(kCyan)
pFront_lateral_width.SetLineColor(kCyan)
pShower_Depth.SetLineColor(kCyan)
pShower_Depth_width.SetLineColor(kCyan)


#code for axion1 file

a1zsegmentation = TH1F("","",3,np.array([-240.,-150.,197.,240.]))
a1sampling1_eta = TH2F("","",3,-240.,240.,480//5,-240.,240.)
a1sampling2_eta = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
a1sampling3_eta = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)

a1Fraction_in_thirdlayer = TH1F("axion1","",100,0,0.01)
Fraction_in_thirdlayer.Add(a1Fraction_in_thirdlayer)

a1Fraction_not_in = TH1F("axion1","",100,0,0.01)
Fraction_not_in.Add(a1Fraction_not_in)

a1Middle_lateral_width = TH1F("axion1","",100,0,100)
Middle_lateral_width.Add(a1Middle_lateral_width)

a1Front_lateral_width = TH1F("axion1","",100,0,100)
Front_lateral_width.Add(a1Front_lateral_width)

a1Shower_Depth = TH1F("axion1","",100,0,1.5)
Shower_Depth.Add(a1Shower_Depth)

a1Shower_Depth_width = TH1F("axion1","",100,0,1.0)
Shower_Depth_width.Add(a1Shower_Depth_width)


for i in range(min(1000,a1tree.GetEntries())):
    a1tree.GetEntry(i)
    if (i%100==0):
        print(i,a1tree.GetEntries())
        pass
    y = "energy"
    exec("%s = %s" % (y,"a1tree.cell_0"))
    total_energy = 0.
    third_layer = 0.
    not_in = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    second_layer_X = 0.
    second_layer_X2 = 0.
    first_layer_X = 0.
    first_layer_X2 = 0.
    front_energy = 0.
    middle_energy = 0.
    for j in range(507):
        exec("%s = %s" % (y,"a1tree.cell_"+str(j)))
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        a1zsegmentation.Fill(a1zsegmentation.GetXaxis().GetBinCenter(zbin+1),energy)
        zvalue = a1zsegmentation.GetXaxis().GetBinCenter(zbin+1)
        yvalue = 0.;
        xvalue = 0.;
        total_energy+=energy
        lateral_depth+=zbin*energy
        lateral_depth2+=zbin*zbin*energy
        if (zbin==2):
            third_layer+=energy
            pass
        if (xbin < 0 or ybin < 0):
            not_in+=energy
            pass
        if (zbin==0):
            a1sampling1_eta.Fill(a1sampling1_eta.GetXaxis().GetBinCenter(xbin+1),a1sampling1_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a1sampling1_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a1sampling1_eta.GetYaxis().GetBinCenter(ybin+1)
        elif (zbin==1):
            a1sampling2_eta.Fill(a1sampling2_eta.GetXaxis().GetBinCenter(xbin+1),a1sampling2_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a1sampling2_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a1sampling2_eta.GetYaxis().GetBinCenter(ybin+1)
        else:
            a1sampling3_eta.Fill(psampling3_eta.GetXaxis().GetBinCenter(xbin+1),a1sampling3_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a1sampling3_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a1sampling3_eta.GetYaxis().GetBinCenter(ybin+1)
            pass
        if (zbin==0):
            first_layer_X += xvalue*energy
            first_layer_X2 += xvalue*xvalue*energy
            front_energy+=energy
        elif (zbin==1):
            second_layer_X += xvalue*energy
            second_layer_X2 += xvalue*xvalue*energy
            middle_energy+=energy
            pass
        pass
    a1Fraction_in_thirdlayer.Fill(third_layer/total_energy)
    a1Fraction_not_in.Fill(not_in/total_energy)
    a1Shower_Depth.Fill(lateral_depth/total_energy)
    if (middle_energy > 0):
        a1Middle_lateral_width.Fill(((second_layer_X2/middle_energy) - (second_layer_X/middle_energy)**2)**0.5)
        pass
    if (front_energy > 0):
        a1Front_lateral_width.Fill((first_layer_X2/front_energy - (first_layer_X/front_energy)**2)**0.5)
        pass
    a1Shower_Depth_width.Fill((lateral_depth2/total_energy - (lateral_depth/total_energy)**2)**0.5)
    pass

a1Fraction_in_thirdlayer.SetLineColor(kOrange)
a1Fraction_not_in.SetLineColor(kOrange)
a1Middle_lateral_width.SetLineColor(kOrange)
a1Front_lateral_width.SetLineColor(kOrange)
a1Shower_Depth.SetLineColor(kOrange)
a1Shower_Depth_width.SetLineColor(kOrange)


#code for axion2 file

a2zsegmentation = TH1F("","",3,np.array([-240.,-150.,197.,240.]))
a2sampling1_eta = TH2F("","",3,-240.,240.,480//5,-240.,240.)
a2sampling2_eta = TH2F("","",480//40,-240.,240.,480//40,-240.,240.)
a2sampling3_eta = TH2F("","",480//40,-240.,240.,480//80,-240.,240.)

a2Fraction_in_thirdlayer = TH1F("axion2","",100,0,0.01)
Fraction_in_thirdlayer.Add(a2Fraction_in_thirdlayer)

a2Fraction_not_in = TH1F("axion2","",100,0,0.01)
Fraction_not_in.Add(a2Fraction_not_in)

a2Middle_lateral_width = TH1F("axion2","",100,0,100)
Middle_lateral_width.Add(a2Middle_lateral_width)

a2Front_lateral_width = TH1F("axion2","",100,0,100)
Front_lateral_width.Add(a2Front_lateral_width)

a2Shower_Depth = TH1F("axion2","",100,0,1.5)
Shower_Depth.Add(a2Shower_Depth)

a2Shower_Depth_width = TH1F("axion2","",100,0,1.0)
Shower_Depth_width.Add(a2Shower_Depth_width)


for i in range(min(1000,a2tree.GetEntries())):
    a2tree.GetEntry(i)
    if (i%100==0):
        print(i,a2tree.GetEntries())
        pass
    y = "energy"
    exec("%s = %s" % (y,"a2tree.cell_0"))
    total_energy = 0.
    third_layer = 0.
    not_in = 0.
    lateral_depth = 0.
    lateral_depth2 = 0.
    second_layer_X = 0.
    second_layer_X2 = 0.
    first_layer_X = 0.
    first_layer_X2 = 0.
    front_energy = 0.
    middle_energy = 0.
    for j in range(507):
        exec("%s = %s" % (y,"a2tree.cell_"+str(j)))
        xbin = get_x(j,get_y(j,get_z(j)),get_z(j))
        ybin = get_y(j,get_z(j))
        zbin = get_z(j)
        a2zsegmentation.Fill(a2zsegmentation.GetXaxis().GetBinCenter(zbin+1),energy)
        zvalue = a2zsegmentation.GetXaxis().GetBinCenter(zbin+1)
        yvalue = 0.;
        xvalue = 0.;
        total_energy+=energy
        lateral_depth+=zbin*energy
        lateral_depth2+=zbin*zbin*energy
        if (zbin==2):
            third_layer+=energy
            pass
        if (xbin < 0 or ybin < 0):
            not_in+=energy
            pass
        if (zbin==0):
            a2sampling1_eta.Fill(a2sampling1_eta.GetXaxis().GetBinCenter(xbin+1),a2sampling1_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a2sampling1_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a2sampling1_eta.GetYaxis().GetBinCenter(ybin+1)
        elif (zbin==1):
            a2sampling2_eta.Fill(a2sampling2_eta.GetXaxis().GetBinCenter(xbin+1),a2sampling2_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a2sampling2_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a2sampling2_eta.GetYaxis().GetBinCenter(ybin+1)
        else:
            a2sampling3_eta.Fill(psampling3_eta.GetXaxis().GetBinCenter(xbin+1),a2sampling3_eta.GetYaxis().GetBinCenter(ybin+1),energy)
            xvalue = a2sampling3_eta.GetXaxis().GetBinCenter(xbin+1)
            yvalue = a2sampling3_eta.GetYaxis().GetBinCenter(ybin+1)
            pass
        if (zbin==0):
            first_layer_X += xvalue*energy
            first_layer_X2 += xvalue*xvalue*energy
            front_energy+=energy
        elif (zbin==1):
            second_layer_X += xvalue*energy
            second_layer_X2 += xvalue*xvalue*energy
            middle_energy+=energy
            pass
        pass
    a2Fraction_in_thirdlayer.Fill(third_layer/total_energy)
    a2Fraction_not_in.Fill(not_in/total_energy)
    a2Shower_Depth.Fill(lateral_depth/total_energy)
    if (middle_energy > 0):
        a2Middle_lateral_width.Fill(((second_layer_X2/middle_energy) - (second_layer_X/middle_energy)**2)**0.5)
        pass
    if (front_energy > 0):
        a2Front_lateral_width.Fill((first_layer_X2/front_energy - (first_layer_X/front_energy)**2)**0.5)
        pass
    a2Shower_Depth_width.Fill((lateral_depth2/total_energy - (lateral_depth/total_energy)**2)**0.5)
    pass

a2Fraction_in_thirdlayer.SetLineColor(kViolet)
a2Fraction_not_in.SetLineColor(kViolet)
a2Middle_lateral_width.SetLineColor(kViolet)
a2Front_lateral_width.SetLineColor(kViolet)
a2Shower_Depth.SetLineColor(kViolet)
a2Shower_Depth_width.SetLineColor(kViolet)



#setting up and drawing pads and stacks
one.SetLogz(0)
two.SetLogz(0)
three.SetLogz(0)
four.SetLogz(0)
five.SetLogz(0)
six.SetLogz(0)

c.cd(1)
Fraction_in_thirdlayer.Draw("nostack") 
one.BuildLegend(0.75,0.7,0.9,0.9)

c.cd(2)
Fraction_not_in.Draw("nostack") 
two.BuildLegend(0.75,0.7,0.9,0.9)

c.cd(3)
Shower_Depth.Draw("nostack")
three.BuildLegend(0.1,0.7,0.25,0.9)

c.cd(4)
Middle_lateral_width.Draw("nostack")
four.BuildLegend(0.75,0.7,0.9,0.9)

c.cd(5)
Front_lateral_width.Draw("nostack") 
five.BuildLegend(0.75,0.7,0.9,0.9)

c.cd(6)
Shower_Depth_width.Draw("nostack")  
six.BuildLegend(0.75,0.7,0.9,0.9)

c.Update()
c.Draw()
