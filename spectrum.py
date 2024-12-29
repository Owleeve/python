#!/usr/bin
import numpy as np
import math
for i in range(0,4):
    inputfile_energy="./SiO2_Bi%s/excitation%s.out" %(i,i)  #The variables in parentheses correspond to the preceding variables
    inputfile_Dipole="./SiO2_Bi%s/Dipole"%i
    outputfile_discrete="./SiO2_Bi%s/spectrum_discrete_Bi%s.dat"%(i,i)
    outputfile_continuous="./SiO2_Bi%s/spectrum_continuous_Bi%s.dat"%(i,i)
    row=1
    start_row=1
    end_row=300
    len_energy=300     #Here the length of excitation_energy and Dipole
    #extract the excitation energy
    inputf_energy=open(inputfile_energy,"r")
    energy=[]
    for inputline in inputf_energy.readlines():
        if row >= start_row and row <= end_row:  #prevent the out_range of list
            inputline=inputline.split()
            energy.append(inputline[5])  # 5 is the energy position in excitation_Bi%s.out
            row=row+1
        else:
            row=row+1
    inputf_energy.close()
    energy=list(map(float,energy))  # change str into float
    for ii in range(len_energy):
        energy[ii]=energy[ii]*13.6056923   #list can't multiply the float correspondly,so we need use for cycle

    #extract the Diple^2
    inputf_Dipole=open(inputfile_Dipole,"r")
    Dipole=[]
    row=1                           
    for inputline in inputf_Dipole.readlines():
        if row >= start_row and row <= end_row:
            inputline=inputline.split()
            Dipole.append(inputline[3])  # 3 is the tatal Mu^2 in Dipole file
            row=row+1
        else:
            row=row+1
    inputf_Dipole.close()
    Dipole=list(map(float,Dipole))

    #get the excitation wavelength
    wavelength=[]
    for E in energy:
        lamb=1240.7/E               #change energy(eV) into wavelengtg(nm)
        wavelength.append(lamb)
    
    #get the spectrum_discrete
    outputf_discrete=open(outputfile_discrete,"a")
    outputf_discrete.write("excitation_energy(eV)  excitation_wavelength(nm) absorption_coefficient \n")
    for j in range(len_energy):
        absorption=energy[j]*Dipole[j]
        outputf_discrete.write(str(energy[j])+" "+str(wavelength[j])+" "+str(absorption)+"\n") #only str is applied to the write(), + can join the str
    outputf_discrete.close()
 
    #get the spectrum_continuous
    Emax=energy[len_energy-1]
    Emin=energy[0]
    wgrid=0.001
    sigma=0.01
    nw=math.floor((Emax-Emin)/wgrid)  #floor(x) is the integer <= x
    Asorption=np.zeros((nw,3))      #construct the 2-dimension array
    for k in range(len_energy):
        for iw in range(nw):
            w=iw*wgrid+Emin-energy[k]
            Asorption[iw][0]=iw*wgrid+Emin
            Asorption[iw][1]=1240.7/(iw*wgrid+Emin)
            Asorption[iw][2]+=sigma/(pow(w,2)+pow(sigma,2))*Dipole[k]*(iw*wgrid+Emin)
    outputf_continuous=open(outputfile_continuous,"a")
    outputf_continuous.write("energy(eV)  wavelength(nm) absorption_coefficient \n")
    for iw in range(nw):
        outputf_continuous.write(str(Asorption[iw][0])+" "+str(Asorption[iw][1])+" "+str(Asorption[iw][2])+"\n")
    outputf_continuous.close()
       

