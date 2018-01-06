# RARS2.R  Model of RARS   Normal and dysplastic haematopoiesis model      3 January 2018
# MDS progenitor population has (optionally) impaired differentiation
# Tabulated output in tableM.txt
# Set Kgx to 2 for wild type or <2 for mutant ras (line 157). Set TET2 to 1.0 (TET2+) or 0.001 (TET2-) (line 165)
#
setwd("~/githubs/R/RARS")
library(deSolve)
# state variables are capitalized
# P     = progenitors (old CFU-GM)
# Q     = quiescent progenitors
# N     = Neutrophils
# TN    = Tissue Neutrophils
# M     = monocytes
# AN    = Activated Neutrophils
# ERYTH = erythroblasts
# MKP   = megalokaryoblasts
# MKcyte= megalokaryocytes
# Lym   = lymphocytes
# RBC   = circulating red cells
# MDP   = CMML progenitors
# MDQ   = quiescent CMML progenitors
# MDM   = dysplastic monocytes
# MDN   = dysplastic neutrophils
# MDT   = dysplastic tissue neutrophils
# MDA   = activated dysplastic cells
# GM    = GM-CSF   
# STAT5 = STAT5 
# STAT3 = STAT3
# Mcl1  = mantle cell lymphoma 1 (antiapoptosis factor)
# GCSF  = G-CSF
# Pla   = platelets
# Retic = reticulocytes
# epo   = erythropoietin
# Fe2   = Fe2+ in nucleus and cytosol
# MFe   = Fe2+ in mitochondria

# Drugs used to manipulate the system

#mo = MOR103 = Ab blocker of GM   (2)
#das= dasatinib                   (3)
#sb = SB272844 = IL8 inhibitor    (4)
#cn = CNDAC = sapacitabine becomes this = SSB producer = DSBs in S via HR  (5)
#seli = seliciclib = cdk9 inhibitor, blocks transcription of Mcl-1 (6)
#imab = infliximab, anti-TNF moab (7)
#ifna = interferon alpha (8)
#antiox = ascorbic acid, antoxidant (9)
#dac  = decitabine (10)

ix=15   # inner loopindex
nlines=12; dim=(nlines+1)*16
Cmatrix <- matrix(1:dim, ncol=16)

colnames(Cmatrix)<-c("Time","P","Retic","WBC","Platelets","RBC","monocytes","epo","MKcytes","mito-Fe","fe2","N","ERYTH","GM-CSF","MKP","G-CSF")
#
SF3B1=1.0;
SF3B1=0.01;
PI=1.11e7; QI=1.43e6; NI=7.12e6; TNI=1.12e5; ANI=12.8; MDPI=0; MDQI=0; MI=4.2e5; Lym=3.0e6;                                                # starting values
MDNI=0; MDTI=0; MDAI=0; GMI=.05673; STAT5I=.01004; STAT3I=26.59; Mcl1I=2.021; ROSI=1.967; ablI=0; il8I=4; GCSFI=.75; il3I=1; 
ERYTHI=1.27e8; RBCI=5.0e9; MKPI=1.423e4; MKcyteI=1.18e5; PlaI=3.45e8; MDMI=0; ReticI=3.53e8; epoI=.05; Fe2I=5.02; MFeI=11.9     # initial iron concentrations in nucleus and mitochondria
weightP=1.0; weightM=1.0; weightE=0.1; weightMK=1.5             # weighting factors for competition function
CFUlimit=PI*weightP+MDPI*weightM+ERYTHI*weightE+MKPI*weightMK
timex=0; timeW=9999; flag1=1

fmyelo <- function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				vX.P = Vx.p*il3/(il3+Kil3)
				vP.P = P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*(1-Qf)
				vP.N = P*Vp.n*GCSF/(GCSF+Kg)
				vP.M = P*Vp.m/(1+(M+MDM)/KiM)
				vN.TN = N*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)
				vN.X = N*Vn.x/(1+Mcl1/Kmcl1)
				vTN.X = TN*Vtn.x/(1+Mcl1/Kmcl1)
				vTN.AN = TN*Vtn.an*STAT5/(STAT5+Kstat5)
				vAN.TN = AN*Van.tn
				vP.Q = P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*Qf
				vQ.P = Q*Vq.p/(1+P/Kp)
				vX.MDP = Vx.mdp*il3/(il3+Kil3)
				vMDP.MDP = MDP*Vmdp.mdp*GM/(GM+Kgx)/(1+cn/Kcn)*(1-Qf)      #  for wild-type Kgx=2; for mutant ras Kgx<2 
				vMDP.MDN = MDP*Vp.n*GCSF/(GCSF+Kg)*dfactor                 # dfactor = MDS differentiation impairment factor (0-1)
    				vMDP.MDM = MDP*Vp.m/(1+(M+MDM)*(TET2+(1-TET2)*dac/Kidac/(1+dac/Kidac))/KiM)
				vMDN.MDT = MDN*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)		
				vMDN.X   = MDN*Vn.x/(1+Mcl1/Kmcl1)
				vMDT.MDA = MDT*Vtn.an*STAT5/(STAT5+Kstat5)
				vMDT.X = MDT*Vtn.x/(1+Mcl1/Kmcl1)
				vMDA.MDT = MDA*Van.tn
				vMDP.MDQ = MDP*Vmdp.mdp*GM/(GM+Kgx)/(1+cn/Kcn)*Qf
				vMDQ.MDP = MDQ*Vq.p/(1+MDP/Kp)
				vX.GM = Vx.gm/(1+(N+M+MDN+MDM)/Kn)/(1+mo/Kmo)
				vGM.X = GM*Vgm.x
				vX.STAT5=Vx.stat5*(1+abl/(Kabl*(1+das/Kdas)))
				vSTAT5.X=STAT5*Vstat5.x
				vX.STAT3=Vx.stat3*(1+AN/Kan)
				vSTAT3.X=STAT3*Vstat3.x
	                        vX.MCL1=Vx.mcl1*(1+abl/(Kabl2*(1+das/Kdas)))/(1+seli/Kseli)/(1+ifna/Kifna)
				vMCL1.X=Mcl1*Vmcl1.x
				vX.ROS=STAT3*Vx.ros/(STAT3+Kstat3)
				vROS.X=ROS*antiox*Kros2
				vX.ABL=tnfa/(1+imab/Kimab)*Vx.abl/(tnfa+Ktnfa)
				vABL.X=abl*Vabl.x
				vIL8.X=Vil8.x*ROS^nHr/(ROS^nHr+kkmax)
				vERYTH.ERYTH=ERYTH*Veryth.eryth *epo/(epo+Kepo)
				vRBC.X=RBC*Krbc.x
				vMKP.MKP=MKP*Vmkp.mkp/(1+Pla/Kpla)
				vMKP.MKcyte=MKP*Vmkp.mkcyte
				vMKcyte.Pla=MKcyte*Kmkcyte.pla
				vPla.x=Pla*Kpla.x
				vM.X=M*Km.x
				vMDM.X=MDM*Km.x
				vX.GCSF=Vx.gcsf/(1+(N+MDN)/Kineut)
				vGCSF.X=GCSF*Kgcsf.x
				vERYTH.RETIC=ERYTH*Veryth.retic
				vRETIC.X=Retic*Kretic.x
				vRETIC.RBC=Retic*Vretic.rbc*Fe2/5.0
				vX.EPO=Vx.epo/(1+RBC/Krbc)
				vEPO.X=epo*Vepo.x
				vX.Fe2 = Vx.Fe2*Fe2/5.0
				vFe2.X = VFe2.x*vRETIC.RBC
				vFe2.MFe = Fe2*VFe2.mfe 
				vMFe.Fe2 = MFe*Vmfe.Fe2*ABCB7
				
				dP = vX.P + vP.P - vP.N - vP.Q + vQ.P - vP.M
				dQ = vP.Q - vQ.P
				dN = vP.N*volM/volP - vN.TN - vN.X
				dTN = vN.TN*volP/volT - vTN.X - vTN.AN + vAN.TN
				dAN = vTN.AN - vAN.TN
			        dMDP = vX.MDP + vMDP.MDP + vMDQ.MDP - vMDP.MDN - vMDP.MDQ - vMDP.MDM
				dMDQ = vMDP.MDQ - vMDQ.MDP
				dMDN = vMDP.MDN*volM/volP - vMDN.MDT - vMDN.X
				dMDT = vMDN.MDT*volP/volT - vMDT.MDA + vMDA.MDT - vMDT.X
				dMDA = vMDT.MDA*volP/volT - vMDA.MDT
				dGM  = vX.GM-vGM.X
				dSTAT5 =vX.STAT5-vSTAT5.X
				dSTAT3 =vX.STAT3-vSTAT3.X
				dMcl1  =vX.MCL1-vMCL1.X
				dROS   = vX.ROS - vROS.X
				dabl   =vX.ABL-vABL.X
				dIL8   = 0
				dGCSF  = vX.GCSF-vGCSF.X
				dil3   =(vX.ABL-vABL.X)/(Km0*(1+das/Kdas))
				dERYTH = vERYTH.ERYTH-vERYTH.RETIC
				dRBC   = vRETIC.RBC-vRBC.X
				dMKP   = vMKP.MKP-vMKP.MKcyte
				dMKcyte= vMKP.MKcyte-vMKcyte.Pla
				dPla   = vMKcyte.Pla*5000-vPla.x
				dM     = vP.M-vM.X
				dMDM   = vMDP.MDM-vMDM.X
				dRetic = vERYTH.RETIC-vRETIC.X-vRETIC.RBC
				depo   = vX.EPO-vEPO.X
				dFe2   = vX.Fe2-vFe2.X-vFe2.MFe+vMFe.Fe2*volmt/volcs
				dMFe   = vFe2.MFe*volcs/volmt-vMFe.Fe2

				return(list(c(dP,dQ,dN,dTN,dAN,dMDP,dMDQ,dMDN,dMDT,dMDA,dGM,dSTAT5,dSTAT3,dMcl1,dROS,dabl,dIL8,dGCSF,dil3,dERYTH,dRBC,dMKP,dMKcyte,dPla,dM,dMDM,dRetic,depo,dFe2,dMFe),
   						c(vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,vRETIC.RBC=vRETIC.RBC,vERYTH.RETIC=vERYTH.RETIC,vRETIC.X=vRETIC.X))
					   )
			}) 
}

pars=c(Vx.p=1.53e4, Vx.mdp=0,Kil3=1, Vp.p=.236, Vmdp.mdp=.236, Kgm=2, Kgx=2, Kcn=0.1, Vp.n=.026, Kg=1, Kmo=1,Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Kstat5=.1, 
Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=1.96e4, Vgm.x=.05,Vx.stat5=.002, Kabl=1, Kdas=5, Vstat5.x=.2, Vx.stat3=1.5, Kan=100, Vstat3.x=.1, Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, 
Vx.ros=980, Kstat3=250, Kros=30, Kgsh=9.5e5, Vp.m=6e-3,Kros2=.104, Kpla=3e9, Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1,Kidac=1.0, Vabl.x=.0145,
Vil8.x=5e-4, tnfa=.01, bact=0, Vtnfa.x=5e-6, Ka2=.5, Km0=.1, dfactor=1.0, TET2=1, Veryth.eryth=.78, Vretic.rbc=8.4e-2, Krbc=5.0e9, Krbc.x=6.3e-3, Vmkp.mkp=.223, Vmkp.mkcyte=.2, 
Kmkcyte.pla=.029, Kpla.x=.05, Km.x=.012, Kras=1.6, Vx.gcsf=.02, Kgcsf.x=.01, Kineut=5e+6, KiM=3.5e+4, Veryth.retic=0.39, Kretic.x=4.05e-2, Vx.epo=0.02, Vepo.x=0.2, Kepo=.05, Vx.Fe2=7.98e-2,
VFe2.x=2.5e-9, VFe2.mfe=1.52e-2, Vmfe.Fe2=6.4e-2, 
Qf=.0373, volP=2900, volM=1400, volT=65700, volcs=0.9, volmt=0.1, ABCB7=SF3B1,    # these are boundary  conditions
bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7,dac=0) # Drug concns

X0=c(P=PI,Q=QI,N=NI,TN=TNI,AN=ANI,MDP=MDPI,MDQ=MDQI,MDN=MDNI,MDT=MDTI,MDA=MDAI,GM=GMI,STAT5=STAT5I,STAT3=STAT3I,Mcl1=Mcl1I,ROS=ROSI,abl=ablI,il8=il8I,GCSF=GCSFI,il3=il3I,ERYTH=ERYTHI,RBC=RBCI,
MKP=MKPI, MKcyte=MKcyteI,Pla=PlaI,M=MI,MDM=MDMI,Retic=ReticI,epo=epoI,Fe2=Fe2I,MFe=MFeI)
                                                 # Starting values for variables

Cmatrix[1,1]<-0;Cmatrix[1,2]<-PI;Cmatrix[1,3]<-ReticI;Cmatrix[1,4]<-NI+MDNI+MI+MDMI+Lym; Cmatrix[1,5]<-PlaI;Cmatrix[1,6]<-RBCI;Cmatrix[1,7]<-MI+MDMI;Cmatrix[1,8]<-epoI;Cmatrix[1,9]=MKcyteI;
Cmatrix[1,10]<-MFeI;Cmatrix[1,11]<-Fe2I;Cmatrix[1,12]<-NI;Cmatrix[1,13]<-ERYTHI;Cmatrix[1,14]=GMI;Cmatrix[1,15]=MKPI;Cmatrix[1,16]=GCSFI;
									   								   
times <- seq(1,ix, by = 1)
out   <- ode(X0, times, fmyelo, pars)
        
	for(i in seq(from=2, to=nlines, by=1)) {
	P=out[ix,2]; Q=out[ix,3]; N=out[ix,4]; TN=out[ix,5]; AN=out[ix,6]; MDP=out[ix,7]; MDQ=out[ix,8]; MDN=out[ix,9]; MDT=out[ix,10]; MDA=out[ix,11]
	GM=out[ix,12]; STAT5=out[ix,13]; STAT3=out[ix,14]; Mcl1=out[ix,15]; ROS=out[ix,16]; abl=out[ix,17]; il8=out[ix,18]; GCSF=out[ix,19]; il3=out[ix,20]; ERYTH=out[ix,21]; RBC=out[ix,22];
        MKP=out[ix,23];MKcyte=out[ix,24];Pla=out[ix,25];M=out[ix,26];MDM=out[ix,27];Retic=out[ix,28];epo=out[ix,29];Fe2=out[ix,30];MFe=out[ix,31]
	X0=c(P,Q,N,TN,AN,MDP,MDQ,MDN,MDT,MDA,GM,STAT5,STAT3,Mcl1,ROS,abl,il8,GCSF,il3,ERYTH,RBC,MKP,MKcyte,Pla,M,MDM,Retic,epo,Fe2,MFe)
                                                                            
Cmatrix[i,1]<-out[ix,1];Cmatrix[i,2]<-P;Cmatrix[i,3]<-Retic;Cmatrix[i,4]<-N+MDN+M+MDM+Lym; Cmatrix[i,5]<-Pla;Cmatrix[i,6]<-RBC; Cmatrix[i,7]=M+MDM; Cmatrix[i,8]<-out[ix,29];Cmatrix[i,9]=MKcyte; 
Cmatrix[i,10]<-out[ix,31]; Cmatrix[i,11]<-out[ix,30]; Cmatrix[i,12]<-N;Cmatrix[i,13]<-ERYTH; Cmatrix[i,14]=GM; Cmatrix[i,15]=out[ix,23]; Cmatrix[i,16]=out[ix,19];

times <- seq(ix*(i-1)+1,ix*i, by = 1)
out   <- ode(X0, times, fmyelo, pars);P=out[ix,2]; MDP=out[ix,7]; ERYTH=out[ix,21]; MKP=out[ix,23];

# Competition function
progenitors=P*weightP+MDP*weightM+ERYTH*weightE+MKP*weightMK; 
if(progenitors>CFUlimit) {ratioN=P*weightP/progenitors;ratioM=MDP*weightM/progenitors;ratioE=ERYTH*weightE/progenitors;ratioMK=MKP*weightMK/progenitors
        P=ratioN*CFUlimit/weightP; MDP=ratioM*CFUlimit/weightM; ERYTH=ratioE*CFUlimit/weightE; MKP=ratioMK*CFUlimit/weightMK; }

timex=out[ix,1]; if(flag1 && (M+MDM)>4e6){flag1=0; timeW=timex}
if(MDP<.5)MDP=0;
                       }
Cmatrix[i+1,1]=timex;Cmatrix[i+1,2]<-out[ix,2];Cmatrix[i+1,3]<-out[ix,28];Cmatrix[i+1,4]<-out[ix,4]+out[ix,9]+out[ix,26]+out[ix,27]+Lym;Cmatrix[i+1,5]=Pla;Cmatrix[i+1,6]<-out[ix,22]
Cmatrix[i+1,7]<-out[ix,26]+out[ix,27]; Cmatrix[i+1,8]<-out[ix,29];Cmatrix[i+1,9]<-out[ix,24]; Cmatrix[i+1,10]<-out[ix,31];Cmatrix[i+1,11]<-out[ix,30];Cmatrix[i+1,12]<-out[ix,4]
Cmatrix[i+1,13]=out[ix,21];Cmatrix[i+1,14]=out[ix,12];Cmatrix[i+1,15]=out[ix,23]; Cmatrix[i+1,16]=out[ix,19]

write.table(Cmatrix, file="tableM.txt", append=FALSE, quote=FALSE, sep=" ", row.names=TRUE, col.names=TRUE); 
Cmatrix
outline<-c("Time to onset", timeW,     "MDA = ", MDA, "GM = ", GM); #print(outline)

matplot(Cmatrix[,1], Cmatrix[,2:7], type = "l", xlab = "time (days)", ylab = "count",
main = "RARS2.R  3 Jan 2018: Normal haematopoiesis", log="y", lwd = 2)
legend(50,7e8, c("progenitor cells", "reticulocytes","WBC", "platelets","RBC","monocytes"), col = 1:6, lty = 1:6,bty="n")




