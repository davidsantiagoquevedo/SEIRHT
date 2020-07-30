#!/usr/bin/python
# -*- coding:utf8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
BBDD_Dir='BBDD/'
"""
En este script se generaliza la implementación del modelo con PCR fijo
distribuido entre IC e IM priorizando IC hasta Agosto 1 y Utilizando
RT para IC desde esta fecha
update:
27/07/2020: se agregá compartimento ICA
"""

########################################################################
#################            SEIRHT CLASS              #################
########################################################################
class SEIRHT():
	def __init__(self):
		self.dummy=True
	
	def var_trans(self,beta_0,beta_1,beta_H,r,A,Movilidad):
		self.beta_0=beta_0
		self.beta_1=beta_1
		self.beta_H=beta_H
		self.r=r
		self.A=A
		self.Movilidad=Movilidad
		
	def beta(self,t):
		if t<92:
			if self.A==0:
				At=self.A
				return (1-At)*self.beta_0 + At*self.beta_1		
			else:
				At=self.Movilidad[int(t)]
				return (1-At)*self.beta_0 + At*self.beta_1	
		if t>=92:
			At=self.A
			return (1-At)*self.beta_0 + At*self.beta_1		
	
	def var_t_estadia(self, omega, gamma_M, sigma_C, sigma_CA, gamma_HR, 
		nu, gamma_R, sigma_HD, sigma_UD):
		self.omega=omega		#T. prom. latencia
		self.gamma_M=gamma_M	#T. prom de recuperación para IM
		self.sigma_C=sigma_C	#T. antes de aislamiento IC
		self.sigma_CA=sigma_CA	#T. en aislamiento ICA
		self.gamma_HR=gamma_HR	#T. prom HR->R
		self.nu=nu				#T. prom HU->IR
		self.gamma_R=gamma_R	#T. prom IR->R
		self.sigma_HD=sigma_HD	#T. prom HD->D
		self.sigma_UD=sigma_UD	#T. prom UD->D
	
	def var_H_UCI(self, delta_M, delta_HR, delta_HD, delta_UR):
		self.delta_M=delta_M
		self.delta_HR=delta_HR
		self.delta_HD=delta_HD
		self.delta_UR=delta_UR
		self.delta_UD=1-delta_HR-delta_HD-delta_UR
		
	def var_testeo(self, alpha, q, positividad, N_PCR, xi_PCR, xi_RT):
		self.alpha=alpha
		self.q=q
		self.positividad=positividad
		self.N_PCR=N_PCR
		self.xi_PCR=xi_PCR
		self.xi_RT=xi_RT
			
	def theta(self, t, IM, ICA):
		#Sin PCR ni RT hasta Junio 1 (día 92)
		if t<92:
			return 0
		
		#3000 PCR desde Junio 1 (día 92) hasta Julio 1 (día 122), se utilizan para IM e IC
		if 92<=t<122:
			if (self.N_PCR[0]/self.positividad[0] - self.sigma_CA*ICA)<0:
				return 0
			else:
				if (self.alpha*IM)<(self.N_PCR[0]/self.positividad[0] - self.sigma_CA*ICA):
					return 1
				else:
					return (self.N_PCR[0]/self.positividad[0] - self.sigma_CA*ICA)/(self.alpha*IM)
		
		#5000 PCR desde Julio 1 (día 122) hasta Agosto 1 (día 153), se utilizan para IM e IC
		if 122<=t<153:
			if (self.N_PCR[1]/self.positividad[1] - self.sigma_CA*ICA)<0:
				return 0
			else:
				if (self.alpha*IM)<(self.N_PCR[1]/self.positividad[1] - self.sigma_CA*ICA):
					return 1
				else:
					return (self.N_PCR[1]/self.positividad[1] - self.sigma_CA*ICA)/(self.alpha*IM)
			
		#15000 PCR desde Agosto 1 (día 153) para IM y RT para IM, hasta Septiembre 1 (día 183)
		if 153<=t<183:
			if (self.N_PCR[2]/self.positividad[2] - (1-self.xi_RT)*self.sigma_CA*ICA)<0:
				return 0
			else:
				if (self.alpha*IM)<(self.N_PCR[2]/self.positividad[2] - (1-self.xi_RT)*self.sigma_CA*ICA):
					return 1
				else:
					return (self.N_PCR[2]/self.positividad[2] - (1-self.xi_RT)*self.sigma_CA*ICA)/(self.alpha*IM)
					
		#10000 tests desde Septiembre 1 (día 183) 
		if 183<=t:
			if (self.N_PCR[3]/self.positividad[3] - (1-self.xi_RT)*self.sigma_CA*ICA)<0:
				return 0
			else:
				if (self.alpha*IM)<(self.N_PCR[3]/self.positividad[3] - (1-self.xi_RT)*self.sigma_CA*ICA):
					return 1
				else:
					return (self.N_PCR[3]/self.positividad[3] - (1-self.xi_RT)*self.sigma_CA*ICA)/(self.alpha*IM)
					
	def var_ini(self, N0, E0=0, IM0=0, IMT0=0, IMNT0=0, IC0=0, ICA0=0, 
				IHR0=0, IUR0=0, IHD0=0, IUD0=0, IR0=0, R0=0, D0=0 ):
		
		self.N0=N0
		self.E0=E0
		
		self.IM0=IM0
		self.IMT0=IMT0
		self.IMNT0=IMNT0
		
		self.IC0=IC0
		self.ICA0=ICA0
				
		self.IHR0=IHR0
		self.IUR0=IUR0
		self.IHD0=IHD0
		self.IUD0=IUD0
		self.IR0=IR0
		self.R0=R0
		self.D0=D0
		self.S0= (self.N0 - self.E0 - self.IM0 - self.IMT0  - self.IMNT0 
				  - self.IC0  - self.ICA0 - self.IHR0 - self.IUR0 - self.IHD0 
				  - self.IUD0 - self.IR0 - self.R0 - self.D0)
	
	def ODES(self,y,t):
		S, E, IM, IMT, IMNT, IC, ICA, IHR, IUR, IHD, IUD, IR, R, D, N = y
		
		dSdt = -S/float(N)*(self.beta(t)*(IM+IC+(1-self.r)*ICA+IMNT+(1-self.q)*IMT)+self.beta_H*(IHR+IUR+IHD+IUD+IR))
		dEdt =  S/float(N)*(self.beta(t)*(IM+IC+(1-self.r)*ICA+IMNT+(1-self.q)*IMT)+self.beta_H*(IHR+IUR+IHD+IUD+IR)) - self.omega*E
		
		dIMdt = self.delta_M*self.omega*E - self.alpha*IM
		dIMTdt = self.theta(t, IM, ICA)*self.alpha*IM - self.gamma_M*IMT
		dIMNTdt = (1-self.theta(t, IM, ICA))*self.alpha*IM - self.gamma_M*IMNT
		
		dICdt = (1-self.delta_M)*self.omega*E - self.sigma_C*IC
		dICAdt = self.sigma_C*IC - self.sigma_CA*ICA
				
		dIHRdt = self.delta_HR*self.sigma_CA*ICA - self.gamma_HR*IHR
		dIURdt = self.delta_UR*self.sigma_CA*ICA - self.nu*IUR
		dIHDdt = self.delta_HD*self.sigma_CA*ICA - self.sigma_HD*IHD
		dIUDdt = self.delta_UD*self.sigma_CA*ICA - self.sigma_UD*IUD
		
		dIRdt = self.nu*IUR - self.gamma_R*IR
		dRdt = self.gamma_HR*IHR + self.gamma_R*IR + self.gamma_M*IMT + self.gamma_M*IMNT
		
		dDdt = self.sigma_HD*IHD + self.sigma_UD*IUD
		dNdt = -self.sigma_HD*IHD - self.sigma_UD*IUD
		
		return [dSdt, dEdt, dIMdt, dIMTdt, dIMNTdt, dICdt, dICAdt,
				dIHRdt, dIURdt, dIHDdt, dIUDdt, dIRdt, dRdt, dDdt, dNdt]
		
	def solve(self,t0,tf):
		self.t0=t0
		self.tf=tf
		y0= [self.S0, self.E0, self.IM0, self.IMT0, self.IMNT0, 
			 self.IC0, self.ICA0, self.IHR0, self.IUR0, 
			 self.IHD0, self.IUD0, self.IR0, self.R0, self.D0, self.N0]
		t_vect= np.linspace(self.t0, self.tf, self.tf*100)
		self.t_=t_vect
		solution= odeint(self.ODES,y0,t_vect)
		
		self.S_vect=solution.T[0]
		self.E_vect=solution.T[1]
		
		self.IM_vect=solution.T[2]
		self.IMT_vect=solution.T[3]
		self.IMNT_vect=solution.T[4]
		
		self.IC_vect=solution.T[5]
		self.ICA_vect=solution.T[6]
				
		self.IHR_vect=solution.T[7]
		self.IUR_vect=solution.T[8]
		self.IHD_vect=solution.T[9]
		self.IUD_vect=solution.T[10]
		
		self.IR_vect=solution.T[11]
		self.R_vect=solution.T[12]
		self.D_vect=solution.T[13]
		self.N_vect=solution.T[14]
		
	def Test_contar(self):
		self.PCR_IC=[]
		self.PCR_IM=[]
		self.PCR_RTO=[]
		self.RT=[]
		for i in range(len(self.t_)):
			if self.t_[i]<122: pos=self.positividad[0] #junio
			if 122<=self.t_[i]<153: pos=self.positividad[1] #julio
			if 153<=self.t_[i]<183: pos=self.positividad[2] #agosto
			if 183<=self.t_[i]: pos=self.positividad[3] #septiembre
			if self.t_[i]<153:
				self.PCR_IC.append(self.sigma_CA*self.ICA_vect[i]*pos)
				self.PCR_IM.append(self.alpha*self.theta(self.t_[i],self.IM_vect[i],self.ICA_vect[i])*self.IM_vect[i]*pos)
				self.RT.append(0)
			if self.t_[i]>=153:
				self.PCR_IC.append((1-self.xi_RT)*self.sigma_CA*self.ICA_vect[i]*pos)
				self.PCR_IM.append(self.alpha*self.theta(self.t_[i],self.IM_vect[i],self.ICA_vect[i])*self.IM_vect[i]*pos)
				self.RT.append(self.sigma_CA*self.ICA_vect[i]*pos)
		
		
########################################################################
#################               PARAMETERS             #################
########################################################################
Indice_Movilidad= pd.read_excel(str(BBDD_Dir)+'/A.xlsx')
A_vect=Indice_Movilidad['Movilidad'][:]

model=SEIRHT()

#1. Transmisión
model.var_trans(beta_0=1.3743,		#R0=3.0
				#beta_0=1.1439,		#R0=2.5
				beta_1=0.019899665, 
				beta_H=0.01,		#Contacto  HR, UR, HD, UD, R
				r=0.65,
				A=0.6,
				Movilidad= A_vect
				)
				
#2. Tiempos de estadía
model.var_t_estadia(omega=1/4.6,		#T. prom. latencia
					gamma_M=1/1.1,		#T. prom IM
					sigma_C=1/3.,		#T. días en IC antes de aislamiento
					sigma_CA=1/4.1,		#T. días en aislamiento ICA
					gamma_HR=1/9.5,		#T. prom HR->R
					nu=1/11.3,			#T. prom HU->IR
					gamma_R=1/3.4,		#T. prom IR->R
					sigma_HD=1/7.6,		#T. prom HD->D
					sigma_UD=1/10.1,	#T. prom UD->D
					)

#3. Prob Transición Hospital y UCI
model.var_H_UCI(delta_M = 0.965578477,
				delta_HR = 0.696594546,
				delta_UR = 0.122566499,
				delta_HD = 0.058272457
				)

#4. Prob. testeo
#N_PCR=[3962,7991,15000,15000]
N_PCR=[2384,2824,15000,15000]
model.var_testeo(alpha=1/1.,		#T. prom. al Test
				 q=0.8,  			#Aislamiento: q
				 N_PCR=N_PCR,
				 positividad=[5.,3.,3.,3.],
				 xi_RT=0.5,		#Sensibilidad RT
				 xi_PCR=0.7		#Sensibilidad PCR
				 )
#4. Condiciones iniciales
#Dist 40% y 60%
model.var_ini(N0= 7592871.0,
			E0= 1018.685,
			IM0= 983.62, 
			IC0= 36.91, 
			IHR0= 107.91, 
			IUR0= 23.5, 
			IHD0= 9.03, 
			IUD0= 23.5,	
			D0  = 49
			)

########################################################################
#################                 PLOT                 #################
########################################################################
T_total=400
T_ini=42
model.solve(t0=T_ini,tf=T_total)
model.Test_contar()

fig = plt.figure()
ax = plt.subplot(111)	
counter=0;
#SS_data= pd.read_excel('Data_SigmaC71_MovMay19.xlsx')
SS_data= pd.read_excel(str(BBDD_Dir)+'/RES_SecSalud_Rev.xlsx',sheet_name='A06R3')
SS_time=SS_data['t']


'''PLOT'''
#nombre='Susceptibles'; Data=model.S_vect; SS=SS_data['Susceptibles']
#nombre='Requieren_HG'; Data=np.array(model.IHD_vect)+np.array(model.IR_vect)+np.array(model.IHR_vect); SS=SS_data['Requieren_HG']
nombre='Requieren_UCI'; Data=np.array(model.IUR_vect)+np.array(model.IUD_vect); SS=SS_data['Requieren_UCI']
#nombre='Sintomas_Leves'; Data=np.array(model.IM_vect)+np.array(model.IMT_vect)+np.array(model.IMNT_vect); SS=SS_data['Sintomas_Leves']
#nombre='Fallecidos'; Data=np.array(model.D_vect); SS=SS_data['Fallecidos']
#nombre='Recuperados'; Data=np.array(model.R_vect); SS=SS_data['Recuperados']
# ~ nombre=nombre1=r'Testeados PCR - $I_M$'; Data=np.array(model.PCR_IM)
# ~ nombre=nombre2=r'Testeados PCR - $I_C$'; Data2=np.array(model.PCR_IC)
# ~ nombre=nombre3='Testeados PCR - Total'; Data3=(Data+Data2)
# ~ nombre=nombre4='Testeados RT'; Data4=np.array(model.RT)
# ~ nombre='Testeados'

t_vect=np.linspace(T_ini,T_total,len(Data))

plt.title(nombre)
plt.axvline(x=92, ymin=0, ymax=1, color='grey', linestyle='--', linewidth=1.0)
ax.text(80,max(Data), str(N_PCR[0])+' tests', weight='bold', fontsize=8, rotation=90)
plt.axvline(x=122, ymin=0, ymax=1, color='grey', linestyle='--', linewidth=1.0)
ax.text(110,max(Data), str(N_PCR[1])+' tests', weight='bold', fontsize=8, rotation=90)
plt.axvline(x=153, ymin=0, ymax=1, color='grey', linestyle='--', linewidth=1.0)
ax.text(141,max(Data), str(N_PCR[2])+' tests', weight='bold', fontsize=8, rotation=90)
plt.axvline(x=183, ymin=0, ymax=1, color='grey', linestyle='--', linewidth=1.0)
ax.text(171,max(Data), str(N_PCR[3])+' tests', weight='bold', fontsize=8, rotation=90)

if 'Testeados' not in nombre:
	dif=round((max(SS)-max(Data))/max(SS)*100,2)
	plt.ylabel(u"Población", weight='bold', fontsize=10)
	ax.plot(SS_time,SS, linestyle='-', color='firebrick', label='Sin Testeo')
	ax.plot(t_vect,Data, linestyle='-', marker=' ', label=r'T. ais= 0.8, $\Delta$: '+str(dif)+r'$\%$')
else:
	plt.ylabel(u"# Tests", weight='bold', fontsize=10)
	plt.axhline(y=N_PCR[0],linestyle='--', color='grey', linewidth=1.0) 
	plt.axhline(y=N_PCR[1],linestyle='--', color='grey', linewidth=1.0) 
	plt.axhline(y=N_PCR[2],linestyle='--', color='grey', linewidth=1.0) 
	plt.axhline(y=N_PCR[3],linestyle='--', color='grey', linewidth=1.0) 
	ax.plot(t_vect,Data4, linestyle='-', marker=' ', label=str(nombre4), color='firebrick')
	ax.plot(t_vect,Data3, linestyle='-', marker=' ', label=str(nombre3), linewidth=2.0, color='green')
	ax.plot(t_vect,Data, linestyle='--', marker=' ', label=str(nombre1), color='royalblue')
	ax.plot(t_vect,Data2, linestyle='--', marker=' ', label=str(nombre2), color='orange')
	
plt.xlabel(u"Días", weight='bold', fontsize=10)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,3))
plt.yticks(size=8)
plt.xticks(size=8)
		
ax.legend(loc='best')
counter+=1
	
plt.savefig(str(nombre)+'.png',dpi = 300)
plt.show()
