# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:13:07 2022

@author: 超偉大的承儫
"""
####與ver12相比，此版本多了一些二維計算畫圖呈現的優化
import numpy as np
import cmath
import matplotlib.pyplot as plt
import copy
import os
import re
import random
import matplotlib.animation as animation
import datetime
import multiprocessing as mp
#from numba import jit
name='grap'
s=0
#net spin
E=0
#Kesoc=[]
#Kesocx=[]
#electric potential -eV/q(e)Å  in y
rdiff=1.42*1
#相鄰 unitcell的間距
dk=np.array([100,100])
#kspce切割

#relat=np.array([[2.10492,1.2152,0],[2.10492,-1.2152,0]])
dirb='o'
relat=np.array([[2.1247,1.22421,0],[2.1247,-1.22421,0]])
t=0.1
rdiffy=1.42*1
soc=0.0
t2=t*soc
ts=0
UC=0.63
UN=11.5
UB=7.8
EC=0
EB=-3.74
EN=-7.25
numberN=0
numberB=0
rcc=1.41530
windnum='o'
berrycal='o'
moisoval=0.3
#isovalue of mo
sitems=10
guess='1'
#dw for inter     up for total
if dk[0]==0 and dk[1]==0:
    folder='E={}_s={}_U={}_soc={}_mono'.format(E,s,UC,soc)
elif dk[0]!=0 and dk[1]==0:
    folder='E={}_s={}_U={}_soc={}_x'.format(E,s,UC,soc)
elif dk[1]!=0 and dk[0]==0:
    folder='E={}_s={}_U={}_soc={}_y'.format(E,s,UC,soc)   
else:
    folder='E={}_s={}_U={}_soc={}_xy'.format(E,s,UC,soc)  
           
try :
    os.mkdir('{}/'.format(name)+folder)
except:
    pass
#圖形比例
strusize=(20,20)
selcons=0.00001
if UC==0:
    newmixrate=1
    selcons=100000
else:
    newmixrate=0.3
def berrycurva(kx,ky):
    matriberry=[]
    for i in range(4):
        matriberry.append(np.zeros(shape=(int(np.round((ndw),1)),int(np.round((ndw),1))),dtype='complex128'))
    for j in range(int(np.round((ndw),1))):
        for k in range(int(np.round((ndw),1))):            
            matriberry[0][j,k]=np.dot(eigveckdw[ky][kx][:,j].conjugate(),eigveckdw[ky][kx+1][:,k])
    for j in range(int(np.round((ndw),1))):
        for k in range(int(np.round((ndw),1))):            
            matriberry[1][j,k]=np.dot(eigveckdw[ky][kx+1][:,j].conjugate(),eigveckdw[ky+1][kx+1][:,k])
    for j in range(int(np.round((ndw),1))):
        for k in range(int(np.round((ndw),1))):            
            matriberry[2][j,k]=np.dot(eigveckdw[ky+1][kx+1][:,j].conjugate(),eigveckdw[ky+1][kx][:,k])
    for j in range(int(np.round((ndw),1))):
        for k in range(int(np.round((ndw),1))):            
            matriberry[3][j,k]=np.dot(eigveckdw[ky+1][kx][:,j].conjugate(),eigveckdw[ky][kx][:,k])
    tempberry=0
    for i in range(4):
        tempberry=tempberry-cmath.log(np.linalg.det(matriberry[i])).imag
    return tempberry
def curvalist(lix,liy):
    res=np.zeros(shape=(lix.shape),dtype='complex128')
    for i in range(lix.shape[0]):
        for j in range(lix.shape[1]):
            res[i,j]=berrycurva(lix[i,j],liy[i,j])
    return res
def plotmoiso(index,locx,locy,spin,moisoval):
    plt.figure(figsize=strusize)
    plt.title('{} orbN={}(k={},{})_E={}_s={}_{}_iso={}'.format(name,index,locx,locy,str(E),str(s),spin,moisoval),fontsize=25)
    if dk[0]!=0 and dk[1]==0:
        for i in range(len(coormax)):
            if coormax[i][0]=='C':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='black',markersize=sitems)
            elif coormax[i][0]=='N':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='red',markersize=sitems)
            elif coormax[i][0]=='B':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coormax)):
            for z in range(len(coormax)):
                if ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5>0.5 and ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5<1.7:
                    plt.plot([coormax[i][1],coormax[z][1]],[coormax[i][2],coormax[z][2]],color='black')
                else:
                    pass
    elif dk[1]!=0 and dk[0]==0:
        for i in range(len(coormax)):
            if coorymax[i][0]=='C':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='black',markersize=sitems)
            elif coormax[i][0]=='N':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='red',markersize=sitems)
            elif coormax[i][0]=='B':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coorymax)):
            for z in range(len(coorymax)):
                if ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5>0.5 and ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5<1.7:
                    plt.plot([coorymax[i][1],coorymax[z][1]],[coorymax[i][2],coorymax[z][2]],color='black')
                else:
                    pass
    elif dk[1]!=0 and dk[0]!=0:
        for i in range(len(coorxy)):
            if coorxy[i][0]=='C':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='black',markersize=sitems)
            elif coorxy[i][0]=='N':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='red',markersize=sitems)
            elif coorxy[i][0]=='B':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coorxy)):
            for z in range(len(coorxy)):
                if ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5>0.5 and ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5<1.7:
                    plt.plot([coorxy[i][1],coorxy[z][1]],[coorxy[i][2],coorxy[z][2]],color='black')
                else:
                    pass
    else:
        for i in range(len(coor1)):
            if re.match('C',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='black',markersize=sitems)
            elif re.match('N',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='red',markersize=sitems)
            elif re.match('B',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coor1)):
            for z in range(len(coor1)):
                if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5>0.5 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5<1.7:
                    plt.plot([coor1[i][1],coor1[z][1]],[coor1[i][2],coor1[z][2]],color='black')
                else:
                    pass

    orden=[]
    if spin=='up':
        for i in range(len(eigveckup[locy][locx])):
            orden.append(0)
            orden[i]+=eigveckup[locy][locx][i,index]
            orden[i]=orden[i]*orden[i].conjugate()
    else:
        for i in range(len(eigveckdw[locy][locx])):
            orden.append(0)
            orden[i]+=eigveckdw[locy][locx][i,index]
            orden[i]=orden[i]*orden[i].conjugate()
    if dk[0]!=0 and dk[1]==0:
        for i in range(len(pltcx)):
            if i>len(coor1)-1:
                if spin=='up':
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckup[locy][locx][i-len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckup[locy][locx][i-len(coor1),index]),4)<cmath.pi:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='black',alpha=0.3) 
                else:
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckdw[locy][locx][i-len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckdw[locy][locx][i-len(coor1),index]),4)<cmath.pi:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='blue',alpha=0.3)
                    else:
                        plt.scatter(pltcx[i],pltcy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='black',alpha=0.3) 

            else:
                if spin=='up':
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckup[locy][locx][i,index]),4)>=0 and np.round(cmath.phase(eigveckup[locy][locx][i,index]),4)<cmath.pi :
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='black',alpha=0.3)
                else:
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckdw[locy][locx][i,index]),4)>=0 and np.round(cmath.phase(eigveckdw[locy][locx][i,index]),4)<cmath.pi :
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcx[i],pltcy[i],s=float(orden[i]/moisoval*4000),c='black',alpha=0.3)


    elif dk[1]!=0 and dk[0]==0:
        for i in range(len(pltcxy)):
            if i>len(coor1)-1:
                if spin=='up':
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckup[locy][locx][i-len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckup[locy][locx][i-len(coor1),index]),4)<cmath.pi:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='black',alpha=0.3)
                else:
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckdw[locy][locx][i-len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckdw[locy][locx][i-len(coor1),index]),4)<cmath.pi:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='blue',alpha=0.3)
                    else:
                        plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i-len(coor1)]/moisoval*4000),c='black',alpha=0.3)
            else:
                if spin=='up':
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckup[locy][locx][i,index]),4)>=0 and np.round(cmath.phase(eigveckup[locy][locx][i,index]),4)<cmath.pi :
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='black',alpha=0.3) 
                else:
                    if locx==0 and locy==0:
                        if np.round(cmath.phase(eigveckdw[locy][locx][i,index]),4)>=0 and np.round(cmath.phase(eigveckdw[locy][locx][i,index]),4)<cmath.pi :
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                        else:
                            plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3) 
                    else:
                        plt.scatter(pltcxy[i],pltcyy[i],s=float(orden[i]/moisoval*4000),c='black',alpha=0.3) 
    elif dk[0]==0 and dk[1]==0:
        for i in range(len(pltcxmono)):
            if spin=='up':
                if np.round((eigveckup[locy][locx][i,index]).real,4)>0:
                    plt.scatter(pltcxmono[i],pltcymono[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                else:
                    plt.scatter(pltcxmono[i],pltcymono[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3)  
            else:
                if np.round((eigveckdw[locy][locx][i,index]).real,4)>0 :
                    plt.scatter(pltcxmono[i],pltcymono[i],s=float(orden[i]/moisoval*4000),c='red',alpha=0.3)    
                else:
                    plt.scatter(pltcxmono[i],pltcymono[i],s=float(orden[i]/moisoval*4000),c='blue',alpha=0.3)         
    else:
        for i in range(len(coorxy)):
            if spin=='up':
                if locx==0 and locy==0:
                    if np.round(cmath.phase(eigveckup[locy][locx][i%len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckup[locy][locx][i%len(coor1),index]),4)<cmath.pi :
                        plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                    else:
                        plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='blue',alpha=0.3) 
                else:
                    plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='black',alpha=0.3) 
            else:
                if locx==0 and locy==0:
                    if np.round(cmath.phase(eigveckdw[locy][locx][i%len(coor1),index]),4)>=0 and np.round(cmath.phase(eigveckdw[locy][locx][i%len(coor1),index]),4)<cmath.pi :
                        plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='red',alpha=0.3)    
                    else:
                        plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='blue',alpha=0.3) 
                else:
                    plt.scatter(coorxy[i][1],coorxy[i][2],s=float(orden[i%len(coor1)]/moisoval*4000),c='black',alpha=0.3) 
    plt.axis('off')   
    plt.savefig('{}/'.format(name)+folder+'/{} orbN={}(k={},{})_E={}_s={}_{}_iso={}.jpg'.format(name,index,locx,locy,str(E),str(s),spin,moisoval)) 
    plt.show()
def bandorb1D():
    if dk[0]!=0 and dk[1]==0:
        plxtemp=np.linspace(0,2*cmath.pi,dk[0])
        plxdw=[]
        plxup=[]
        for i in range(dk[0]):
            for k in range(len(coor1)):
                plxdw.append(plxtemp[i])
                plxup.append(plxtemp[i])
        plyup=[]
        for i in range(dk[0]):
            for z in range(len(eigekup[0][i])):
                plyup.append(eigekup[0][i][z]-fermiup)
        plydw=[]
        for i in range(dk[0]):
            for z in range(len(eigekdw[0][i])):
                plydw.append(eigekdw[0][i][z]-fermidw)
        plt.figure(figsize=(5,15))
        plt.title('{} bandx of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,2*cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
        plt.scatter(plxup,plyup,s=0.1,c='black')
        plt.xlabel('k(1/a)')
        plt.ylabel('energy (eV)')
        plt.xlim(0,2*cmath.pi)
        loax=list(np.arange(0,2.5*cmath.pi,0.5*cmath.pi))
        xst=['0','0.5π','π','1.5π','2π']
        plt.xticks(loax,xst)
        plt.ylim(-4,4)
        plt.savefig('{}/'.format(name)+folder+'/{}spin up band_E={}_s={}.jpg'.format(name,str(E),str(s)))
        
        
        plt.figure(figsize=(5,15))
        plt.title('{} bandx of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
        plt.scatter(plxdw,plydw,s=0.1,c='black')
        plt.xlabel('k(1/a)')
        plt.ylabel('energy (eV)')
        plt.xlim(0,cmath.pi)
        loax=list(np.arange(0,1.5*cmath.pi,0.5*cmath.pi))
        xst=['0','0.5π','π']
        plt.xticks(loax,xst)
        plt.ylim(-4,4)
        plt.savefig('{}/'.format(name)+folder+'/{}spin dw band_E={}_s={}.jpg'.format(name,str(E),str(s)))
    elif dk[1]!=0 and dk[0]==0:
        plxtemp=np.linspace(0,2*cmath.pi,dk[1])
        plxdw=[]
        plxup=[]
        for i in range(dk[1]):
            for k in range(len(coor1)):
                plxdw.append(plxtemp[i])
                plxup.append(plxtemp[i])
        plyup=[]
        for i in range(dk[1]):
            for z in range(len(eigekup[i][0])):
                plyup.append(eigekup[i][0][z]-fermiup)
        plydw=[]
        for i in range(dk[1]):
            for z in range(len(eigekdw[i][0])):
                plydw.append(eigekdw[i][0][z]-fermidw)
        plt.figure(figsize=(5,15))
        plt.title('{} bandy of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
        plt.scatter(plxup,plyup,s=0.1,c='black')
        plt.xlabel('k(1/a)')
        plt.ylabel('energy (eV)')
        plt.xlim(0,cmath.pi)
        loax=list(np.arange(0,1.5*cmath.pi,0.5*cmath.pi))
        xst=['0','0.5π','π']
        plt.xticks(loax,xst)
        plt.ylim(-1,1)
        plt.savefig('{}/'.format(name)+folder+'/{}spin up band_E={}_s={}.jpg'.format(name,str(E),str(s)))
        
        
        plt.figure(figsize=(5,15))
        plt.title('{} bandy of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
        plt.scatter(plxdw,plydw,s=0.1,c='black')
        plt.xlabel('k(1/a)')
        plt.ylabel('energy (eV)')
        plt.xlim(0,cmath.pi)
        loax=list(np.arange(0,1.5*cmath.pi,0.5*cmath.pi))
        xst=['0','0.5π','π']
        plt.xticks(loax,xst)
        plt.ylim(-4,4)
        plt.savefig('{}/'.format(name)+folder+'/{}spin dw band_E={}_s={}.jpg'.format(name,str(E),str(s)))
    elif dk[1]==0 and dk[0]==0:
        plydw=[]
        for i in range(len(eigekdw[0][0])):
            plydw.append(eigekdw[0][0][i]-fermidw)
        plyup=[]
        for i in range(len(eigekup[0][0])):
            plyup.append(eigekup[0][0][i]-fermiup)
        plt.figure(figsize=(5,15))
        plt.title('{} orbital of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
    
    
        for i in range(len(plyup)):
            plt.axhline(y=plyup[i],xmin=0.25,xmax=0.75,c='black')
    
        plt.ylabel('energy (eV)')
        plt.xlim(0,cmath.pi)
    
        plt.xticks([],[])
        plt.ylim(-1,1)
        plt.savefig('{}/'.format(name)+folder+'/{}spin up band_E={}_s={}.jpg'.format(name,str(E),str(s)))
        
        
        plt.figure(figsize=(5,15))
        plt.title('{} orbital of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
        plt.plot([0,cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
        for i in range(len(plydw)):
            plt.axhline(y=plydw[i],xmin=0.25,xmax=0.75,c='black')
    
        plt.ylabel('energy (eV)')
        plt.xlim(0,cmath.pi)
    
        plt.xticks([],[])
        plt.ylim(-4,4)
        plt.savefig('{}/'.format(name)+folder+'/{}spin dw band_E={}_s={}.jpg'.format(name,str(E),str(s)))       
    else:
        pass  
def band2D():
    plxtempx=np.linspace(0,cmath.pi,int((dk[0])/2))
    plxtempy=np.linspace(0,cmath.pi,int((dk[1])/2))
    plx=[]
    for i in range(int((dk[1])/2)):
        for z in range(len(coor1)):
            plx.append(plxtempy[i])
    for i in range(int((dk[0])/2)):
        for z in range(len(coor1)):
            plx.append(plxtempx[i]+cmath.pi)
    for i in range(int((dk[1])/2)):
        for z in range(len(coor1)):
            plx.append(plxtempy[i]+2*cmath.pi)
    plyup=[]
    for i in range(int((dk[1])/2)):
        for z in range(len(coor1)):
            plyup.append(eigekup[int((dk[1]-1)/2)-i][0][z]-fermiup)
    for i in range(int((dk[0])/2)):
        for z in range(len(coor1)):
            plyup.append(eigekup[0][i][z]-fermiup)
    for i in range(int((dk[1])/2)):
        for z in range(len(eigekup[i][int((dk[0]-1)/2)])):
            plyup.append(eigekup[i][int((dk[0]-1)/2)][z]-fermiup)            
    plydw=[]
    for i in range(int((dk[1])/2)):
        for z in range(len(coor1)):
            plydw.append(eigekdw[int((dk[1]-1)/2)-i][0][z]-fermidw)
    for i in range(int((dk[0])/2)):
        for z in range(len(coor1)):
            plydw.append(eigekdw[0][i][z]-fermidw)
    for i in range(int((dk[1])/2)):
        for z in range(len(coor1)):
            plydw.append(eigekdw[i][int((dk[0]-1)/2)][z]-fermidw)   
    plt.figure(figsize=(5,15))
    plt.title('{} bandy of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
    plt.plot([0,3*cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
    plt.axvline(x=cmath.pi, ymin=-4, ymax=4)
    plt.axvline(x=2*cmath.pi, ymin=-4, ymax=4)
    plt.scatter(plx,plyup,s=0.1,c='black')
    plt.xlabel('k(1/a)')
    plt.ylabel('energy (eV)')
    plt.xlim(0,3*cmath.pi)
    loax=list(np.arange(0,4*cmath.pi,1*cmath.pi))
    xst=['Y','Gamma','X','M']
    plt.xticks(loax,xst)
    plt.ylim(-6,6)
    plt.savefig('{}/'.format(name)+folder+'/{}spin up band_E={}_s={}.jpg'.format(name,str(E),str(s)))
        
    plt.figure(figsize=(5,15))
    plt.title('{} bandy of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=15)
    plt.plot([0,3*cmath.pi],[0,0],marker='.',color='red', linewidth=2, linestyle=":")
    plt.axvline(x=cmath.pi, ymin=-4, ymax=4)
    plt.axvline(x=2*cmath.pi, ymin=-4, ymax=4)
    plt.scatter(plx,plydw,s=0.1,c='black')
    plt.xlabel('k(1/a)')
    plt.ylabel('energy (eV)')
    plt.xlim(0,3*cmath.pi)
    loax=list(np.arange(0,4*cmath.pi,1*cmath.pi))
    xst=['Y','Gamma','X','M']
    plt.xticks(loax,xst)
    plt.ylim(-4,4)
    plt.savefig('{}/'.format(name)+folder+'/{}spin dw band_E={}_s={}.jpg'.format(name,str(E),str(s)))
def cauwind2D():
    fig=plt.figure(figsize=(15,15))
    ax = fig.gca()
    ax.set_title('{} winding number X per Y of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    windpltxup=[]
    windpltyup=[]
    for i in range(dk[1]):
        windpltxup.append([])
        windpltyup.append([])
    for i in range(len(windup)):
        windpltxup[i//dk[1]].append(windup[i].real/10**3)
        windpltyup[i//dk[1]].append(windup[i].imag/10**3)
    winvec=ax.scatter(windpltxup[0],windpltyup[0],s=100,color='blue')
    #def init():
    #    ax.scatter(windpltxup[0],windpltyup[0],s=100)
    def upfun(i):
        #winvec=ax.scatter(windpltxup[i],windpltyup[i],color='blue',s=100)
        winvec.set_offsets(np.c_[windpltxup[i],windpltyup[i]])
    ani=animation.FuncAnimation(fig,upfun,frames=len(windpltxup),interval=1)
    ax.set_xlim(-20,10)
    ax.set_ylim(-15,15)
    ax.scatter(0,0,s=100,c='red')
    ax.set_xlabel('real',fontsize=25)
    ax.set_ylabel('image',fontsize=25)
    ani.save('{}/'.format(name)+folder+'/{}X per Y_upwindanimation.gif'.format(name),fps=40)
    fig.savefig('{}/'.format(name)+folder+'/{} X per Y winding number of spin up_E={}_s={}.jpg'.format(name,str(E),str(s)))
    fig=plt.figure(figsize=(15,15))
    ax = fig.gca()
    ax.set_title('{} winding number Y per X of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    windpltxup=[]
    windpltyup=[]
    for i in range(dk[0]):
        windpltxup.append([])
        windpltyup.append([])
    for i in range(len(windup)):
        windpltxup[i%dk[0]].append(windup[i].real/10**3)
        windpltyup[i%dk[0]].append(windup[i].imag/10**3)
    winvec2=ax.scatter(windpltxup[0],windpltyup[0],s=100,color='blue')
    #def init():
    #    ax.scatter(windpltxup[0],windpltyup[0],s=100)
    def upfun2(i):
        #winvec=ax.scatter(windpltxup[i],windpltyup[i],color='blue',s=100)
        winvec2.set_offsets(np.c_[windpltxup[i],windpltyup[i]])
    ani=animation.FuncAnimation(fig,upfun2,frames=len(windpltxup),interval=1)
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    ax.scatter(0,0,s=100,c='red')
    ax.set_xlabel('real',fontsize=25)
    ax.set_ylabel('image',fontsize=25)
    ani.save('{}/'.format(name)+folder+'/{}Y per X_upwindanimation.gif'.format(name),fps=40)
    fig.savefig('{}/'.format(name)+folder+'/{} Y per X winding number of spin up_E={}_s={}.jpg'.format(name,str(E),str(s)))
def cauwind1D():    
    fig=plt.figure(figsize=(15,15))
    ax = fig.gca()
    if dk[0]!=0 and dk[1]==0:
        ax.set_title('{} winding numberX of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    elif dk[1]!=0 and dk[0]==0:
        ax.set_title('{} winding numberY of spin up_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    windpltxup=[]
    windpltyup=[]
    winvec,=ax.plot([],[],color='blue',linewidth=3)
    for i in range(len(windup)):
        windpltxup.append(windup[i].real)
        windpltyup.append(windup[i].imag)
    def upfun(i):
        wx=[windpltxup[i],0]
        wy=[windpltyup[i],0]
        winvec.set_data(wx,wy)
    ani=animation.FuncAnimation(fig=fig,func=upfun,frames=len(windpltxup),interval=1)
    ax.scatter(0,0,s=100,c='red')
    ax.scatter(windpltxup,windpltyup,s=10,c='black')
    ax.set_xlabel('real',fontsize=25)
    ax.set_ylabel('image',fontsize=25)
    ani.save('{}/'.format(name)+folder+'/{}_upwindanimation.gif'.format(name),fps=40)
    fig.savefig('{}/'.format(name)+folder+'/{}winding number of spin up_E={}_s={}.jpg'.format(name,str(E),str(s)))
    
    
    fig=plt.figure(figsize=(15,15))
    ax = fig.gca()
    if dk[0]!=0 and dk[1]==0:
        ax.set_title('{} winding numberX of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    elif dk[1]!=0 and dk[0]==0:
        ax.set_title('{} winding numberY of spin dw_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    windpltxdw=[]
    windpltydw=[]
    #winvec,=ax.plot([],[],color='blue',linewidth=3)
    wnumber=0
    for i in range(len(winddw)):
        windpltxdw.append(winddw[i].real)
        windpltydw.append(winddw[i].imag)
        if i==0:
            pass
        else:
            if cmath.phase(winddw[i])-cmath.phase(winddw[i-1])>0.1 or cmath.phase(winddw[i])-cmath.phase(winddw[i-1])<-0.1:
                pass
            else:
                wnumber=wnumber+(cmath.phase(winddw[i])-cmath.phase(winddw[i-1]))
    print('wind number={}'.format(np.round(wnumber/(2*cmath.pi),1)))
    #def upfun(i):
    #    wx=[windpltxdw[i],0]
    #    wy=[windpltydw[i],0]
    #    winvec.set_data(wx,wy)
    #    return winvec,
    #ani=animation.FuncAnimation(fig=fig,func=upfun,frames=len(windpltxdw),blit=True,interval=1)
    ax.scatter(0,0,s=100,c='red')
    ax.scatter(windpltxdw,windpltydw,s=10,c='black')
    ax.set_xlabel('real',fontsize=25)
    ax.set_ylabel('image',fontsize=25)
    #ani.save('{}/{}_dwwindanimation.gif'.format(name,name), fps=40)
    fig.savefig('{}/'.format(name)+folder+'/{}winding number of spin dw_E={}_s={}.jpg'.format(name,str(E),str(s)))
def berryca(kpoint):
    if kpoint[1]==0:
        berrydata=open('{}/'.format(name)+folder+'/{}_berryphaseX_E={}_s={}.txt'.format(name,E,s),'w+')
    
        matriberry=[]
        for i in range(dk[0]-2):
            matriberry.append(np.zeros(shape=(int(np.round((ndw),1)),int(np.round((ndw),1))),dtype='complex128'))        
            for j in range(int(np.round((ndw),1))):
                for k in range(int(np.round((ndw),1))):            
                    matriberry[i][j,k]=np.dot(eigveckdw[0][i][:,j].conjugate(),eigveckdw[0][i+1][:,k])
                    
                    
                       
        matriberry.append(np.zeros(shape=(int(np.round((ndw),1)),int(np.round((ndw),1))),dtype='complex128'))    
        for j in range(int(np.round((ndw),1))):
           for k in range(int(np.round((ndw),1))):            
                matriberry[-1][j,k]=np.dot(eigveckdw[0][-2][:,j].conjugate(),eigveckdw[0][0][:,k])
        
        tempberry=0
        for i in range(dk[0]-1):
            tempberry=tempberry-cmath.log(np.linalg.det(matriberry[i])).imag
        print('Z2 = {} spin dw'.format(np.round((tempberry/cmath.pi),4)%2))
        #print(np.round((tempberry/cmath.pi),4))
        berrydata.write('Z2 = {} spin dw \n'.format(np.round((tempberry/cmath.pi),4)%2))
        matriberry=[]
        for i in range(dk[0]-2):
            matriberry.append(np.zeros(shape=(int(np.round((nup),1)),int(np.round((nup),1))),dtype='complex128'))        
            for j in range(int(np.round((nup),1))):
                for k in range(int(np.round((nup),1))):            
                    matriberry[i][j,k]=np.dot(eigveckup[0][i][:,j].conjugate(),eigveckup[0][i+1][:,k])   
        matriberry.append(np.zeros(shape=(int(np.round((nup),1)),int(np.round((nup),1))),dtype='complex128'))    
        for j in range(int(np.round((nup),1))):
           for k in range(int(np.round((nup),1))):            
               matriberry[-1][j,k]=np.dot(eigveckup[0][-2][:,j].conjugate(),eigveckup[0][0][:,k])
        tempberry=0
        for i in range(dk[0]-1):
            tempberry=tempberry-cmath.log(np.linalg.det(matriberry[i])).imag
        print('Z2 = {} spin up'.format(np.round(tempberry/cmath.pi,4)%2))
        #print(np.round((tempberry/cmath.pi),4))
        berrydata.write('Z2 = {} spin up \n'.format(np.round(tempberry/cmath.pi,4)%2))
        berrydata.close()
    elif kpoint[0]==0:
        berrydata=open('{}/'.format(name)+folder+'/{}_berryphaseY_E={}_s={}.txt'.format(name,E,s),'w+')
    
        matriberry=[]
        for i in range(dk[1]-2):
            matriberry.append(np.zeros(shape=(int(np.round((ndw),1)),int(np.round((ndw),1))),dtype='complex128'))        
            for j in range(int(np.round((ndw),1))):
                for k in range(int(np.round((ndw),1))):            
                    matriberry[i][j,k]=np.dot(eigveckdw[i][0][:,j].conjugate(),eigveckdw[i+1][0][:,k])
                    
                    
                       
        matriberry.append(np.zeros(shape=(int(np.round((ndw),1)),int(np.round((ndw),1))),dtype='complex128'))    
        for j in range(int(np.round((ndw),1))):
           for k in range(int(np.round((ndw),1))):            
                matriberry[-1][j,k]=np.dot(eigveckdw[-2][0][:,j].conjugate(),eigveckdw[0][0][:,k])
        
        tempberry=0
        for i in range(dk[1]-1):
            tempberry=tempberry-cmath.log(np.linalg.det(matriberry[i])).imag
        print('Z2 = {} spin dw'.format(np.round((tempberry/cmath.pi),4)%2))
        berrydata.write('Z2 = {} spin dw \n'.format(np.round((tempberry/cmath.pi),4)%2))
        matriberry=[]
        for i in range(dk[1]-2):
            matriberry.append(np.zeros(shape=(int(np.round((nup),1)),int(np.round((nup),1))),dtype='complex128'))        
            for j in range(int(np.round((nup),1))):
                for k in range(int(np.round((nup),1))):            
                    matriberry[i][j,k]=np.dot(eigveckup[i][0][:,j].conjugate(),eigveckup[i+1][0][:,k])   
        matriberry.append(np.zeros(shape=(int(np.round((nup),1)),int(np.round((nup),1))),dtype='complex128'))    
        for j in range(int(np.round((nup),1))):
           for k in range(int(np.round((nup),1))):            
               matriberry[-1][j,k]=np.dot(eigveckup[-2][0][:,j].conjugate(),eigveckup[0][0][:,k])
        tempberry=0
        for i in range(dk[1]-1):
            tempberry=tempberry-cmath.log(np.linalg.det(matriberry[i])).imag
        print('Z2 = {} spin up'.format(np.round(tempberry/cmath.pi,4)%2))
        berrydata.write('Z2 = {} spin up \n'.format(np.round(tempberry/cmath.pi,4)%2))
        berrydata.close()
def hmatrix(kpoint):
    matxup=np.zeros(shape=(len(coor1),len(coor1)),dtype='complex128')
    for i in range(len(coor1)):
        if coor1[i] in lib_a:
            privec=priveca
        else:
           privec=privecb
        if dk[0]==0 and dk[1]==0:
            for z in range(len(coor1)):
                if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                    matxup[i,z]+=t
                elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.40:
                    matxup[i,z]+=ts
                    if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                elif  i==z:
                    if re.match('C',coor1[i][0]):
                        matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EC
                    elif re.match('N',coor1[i][0]):
                        matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EN
                    elif re.match('B',coor1[i][0]):
                        matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EB
                    else:
                        pass
                else:
                    pass
        elif dk[0]!=0 and dk[1]==0:
            if coor1[i][1]>(max(coorx)+min(coorx))/2 :
                for z in range(len(coormax)):
                    if ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxup[i,z%len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[0]))
                        else:
                            matxup[i,z]+=t   
                    elif ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxup[i,z%len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            if (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                        else:
                            matxup[i,z]+=ts
                            if (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass   
            elif coor1[i][1]<(max(coorx)+min(coorx))/2 :
                for z in range(len(coormin)):
                    if ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                        else:
                            matxup[i,z]+=t
                    elif ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            if (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                        else:
                            matxup[i,z]+=ts
                            if (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass   
            else:
                for z in range(len(coor1)):
                    if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                        matxup[i,z]+=t
                    elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.43:
                        matxup[i,z]+=ts
                        if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass
        elif dk[1]!=0 and dk[0]==0:
            if coor1[i][2]>(max(coory)+min(coory))/2 :
                for z in range(len(coorymax)):
                    if ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[1]))
                        else:
                            matxup[i,z]+=t   
                    elif ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            if (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                        else:
                            matxup[i,z]+=ts
                            if (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass   
            elif coor1[i][2]<(max(coory)+min(coory))/2 :
                for z in range(len(coorymin)):
                    if ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                        else:
                            matxup[i,z]+=t
                    elif ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxup[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            if (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z-len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                        else:
                            matxup[i,z]+=ts
                            if (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass   
            else:
                for z in range(len(coor1)):
                    if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                        matxup[i,z]+=t
                    elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.43:
                        matxup[i,z]+=ts
                        if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxup[i,z]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxup[z,z]+=UC*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxup[z,z]+=UN*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxup[z,z]+=UB*(densitydw[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass
        else:
            for z in range(len(coorxy)):
                rlaccon=[0,0]
                if z//len(coor1)==8:
                    rlaccon=[1,1]
                    
                elif z//len(coor1)==7:
                    rlaccon=[0,1]
                    
                elif z//len(coor1)==6:
                    rlaccon=[-1,1]
                    
                elif z//len(coor1)==5:
                    rlaccon=[1,0]
                    
                elif z//len(coor1)==4:
                    rlaccon=[0,0]
                    
                elif z//len(coor1)==3:
                    rlaccon=[-1,0]
                    
                elif z//len(coor1)==2:
                    rlaccon=[1,-1]
                    
                elif z//len(coor1)==1:
                    rlaccon=[0,-1]
                    
                elif z//len(coor1)==0:   
                    rlaccon=[-1,-1]    
                    
                 
                if ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5>0.3:
                    matxup[i,z%len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                elif ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5>2.43:
                    matxup[i,z%len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    if (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxup[i,z%len(coor1)]+=1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                elif  i==z:
                    if re.match('C',coor1[i][0]):
                        matxup[z,z]+=UC*(densitydw[i])+EC
                    elif re.match('N',coor1[i][0]):
                        matxup[z,z]+=UN*(densitydw[i])+EN
                    elif re.match('B',coor1[i][0]):
                        matxup[z,z]+=UB*(densitydw[i])+EB
                    else:
                        pass
                else:
                    pass       
    matxdw=np.zeros(shape=(len(coor1),len(coor1)),dtype='complex128')
    for i in range(len(coor1)):
        if coor1[i] in lib_a:
            privec=priveca
        else:
            privec=privecb
        if dk[0]==0 and dk[1]==0:
            for z in range(len(coor1)):
                if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                    matxdw[i,z]+=t
                elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.43:
                    matxdw[i,z]+=ts
                    if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                    elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                elif  i==z:
                    if re.match('C',coor1[i][0]):
                        matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EC
                    elif re.match('N',coor1[i][0]):
                        matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EN
                    elif re.match('B',coor1[i][0]):
                        matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EB
                    else:
                        pass
                else:
                    pass
        elif dk[0]!=0 and dk[1]==0:
            if coor1[i][1]>(max(coorx)+min(coorx))/2 :
                for z in range(len(coormax)):
                    if ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[0]))
                        else:
                            matxdw[i,z]+=t   
                    elif ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coormax[z][1])**2+(coor1[i][2]-coormax[z][2])**2+(coor1[i][3]-coormax[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            if (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]))
                        else:
                            matxdw[i,z]+=ts
                            if (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormax[z][1],coor1[i][2]-coormax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass   
            elif coor1[i][1]<(max(coorx)+min(coorx))/2 :
                for z in range(len(coormin)):
                    if ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                        else:
                            matxdw[i,z]+=t
                    elif ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coormin[z][1])**2+(coor1[i][2]-coormin[z][2])**2+(coor1[i][3]-coormin[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            if (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[0]))
                        else:
                            matxdw[i,z]+=ts
                            if (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass   
            else:
                for z in range(len(coor1)):
                    if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                        matxdw[i,z]+=t
                    elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.43:
                        matxdw[i,z]+=ts
                        if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coory)-coor1[i][2])+EB
                        else:
                            pass
                    else:
                        pass
        elif dk[1]!=0 and dk[0]==0:
            if coor1[i][2]>(max(coory)+min(coory))/2 :
                for z in range(len(coorymax)):
                    if ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[1]))
                        else:
                            matxdw[i,z]+=t   
                    elif ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorymax[z][1])**2+(coor1[i][2]-coorymax[z][2])**2+(coor1[i][3]-coorymax[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            if (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[1]))
                        else:
                            matxdw[i,z]+=ts
                            if (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymax[z][1],coor1[i][2]-coorymax[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass   
            elif coor1[i][2]<(max(coory)+min(coory))/2 :
                for z in range(len(coorymin)):
                    if ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5>0.3:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                        else:
                            matxdw[i,z]+=t
                    elif ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorymin[z][1])**2+(coor1[i][2]-coorymin[z][2])**2+(coor1[i][3]-coorymin[z][3])**2)**0.5>2.43:
                        if z>len(coor1)-1:
                            matxdw[i,z-len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            if (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z-len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,-relat[1]))
                        else:
                            matxdw[i,z]+=ts
                            if (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                            elif (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorymin[z][1],coor1[i][2]-coorymin[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                                matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass   
            else:
                for z in range(len(coor1)):
                    if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>0.3:
                        matxdw[i,z]+=t
                    elif ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2+(coor1[i][3]-coor1[z][3])**2)**0.5>2.43:
                        matxdw[i,z]+=ts
                        if (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coormin[z][1],coor1[i][2]-coormin[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])
                        elif (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coor1[z][1],coor1[i][2]-coor1[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                            matxdw[i,z]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2]) 
                    elif  i==z:
                        if re.match('C',coor1[i][0]):
                            matxdw[z,z]+=UC*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EC
                        elif re.match('N',coor1[i][0]):
                            matxdw[z,z]+=UN*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EN
                        elif re.match('B',coor1[i][0]):
                            matxdw[z,z]+=UB*(densityup[i])+E*abs(max(coorx)-coor1[i][1])+EB
                        else:
                            pass
                    else:
                        pass
        else:
            for z in range(len(coorxy)):
                rlaccon=[0,0]
                if z//len(coor1)==8:
                    rlaccon=[1,1]
                    
                elif z//len(coor1)==7:
                    rlaccon=[0,1]
                    
                elif z//len(coor1)==6:
                    rlaccon=[-1,1]
                    
                elif z//len(coor1)==5:
                    rlaccon=[1,0]
                    
                elif z//len(coor1)==4:
                    rlaccon=[0,0]
                    
                elif z//len(coor1)==3:
                    rlaccon=[-1,0]
                    
                elif z//len(coor1)==2:
                    rlaccon=[1,-1]
                    
                elif z//len(coor1)==1:
                    rlaccon=[0,-1]
                    
                elif z//len(coor1)==0:   
                    rlaccon=[-1,-1]    
                    
                
                if ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5<rcc+0.3 and ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5>0.3:
                    matxdw[i,z%len(coor1)]+=t*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                elif ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5<2.50 and ((coor1[i][1]-coorxy[z][1])**2+(coor1[i][2]-coorxy[z][2])**2+(coor1[i][3]-coorxy[z][3])**2)**0.5>2.43:
                    matxdw[i,z%len(coor1)]+=ts*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    if (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[0]-privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[0]-privec[1])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[0]+privec[1])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[0]+privec[1])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[0]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[0]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[0],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[0]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[0]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[0])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(privec[1]-privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(privec[1]-privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[1],-privec[2])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                    elif (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])<=(-privec[1]+privec[2])+np.array([0.5,0.5,0])).all() and (np.array([coor1[i][1]-coorxy[z][1],coor1[i][2]-coorxy[z][2],0])>=(-privec[1]+privec[2])-np.array([0.5,0.5,0])).all():
                        matxdw[i,z%len(coor1)]+=-1j*t2*np.sign(np.cross(privec[2],-privec[1])[2])*cmath.exp(1j*np.dot(kpoint,relat[0]*rlaccon[0]+relat[1]*rlaccon[1]))
                elif  i==z:
                    if re.match('C',coor1[i][0]):
                        matxdw[z,z]+=UC*(densityup[i])+EC
                    elif re.match('N',coor1[i][0]):
                        matxdw[z,z]+=UN*(densityup[i])+EN
                    elif re.match('B',coor1[i][0]):
                        matxdw[z,z]+=UB*(densityup[i])+EB
                    else:
                        pass
                else:
                    pass   
    return matxup,matxdw





#@jit(nopython=True)
def multcore(mkx,mky,eigekup,eigekdw,eigveckup,eigveckdw,windup,winddw,kspacex,kspacey):
    matxup,matxdw=hmatrix(kspacex[mkx]+kspacey[mky])      
    eigeup,eigvecup=np.linalg.eig(matxup)
    eigedw,eigvecdw=np.linalg.eig(matxdw)
    eigekup[mky,mkx]=copy.deepcopy(eigeup)
    eigekdw[mky,mkx]=copy.deepcopy(eigedw)
    eigveckup[mky,mkx]=copy.deepcopy(eigvecup)
    eigveckdw[mky,mkx]=copy.deepcopy(eigvecdw)
    if windnum=='on':
        tempmatup=np.zeros(shape=(len(lib_a),len(lib_a)),dtype='complex128')
        tempmatdw=np.zeros(shape=(len(lib_a),len(lib_a)),dtype='complex128')
        for i in range(len(lib_a)):
            for j in range(len(lib_a)):
                if matxup[i][j]==0 :
                    pass
                else:
                    print('! not chiral symmetry')
        for i in range(len(lib_a)):
            for j in range(len(lib_b)):
                tempmatup[i][j]=matxup[i+len(lib_a)][j]
                tempmatdw[i][j]=matxdw[i+len(lib_a)][j]
        windup[mky*dk[0]+mkx]=np.linalg.det(tempmatup)
        winddw[mky*dk[0]+mkx]=np.linalg.det(tempmatdw)
    else:
        pass
com=open('{}/{}.com'.format(name,name),'r+')
commax=open('{}/{}_onlyC.com'.format(name,name),'w+')
#commax1=open('{}max1.com'.format(name),'w+')
densityrecord=open('{}/{}density_'.format(name,name)+folder+'.txt','a+')
densityrecord.close()
densityrecordread=open('{}/{}density_'.format(name,name)+folder+'.txt','r+')
denrec=densityrecordread.readlines()
coor=com.readlines()
densityrecordread.close()
coor1=[]



#1:random

for i in range(len(coor)):
    coor2=coor[i].split() 
    coor1.append(coor2)
i=0
while i<len(coor1):
    if len(coor1[i])<1:
        coor1.pop(i)
    else:
        if re.match('C',coor1[i][0]) or re.match('B',coor1[i][0]) or re.match('N',coor1[i][0]):
            i=i+1
        else:
            coor1.pop(i)
            
for i in range(len(coor1)):
    coor1[i][1]=np.round(float(coor1[i][1]),4)
    coor1[i][2]=np.round(float(coor1[i][2]),4)
    coor1[i][3]=np.round(float(coor1[i][3]),4)
    
lib_b=[]
lib_b.append(coor1[0])
countb=0
while True:      
    if countb+1>len(lib_b):
        break
    else:            
        for j in range(len(coor1)):
            if ((coor1[j][1]-lib_b[countb][1])**2+(coor1[j][2]-lib_b[countb][2])**2)**0.5>2.44 and ((coor1[j][1]-lib_b[countb][1])**2+(coor1[j][2]-lib_b[countb][2])**2)**0.5<2.51:
                if coor1[j] not in lib_b:
                    lib_b.append(coor1[j])
                else:
                    pass
            else:
                pass
        countb=countb+1
lib_bnum=len(lib_b)
lib_anum=len(coor)-lib_bnum
lib_a=[]
for i in range(len(coor1)):
    if coor1[i] not in lib_b:
        lib_a.append(coor1[i])
    else:
        pass
coor2=copy.deepcopy(coor1)
coor1=[]
for i in range(len(lib_b)):
    coor1.append(lib_b[i])    
for i in range(len(lib_a)):
    coor1.append(lib_a[i])
coorx=[]
for i in range(len(coor1)):
    coorx.append(coor1[i][1])
coory=[]
for i in range(len(coor1)):
    coory.append(coor1[i][2])
rlac=np.round(max(coorx)-min(coorx)+rdiff,4)
rlacy=np.round(max(coory)-min(coory)+rdiffy,4)
bzlat=[]
if dk[0]==0 and dk[1]==0:
    relat=np.array([[0,0,0],[0,0,0]])
    bzlat.append(np.array([0,0,0]))
    bzlat.append(np.array([0,0,0]))
elif dk[0]!=0 and dk[1]==0:
    if dirb=='on':
        relat=np.array([[rlac,0,0],[0,0,0]])
    else:
        pass
    bzlat.append(2*cmath.pi*relat[0]/(np.linalg.norm(relat[0])**2))
    bzlat.append(np.array([0,0,0]))
elif dk[1]!=0 and dk[0]==0:
    if dirb=='on':
        relat=np.array([[0,0,0],[0,rlacy,0]])
    else:
        pass
    bzlat.append(np.array([0,0,0]))
    bzlat.append(2*cmath.pi*relat[1]/(np.linalg.norm(relat[1])**2))
    
else:
    if dirb=='on':
        relat=np.array([[rlac,0,0],[0,rlacy,0]])
    else:
        pass
    bzlat.append(2*cmath.pi*np.cross(relat[1],np.array([0,0,1]))/np.linalg.norm(np.cross(relat[0],relat[1])))
    bzlat.append(2*cmath.pi*np.cross(np.array([0,0,1]),relat[0])/np.linalg.norm(np.cross(relat[0],relat[1])))

coormax=copy.deepcopy(coor1)
coormin=copy.deepcopy(coor1)
cooa=copy.deepcopy(coor1)
bcoo=copy.deepcopy(coor1)
for i in range(len(coor1)):
    cooa[i][1]=np.round(cooa[i][1]+relat[0][0],4)
    cooa[i][2]=np.round(cooa[i][2]+relat[0][1],4)
for i in range(len(cooa)):
    coormax.append(cooa[i])
for i in range(len(coor1)):
    bcoo[i][1]=np.round(bcoo[i][1]-relat[0][0],4)
    bcoo[i][2]=np.round(bcoo[i][2]-relat[0][1],4)
for i in range(len(bcoo)):
    coormin.append(bcoo[i])

com.close()
pltcx=[]
pltcy=[]
pltcxmono=[]
pltcymono=[]
pltcxy=[]
pltcyy=[]
for i in range(len(coormax)):
    pltcx.append(coormax[i][1])
    pltcy.append(coormax[i][2])
for i in range(len(coor1)):
    pltcxmono.append(coor1[i][1])
    pltcymono.append(coor1[i][2])

coorymax=copy.deepcopy(coor1)
coorymin=copy.deepcopy(coor1)
cooay=copy.deepcopy(coor1)
bcooy=copy.deepcopy(coor1)
for i in range(len(coor1)):
    cooay[i][1]=np.round(cooay[i][1]+relat[1][0],4)
    cooay[i][2]=np.round(cooay[i][2]+relat[1][1],4)
for i in range(len(cooay)):
    coorymax.append(cooay[i])
for i in range(len(coor1)):
    bcooy[i][1]=np.round(bcooy[i][1]-relat[1][0],4)
    bcooy[i][2]=np.round(bcooy[i][2]-relat[1][1],4)
for i in range(len(bcooy)):
    coorymin.append(bcooy[i])    
    
for i in range(len(coorymax)):
    pltcxy.append(coorymax[i][1])
    pltcyy.append(coorymax[i][2])
    

    
coorxy=[]
for i in range(3):
    for k in range(3):
        cooraxy=copy.deepcopy(coor1)
        for j in range(len(coor1)):
            cooraxy[j][1]=np.round(cooraxy[j][1]+(i-1)*relat[1][0],4)
            cooraxy[j][2]=np.round(cooraxy[j][2]+(i-1)*relat[1][1],4)
            cooraxy[j][1]=np.round(cooraxy[j][1]+(k-1)*relat[0][0],4)
            cooraxy[j][2]=np.round(cooraxy[j][2]+(k-1)*relat[0][1],4)
            coorxy.append(cooraxy[j])
nearestnb=[]    
indexvec=0
while True:
    nearestnb=[]    
    for i in range(len(coorxy)):
        if ((lib_a[indexvec][1]-coorxy[i][1])**2+(lib_a[indexvec][2]-coorxy[i][2])**2)**0.5<rcc+0.3 and ((lib_a[indexvec][1]-coorxy[i][1])**2+(lib_a[indexvec][2]-coorxy[i][2])**2)**0.5>0.3:
            nearestnb.append(coorxy[i])
        else:
            pass
    if len(nearestnb)<2:
        indexvec=indexvec+1
    else:
        break
priveca=[]
for i in range(len(nearestnb)):
    priveca.append(np.array([lib_a[indexvec][1]-nearestnb[i][1],lib_a[indexvec][2]-nearestnb[i][2],0]))
if len(priveca)<3:
    priveca.append(-(priveca[0]+priveca[1])*np.linalg.norm(priveca[0])/np.linalg.norm(priveca[0]+priveca[1]))
privecb=[]
for i in range(len(priveca)):
    privecb.append(priveca[i]*-1)

for i in coor1:
    i[1]=str(i[1])
    i[2]=str(i[2])
    i[3]=str(i[3])
for i in range(len(coor1)):
    w='   '.join(coor1[i])
    commax.write(w+'\n')    
    
commax.close()

for i in range(len(coor1)):
    coor1[i][1]=np.round(float(coor1[i][1]),4)
    coor1[i][2]=np.round(float(coor1[i][2]),4)
    coor1[i][3]=np.round(float(coor1[i][3]),4)




if denrec==[]:
    densityup=[]
    densitydw=[]    
    for i in range(len(coor1)):
        densityup.append(0)
    for i in range(len(coor1)):
        densitydw.append(0)
    #for i in range(len(eigek)):
    #    for z in range(len(eigek[i])):
    #        if eigek[i][z]>fermi:
    #            pass
    #        else:
    #            for j in range(len(coor1)):
    #                density[z]=np.round(density[z]+(eigveck[i][j][z]*eigveck[i][j][z].conjugate()).real,4)
    #for i in range(len(density)):
    #    density[i]=density[i]/len(eigek)        
    if guess=='1':
        coundw=0
        counup=0
        randlist=[]
        while True:
            if len(randlist)>=s:
                break
            else:
                tempran=random.randint(0,int(len(coor1)/2))
                if tempran in randlist:
                    pass
                else:
                    randlist.append(tempran)
                    densityup[tempran]=1
                    densitydw[tempran]=0
            
        if UC!=0:
            for i in range(len(coor1)):   
                if coundw<np.round((len(coor1)-s+numberN-numberB)/2,1) and counup<np.round((len(coor1)+s+numberN-numberB)/2,1):       
                    if i in randlist:
                        counup+=1
                    else:
                        if random.choice([True,False]):
                            densityup[i]=1
                            densitydw[i]=0
                            counup+=1
                        else:
                            densityup[i]=0
                            densitydw[i]=1
                            coundw+=1
                else:
                    break
            for i in range(len(coor1)-counup-coundw):
                if i+counup+coundw in randlist:
                    pass
                else:
                    if coundw!=(len(coor1)-s+numberN-numberB)/2 and counup==(len(coor1)+s+numberN-numberB)/2:
                        densityup[i+counup+coundw]=0
                        densitydw[i+counup+coundw]=1
                    else:
                        densityup[i+counup+coundw]=1
                        densitydw[i+counup+coundw]=0
        else:
            for i in range(len(coor1)):
                densityup[i]=0.5
                densitydw[i]=0.5
    elif guess=='4':
        for i in range(len(coor1)):
            if coor1[i][0]=='C1':
                densityup[i]=1
                densitydw[i]=0
            elif coor1[i][0]=='C-1':
                densityup[i]=0
                densitydw[i]=1
            elif coor1[i][0]=='C':
                densityup[i]=0.5
                densitydw[i]=0.5
            else:
                pass
    else:
        print('error')
    nup=0
    ndw=0
    for i in range(len(densityup)):
        nup=nup+densityup[i]
        ndw=ndw+densitydw[i]
else:
    index=0
    for i in reversed(range(len(denrec))):
        if len(denrec[i].split())<3:
            pass
        elif denrec[i].split()[0]!='densitydown:':
            pass
        else:
            densitydw=denrec[i].split()
            index=i
            break
    densitydw.pop(0)
    for i in range(len(densitydw)):
       densitydw[i]=float(densitydw[i])
    densityup=denrec[index-1].split()
    densityup.pop(0)
    for i in range(len(densityup)):
       densityup[i]=float(densityup[i])
    nup=0
    ndw=0
    for i in range(len(densityup)):
        nup=nup+densityup[i]
        ndw=ndw+densitydw[i] 









if __name__=='__main__':
    startime=datetime.datetime.now()
    print('start time : {}'.format(startime))
    if dk[0]!=0 and dk[1]==0:
        plt.figure(figsize=strusize)
        plt.title('{} 2unit structure_s={}_E={}'.format(name,s,E),fontsize=25)
        for i in range(len(coormax)):
            if coormax[i][0]=='C':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='black',markersize=sitems)
            elif coormax[i][0]=='N':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='red',markersize=sitems)
            elif coormax[i][0]=='B':
                plt.plot(coormax[i][1],coormax[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coormax)):
            for z in range(len(coormax)):
                if ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5>0.5 and ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5<1.7:
                    plt.plot([coormax[i][1],coormax[z][1]],[coormax[i][2],coormax[z][2]],color='black')
                else:
                    pass
        plt.axis('off') 
        plt.show()
    elif dk[0]==0 and dk[1]==0:
        plt.figure(figsize=strusize)
        plt.title('{} molecule structure_s={}_E={}'.format(name,s,E),fontsize=25)
        for i in range(len(coor1)):
            if re.match('C',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='black',markersize=sitems)
            elif re.match('N',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='red',markersize=sitems)
            elif re.match('B',coor1[i][0]):
                plt.plot(coor1[i][1],coor1[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coor1)):
            for z in range(len(coor1)):
                if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5>0.5 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5<1.7:
                    plt.plot([coor1[i][1],coor1[z][1]],[coor1[i][2],coor1[z][2]],color='black')
                else:
                    pass
        plt.axis('off') 
        plt.show()
    elif dk[1]!=0 and dk[0]==0:
        plt.figure(figsize=strusize)
        plt.title('{} 2unit_y structure_s={}_E={}'.format(name,s,E),fontsize=25)
        for i in range(len(coorymax)):
            if coorymax[i][0]=='C':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='black',markersize=sitems)
            elif coorymax[i][0]=='N':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='red',markersize=sitems)
            elif coorymax[i][0]=='B':
                plt.plot(coorymax[i][1],coorymax[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coorymax)):
            for z in range(len(coorymax)):
                if ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5>0.5 and ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5<1.7:
                    plt.plot([coorymax[i][1],coorymax[z][1]],[coorymax[i][2],coorymax[z][2]],color='black')
                else:
                    pass
        plt.axis('off') 
        plt.show()    
        
    else:
        plt.figure(figsize=strusize)
        plt.title('{} unit_xy structure_s={}_E={}'.format(name,s,E),fontsize=25)
        for i in range(len(coorxy)):
            if coorxy[i][0]=='C':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='black',markersize=sitems)
            elif coorxy[i][0]=='N':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='red',markersize=sitems)
            elif coorxy[i][0]=='B':
                plt.plot(coorxy[i][1],coorxy[i][2],'o',color='green',markersize=sitems)
            else:
                pass
        for i in range(len(coorxy)):
            for z in range(len(coorxy)):
                if ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5>0.5 and ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5<1.7:
                    plt.plot([coorxy[i][1],coorxy[z][1]],[coorxy[i][2],coorxy[z][2]],color='black')
                else:
                    pass
        plt.axis('off') 
        plt.show()  
     
    
    
    if denrec==[]:
        print('load spinupdensity={}'.format(densityup))
        print('load spindowndensity={}'.format(densitydw))
        print('guess number of spinup e- : {}'.format(np.round((nup),2)))
        print('guess number of spindw e- : {}'.format(np.round((ndw),2)))
    else:
        print('load spinupdensity={}'.format(densityup))
        print('load spindowndensity={}'.format(densitydw))
        print('guess number of spinup e- : {}'.format(np.round((nup),2)))
        print('guess number of spindw e- : {}'.format(np.round((ndw),2)))
    count=1
    realdk=dk
    while True:
        if count<50 and denrec==[] and UC!=0:
            dk=realdk*0.05
        else:
            dk=realdk
            
        #manager=mp.Manager()

        kspacex1=np.linspace(0,bzlat[0][0],int(dk[0]))
        kspacex2=np.linspace(0,bzlat[0][1],int(dk[0]))
        kspacey1=np.linspace(0,bzlat[1][0],int(dk[1]))
        kspacey2=np.linspace(0,bzlat[1][1],int(dk[1]))
        kspacex=[]
        kspacey=[]
        for i in range(len(kspacex1)):
            kspacex.append(np.array([kspacex1[i],kspacex2[i],0]))
        for i in range(len(kspacey1)):
            kspacey.append(np.array([kspacey1[i],kspacey2[i],0]))
        if dk[0]==0:
            kspacex=np.array([[0,0,0]])
        if dk[1]==0:
            kspacey=np.array([[0,0,0]])
            
            
        windup=mp.Manager().dict()
        winddw=mp.Manager().dict()
        eigveckdwtemp=mp.Manager().dict()
        eigveckuptemp=mp.Manager().dict()
        eigekuptemp=mp.Manager().dict()
        eigekdwtemp=mp.Manager().dict()
        
        
        poollist=[]
        for ky in range(np.shape(kspacey)[0]):
            for k in range(np.shape(kspacex)[0]):
                poollist.append((k,ky,eigekuptemp,eigekdwtemp,eigveckuptemp,eigveckdwtemp,windup,winddw,kspacex,kspacey))
            
        
        mypool = mp.Pool(int(mp.cpu_count()/2))

        mypool.starmap(multcore,poollist)
        
        mypool.close()
        mypool.join()
        
        eigveckdw=[]
        eigveckup=[]
        eigekup=[]
        eigekdw=[]
        
        for ky in range(np.shape(kspacey)[0]):
            eigveckdw.append([])
            eigveckup.append([])
            eigekup.append([])
            eigekdw.append([])
            for k in range(np.shape(kspacex)[0]):
                eigekup[ky].append([])
                eigveckup[ky].append([])
                eigekdw[ky].append([])
                eigveckdw[ky].append([])
        
        for ky in range(np.shape(kspacey)[0]):
            for k in range(np.shape(kspacex)[0]):
                eigveckdw[ky][k]=eigveckdwtemp[ky,k]
                eigveckup[ky][k]=eigveckuptemp[ky,k]
                eigekdw[ky][k]=eigekdwtemp[ky,k]
                eigekup[ky][k]=eigekuptemp[ky,k]
        
        
        fermilistup=[]
        fermilistdw=[]
        newdensityup=[]
        newdensitydw=[]
        for z in range(len(kspacey)):
            for i in range(len(kspacex)):
                for j in range(len(eigekup[z][i])):
                    fermilistup.append(eigekup[z][i][j])
        fermilistup.sort()
        fermiup=(fermilistup[int(((len(coor1)+s+numberN-numberB)*len(kspacey)*len(kspacex))/2)]+fermilistup[int(((len(coor1)+s+numberN-numberB)*len(kspacey)*len(kspacex))/2)-1])/2
        for z in range(len(kspacey)):
            for i in range(len(kspacex)):
                for j in range(len(eigekdw[z][i])):
                    fermilistdw.append(eigekdw[z][i][j])
        fermilistdw.sort()
        fermidw=(fermilistdw[int(((len(coor1)-s+numberN-numberB)*len(kspacey)*len(kspacex))/2)]+fermilistdw[int(((len(coor1)-s+numberN-numberB)*len(kspacey)*len(kspacex))/2)-1])/2
        for i in range(len(coor1)):
            newdensityup.append(0)
            newdensitydw.append(0)  
        for ky in range(len(kspacey)):
            for i in range(len(kspacex)):
                for z in range(len(eigekup[ky][i])):
                    if eigekup[ky][i][z]>fermiup:
                        pass
                    else:
                        for j in range(len(eigveckup[ky][i])):                    
                            newdensityup[j]=newdensityup[j]+(eigveckup[ky][i][j][z]*eigveckup[ky][i][j][z].conjugate()).real
        for ky in range(len(kspacey)):
            for i in range(len(kspacex)):
                for z in range(len(eigekdw[ky][i])):
                    if eigekdw[ky][i][z]>fermidw:
                        pass
                    else:
                        for j in range(len(coor1)):
                            newdensitydw[j]=newdensitydw[j]+(eigveckdw[ky][i][j][z]*eigveckdw[ky][i][j][z].conjugate()).real   
        for i in range(len(newdensityup)):
            newdensityup[i]=newdensityup[i]/(len(kspacex)*len(kspacey)) 
            newdensitydw[i]=newdensitydw[i]/(len(kspacex)*len(kspacey)) 
        diff=0
        for i in range(len(densityup)):
            diff=diff+abs(densityup[i]-newdensityup[i])+abs(densitydw[i]-newdensitydw[i])
        if diff<selcons:
            for i in range(len(densitydw)):
                densitydw[i]=newmixrate*newdensitydw[i]+(1-newmixrate)*densitydw[i]
                densityup[i]=newmixrate*newdensityup[i]+(1-newmixrate)*densityup[i]
            print('step{} : diff = {} < {}\njob done  （<ゝω・）☆'.format(count,np.round(diff,9),selcons))
            with open('{}/{}density_'.format(name,name)+folder+'.txt','a+') as densityrecord:
                densityrecord.write('densityup: ')
                for i in range(len(newdensityup)):
                    newdensityup[i]=np.round(newdensityup[i],8)
                    newdensityup[i]=str(newdensityup[i])
                w='   '.join(newdensityup)
                densityrecord.write(w+'\n')   
                densityrecord.write('densitydown: ')
                for i in range(len(newdensitydw)):
                    newdensitydw[i]=np.round(newdensitydw[i],8)
                    newdensitydw[i]=str(newdensitydw[i])
                w='   '.join(newdensitydw)
                densityrecord.write(w+'\n')
                densityrecord.write('\n')
            break
        else:
            for i in range(len(densitydw)):
                densitydw[i]=newmixrate*newdensitydw[i]+(1-newmixrate)*densitydw[i]
                densityup[i]=newmixrate*newdensityup[i]+(1-newmixrate)*densityup[i]
            print('step{} : diff = {} > {}'.format(count,np.round(diff,6),selcons))
            count=count+1
            with open('{}/{}density_'.format(name,name)+folder+'.txt','a+') as densityrecord:
                densityrecord.write('densityup: ')
                for i in range(len(newdensityup)):
                    newdensityup[i]=np.round(newdensityup[i],8)
                    newdensityup[i]=str(newdensityup[i])
                w='   '.join(newdensityup)
                densityrecord.write(w+'\n')   
                densityrecord.write('densitydown: ')
                for i in range(len(newdensitydw)):
                    newdensitydw[i]=np.round(newdensitydw[i],8)
                    newdensitydw[i]=str(newdensitydw[i])
                w='   '.join(newdensitydw)
                densityrecord.write(w+'\n')
                densityrecord.write('\n')

            
            
            
    dk=realdk
    for z in range(len(kspacey)):
        for i in range(len(kspacex)):
            eigveckup[z][i]=eigveckup[z][i][:,eigekup[z][i].argsort()]
            eigveckdw[z][i]=eigveckdw[z][i][:,eigekdw[z][i].argsort()]
            eigekup[z][i]=np.round(np.sort(eigekup[z][i]),8)    
            eigekdw[z][i]=np.round(np.sort(eigekdw[z][i]),8)
        
    nup=0
    ndw=0
    for i in range(len(densityup)):
        nup=nup+densityup[i]
        ndw=ndw+densitydw[i]
    print('final number of spinup e-:{}'.format(np.round((nup),2)))
    print('final number of spindw e-:{}'.format(np.round((ndw),2)))
    energy=0
    for i in range(len(fermilistup)):
        if fermilistup[i]<fermiup:
            energy=energy+fermilistup[i]
        else:
            pass
    for i in range(len(fermilistdw)):
        if fermilistdw[i]<fermidw:
            energy=energy+fermilistdw[i]
        else:
            pass
    energy=np.round((energy/(len(kspacex)*len(kspacey))).real,4)
    for i in range(len(coor1)):
        energy=energy+0.5*UC*(densityup[i]*densitydw[i])
    
    if berrycal=='on' :
        berryca(dk)
        
    else:
        pass
    
    
    if windnum=='on':
        if dk[0]!=0 and dk[1]!=0:
            cauwind2D()
        else:
            cauwind1D()
    else:
        pass
    resu=open('{}/'.format(name)+folder+'/{}_energy_E={}_s={}.txt'.format(name,E,s),'w+')
    print('{}_E={}_s={} Energy = {}'.format(name,str(E),str(s),energy))
    resu.write('Energy = {}'.format(energy))
    resu.close()
    
    totden=[]
    tot=[]
    if dk[1]!=0 and dk[0]!=0:
        band2D()
    else:
        bandorb1D()
    for i in range(len(densitydw)):       
        totden.append(0)
        tot.append(0)
        totden[i]=densitydw[i]-densityup[i]
        tot[i]=densitydw[i]+densityup[i]    
    abtotden=[]
    for i in range(len(densitydw)):       
        abtotden.append(0)
        abtotden[i]=abs(totden[i])
    sumtotden=0
    for i in range(len(densitydw)):
        sumtotden=sumtotden+abs(totden[i])
    if UC==0:
        pass
    else:
        plt.figure(figsize=strusize)
        plt.title('{} spindensity_E={}_s={}_(n{}){}={}'.format(name,str(E),str(s),'$_{spin}$','$_{total}$',np.round(sumtotden,4)),fontsize=25)
        if dk[0]!=0 and dk[1]==0 :
            for i in range(len(coormax)):
                if coormax[i][0]=='C':
                    plt.plot(coormax[i][1],coormax[i][2],'o',color='black',markersize=sitems)
                elif coormax[i][0]=='N':
                    plt.plot(coormax[i][1],coormax[i][2],'o',color='red',markersize=sitems)
                elif coormax[i][0]=='B':
                    plt.plot(coormax[i][1],coormax[i][2],'o',color='green',markersize=sitems)
                else:
                    pass
        elif dk[1]!=0 and dk[0]==0:
            for i in range(len(coorymax)):
                if coorymax[i][0]=='C':
                    plt.plot(coorymax[i][1],coorymax[i][2],'o',color='black',markersize=sitems)
                elif coorymax[i][0]=='N':
                    plt.plot(coorymax[i][1],coorymax[i][2],'o',color='red',markersize=sitems)
                elif coorymax[i][0]=='B':
                    plt.plot(coorymax[i][1],coorymax[i][2],'o',color='green',markersize=sitems)
                else:
                    pass
        elif dk[1]!=0 and dk[0]!=0:
            for i in range(len(coorxy)):
                if coorxy[i][0]=='C':
                    plt.plot(coorxy[i][1],coorxy[i][2],'o',color='black',markersize=sitems)
                elif coorxy[i][0]=='N':
                    plt.plot(coorxy[i][1],coorxy[i][2],'o',color='red',markersize=sitems)
                elif coorxy[i][0]=='B':
                    plt.plot(coorxy[i][1],coorxy[i][2],'o',color='green',markersize=sitems)
                else:
                    pass
            for i in range(len(coorxy)):
                for z in range(len(coorxy)):
                    if ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5>0.5 and ((coorxy[i][1]-coorxy[z][1])**2+(coorxy[i][2]-coorxy[z][2])**2)**0.5<1.7:
                        plt.plot([coorxy[i][1],coorxy[z][1]],[coorxy[i][2],coorxy[z][2]],color='black')
                    else:
                        pass
        else:
            for i in range(len(coor1)):
                if re.match('C',coor1[i][0]):
                    plt.plot(coor1[i][1],coor1[i][2],'o',color='black',markersize=sitems)
                elif re.match('N',coor1[i][0]):
                    plt.plot(coor1[i][1],coor1[i][2],'o',color='red',markersize=sitems)
                elif re.match('B',coor1[i][0]):
                    plt.plot(coor1[i][1],coor1[i][2],'o',color='green',markersize=sitems)
                else:
                    pass
                
        if dk[0]!=0 and dk[1]==0 :
            for i in range(len(coormax)):
                for z in range(len(coormax)):
                    if ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5>0.5 and ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5<1.7:
                        plt.plot([coormax[i][1],coormax[z][1]],[coormax[i][2],coormax[z][2]],color='black')
                    else:
                        pass
            
            for i in range(len(pltcx)):
                if i>len(coor1)-1:
                    if totden[i-len(coor1)]<0:
                            plt.scatter(pltcx[i],pltcy[i],s=3000,c='red',alpha=(abs(totden[i-len(coor1)])/max(abtotden)))
                    else:
                            plt.scatter(pltcx[i],pltcy[i],s=3000,c='blue',alpha=(abs(totden[i-len(coor1)])/max(abtotden)))
                            
                            
                else:
                    if totden[i]<0:
                        plt.scatter(pltcx[i],pltcy[i],s=3000,c='red',alpha=(abs(totden[i])/max(abtotden)))
                    else:
                        plt.scatter(pltcx[i],pltcy[i],s=3000,c='blue',alpha=(abs(totden[i])/max(abtotden)))
        elif dk[1]!=0 and dk[0]==0 :
            for i in range(len(coorymax)):
                for z in range(len(coorymax)):
                    if ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5>0.5 and ((coorymax[i][1]-coorymax[z][1])**2+(coorymax[i][2]-coorymax[z][2])**2)**0.5<1.7:
                        plt.plot([coorymax[i][1],coorymax[z][1]],[coorymax[i][2],coorymax[z][2]],color='black')
                    else:
                        pass
            
            for i in range(len(pltcxy)):
                if i>len(coor1)-1:
                    if totden[i-len(coor1)]<0:
                            plt.scatter(pltcxy[i],pltcyy[i],s=3000,c='red',alpha=(abs(totden[i-len(coor1)])/max(abtotden)))
                    else:
                            plt.scatter(pltcxy[i],pltcyy[i],s=3000,c='blue',alpha=(abs(totden[i-len(coor1)])/max(abtotden)))
                            
                            
                else:
                    if totden[i]<0:
                        plt.scatter(pltcxy[i],pltcyy[i],s=3000,c='red',alpha=(abs(totden[i])/max(abtotden)))
                    else:
                        plt.scatter(pltcxy[i],pltcyy[i],s=3000,c='blue',alpha=(abs(totden[i])/max(abtotden)))
        elif dk[1]!=0 and dk[0]!=0 :
            
            for i in range(len(coorxy)):                        
                if totden[i]<0:
                    plt.scatter(coorxy[i][1],coorxy[i][2],s=3000,c='red',alpha=(abs(totden[i%len(coor1)])/max(abtotden)))
                else:
                    plt.scatter(pltcxy[i][1],pltcyy[i][2],s=3000,c='blue',alpha=(abs(totden[i%len(coor1)])/max(abtotden)))
        else:
            for i in range(len(coor1)):
                for z in range(len(coor1)):
                    if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5>0.5 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5<1.7:
                        plt.plot([coor1[i][1],coor1[z][1]],[coor1[i][2],coor1[z][2]],color='black')
                    else:
                        pass
            
            for i in range(len(pltcxmono)):
                if totden[i]<0:
                    plt.scatter(pltcxmono[i],pltcymono[i],s=3000,c='red',alpha=(abs(totden[i])/max(abtotden)))
                else:
                    plt.scatter(pltcxmono[i],pltcymono[i],s=3000,c='blue',alpha=(abs(totden[i])/max(abtotden)))        
        plt.axis('off')    
        plt.savefig('{}/'.format(name)+folder+'/{}spindensity_E={}_s={}.jpg'.format(name,str(E),str(s)))
        plt.show()
    #plt.figure(figsize=strusize)
    #plt.title('{} totaldensity_E={}_s={}'.format(name,str(E),str(s)),fontsize=25)
    #if dk[0]!=0:
    #    for i in range(len(coormax)):
    #        if coormax[i][0]=='C':
    #            plt.plot(coormax[i][1],coormax[i][2],'o',color='black',markersize=sitems)
    #        elif coormax[i][0]=='N':
    #            plt.plot(coormax[i][1],coormax[i][2],'o',color='red',markersize=sitems)
    #        elif coormax[i][0]=='B':
    #            plt.plot(coormax[i][1],coormax[i][2],'o',color='green',markersize=sitems)
    #        else:
    #            pass
    #else:
    #    for i in range(len(coor1)):
    #        if re.match('C',coor1[i][0]):
    #            plt.plot(coor1[i][1],coor1[i][2],'o',color='black',markersize=sitems)
    #        elif re.match('N',coor1[i][0]):
    #            plt.plot(coor1[i][1],coor1[i][2],'o',color='red',markersize=sitems)
    #        elif re.match('B',coor1[i][0]):
    #            plt.plot(coor1[i][1],coor1[i][2],'o',color='green',markersize=sitems)
    #        else:
    #            pass
    #if dk[0]!=0:
    #    for i in range(len(coormax)):
    #        for z in range(len(coormax)):
    #            if ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5>0.5 and ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5<1.7:
    #                plt.plot([coormax[i][1],coormax[z][1]],[coormax[i][2],coormax[z][2]],color='black')
    #            else:
    #                pass
    #    for i in range(len(pltcx)):
    #        if i>len(coor1)-1:
    #            plt.scatter(pltcx[i],pltcy[i],s=2000,c='grey',alpha=(tot[i-len(coor1)]/max(tot)))                
    #        else:
    #            plt.scatter(pltcx[i],pltcy[i],s=2000,c='grey',alpha=((tot[i]/max(tot))))
    #else:
    #    for i in range(len(coor1)):
    #        for z in range(len(coor1)):
    #            if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5>0.5 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5<1.7:
    #                plt.plot([coor1[i][1],coor1[z][1]],[coor1[i][2],coor1[z][2]],color='black')
    #            else:
    #                pass
    #    for i in range(len(pltcxmono)):
    #        plt.scatter(pltcxmono[i],pltcymono[i],s=2000,c='grey',alpha=((tot[i]/max(tot))))    
        
    #plt.axis('off')   
    #plt.savefig('{}/E={}_s={}_U={}/{}totaldensity_E={}_s={}.jpg'.format(name,E,s,UC,name,str(E),str(s))) 
    #plt.show()
    #def plotmo(typ,inde,loc,spin):
    #    plt.figure(figsize=strusize)
    #    if inde!=0:
    #        if inde<0:
    #            plt.title('{} {}{}(k={})_E={}_s={}_{}'.format(name,typ,inde,loc,str(E),str(s),spin),fontsize=25)
    #        else:
    #            plt.title('{} {}+{}(k={})_E={}_s={}_{}'.format(name,typ,inde,loc,str(E),str(s),spin),fontsize=25)
    #    else:
    #        plt.title('{} {}(k={})_E={}_s={}_{}'.format(name,typ,loc,str(E),str(s),spin),fontsize=25)
    #    if dk[0]!=0:
    #        for i in range(len(coormax)):
    #            if coormax[i][0]=='C':
    #                plt.plot(coormax[i][1],coormax[i][2],'o',color='black',markersize=sitems)
    #            elif coormax[i][0]=='N':
    #                plt.plot(coormax[i][1],coormax[i][2],'o',color='red',markersize=sitems)
    #            elif coormax[i][0]=='B':
    #                plt.plot(coormax[i][1],coormax[i][2],'o',color='green',markersize=sitems)
    #            else:
    #                pass
    #        for i in range(len(coormax)):
    #            for z in range(len(coormax)):
    #                if ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5>0.5 and ((coormax[i][1]-coormax[z][1])**2+(coormax[i][2]-coormax[z][2])**2)**0.5<1.7:
    #                    plt.plot([coormax[i][1],coormax[z][1]],[coormax[i][2],coormax[z][2]],color='black')
    #                else:
    #                    pass
    #    else:
    #        for i in range(len(coor1)):
    #            if re.match('C',coor1[i][0]):
    #                plt.plot(coor1[i][1],coor1[i][2],'o',color='black',markersize=sitems)
    #            elif re.match('N',coor1[i][0]):
    #                plt.plot(coor1[i][1],coor1[i][2],'o',color='red',markersize=sitems)
    #            elif re.match('B',coor1[i][0]):
    #                plt.plot(coor1[i][1],coor1[i][2],'o',color='green',markersize=sitems)
    #            else:
    #                pass
    #        for i in range(len(coor1)):
    #            for z in range(len(coor1)):
    #                if ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5>0.5 and ((coor1[i][1]-coor1[z][1])**2+(coor1[i][2]-coor1[z][2])**2)**0.5<1.7:
    #                    plt.plot([coor1[i][1],coor1[z][1]],[coor1[i][2],coor1[z][2]],color='black')
    #                else:
    #                    pass
    #    indexo=0
    #    if loc=='pi':
    #        if spin=='up':
    #
    #            for i in range(len(eigekup[int((dk-1)/2)])):
    #                if eigekup[int((dk-1)/2)][i]<fermiup:
    #                    pass
    #                else:
    #                    if typ=='homo':
    #                        indexo=i-1
    #                        break
    #                    else:
    #                        indexo=i
    #                        break
    #        else:
    #
    #            for i in range(len(eigekdw[int((dk-1)/2)])):
    #                if eigekdw[int((dk-1)/2)][i]<fermidw:
    #                    pass
    #                else:
    #                    if typ=='homo':
    #                        indexo=i-1
    #                
    #                        break
    #                    else:
    #                        indexo=i
    #                        break
    #    else:
    #        if spin=='up':
    #
    #            for i in range(len(eigekup[0])):
    #                if eigekup[0][i]<fermiup:
    #                    pass
    #                else:
    #                    if typ=='homo':
    #                        indexo=i-1
    #                        break
    #                    else:
    #                        indexo=i
    #                        break
    #        else:
    #
    #            for i in range(len(eigekdw[0])):
    #                if eigekdw[0][i]<fermidw:
    #                    pass
    #                else:
    #                    if typ=='homo':
    #                        indexo=i-1
    #                        break
    #                    else:
    #                        indexo=i
    #                        break
    #    orden=[]
    #    indexo=indexo+inde
    #    if spin=='up':
    #        if loc=='pi':
    #            for i in range(len(eigveckup[int((dk-1)/2)])):
    #                orden.append(0)
    #                orden[i]+=eigveckup[int((dk-1)/2)][i,indexo]
    #                orden[i]=abs(orden[i])
    #        else:
    #            for i in range(len(eigveckup[0])):
    #                orden.append(0)
    #                orden[i]+=eigveckup[0][i,indexo]
    #                orden[i]=abs(orden[i])
    #    else:
    #        if loc=='pi':
    #            for i in range(len(eigveckdw[int((dk-1)/2)])):
    #                orden.append(0)
    #                orden[i]+=eigveckdw[int((dk-1)/2)][i,indexo]
    #                orden[i]=abs(orden[i])
    #        else:
    #            for i in range(len(eigveckdw[0])):
    #                orden.append(0)
    #                orden[i]+=eigveckdw[0][i,indexo]
    #                orden[i]=abs(orden[i])
    #    if dk[0]!=0:
    #        for i in range(len(pltcx)):
    #            if i>len(coor1)-1:
    #                if loc=='pi':
    #                    if spin=='up':
    #                        if np.round(cmath.phase(cmath.exp(1j*cmath.pi)*eigveckup[int((dk-1)/2)][i-len(coor1),indexo]),4)>=0 and np.round(cmath.exp(1j*cmath.pi)*cmath.phase(eigveckup[int((dk-1)/2)][i-len(coor1),indexo]),4)<=cmath.pi:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='blue',alpha=0.3) 
    #                    else:
    #                        if np.round(cmath.phase(cmath.exp(1j*cmath.pi)*eigveckdw[int((dk-1)/2)][i-len(coor1),indexo]),4)>=0 and np.round(cmath.exp(1j*cmath.pi)*cmath.phase(eigveckdw[int((dk-1)/2)][i-len(coor1),indexo]),4)<=cmath.pi:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='blue',alpha=0.3)
    #                else:
    #                    if spin=='up':
    #                        if eigveckup[0][i-len(coor1),indexo].real>0:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='blue',alpha=0.3) 
    #                    else:
    #                        if eigveckdw[0][i-len(coor1),indexo].real>0:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i-len(coor1)]/max(orden)*4000,c='blue',alpha=0.3)
    #            else:
    #                if spin=='up':
    #                    if loc=='pi':
    #                        if np.round(cmath.phase(eigveckup[int((dk-1)/2)][i,indexo]),4)>=0 and np.round(cmath.phase(eigveckup[int((dk-1)/2)][i,indexo]),4)<=cmath.pi:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #                    else:
    #                        if eigveckup[0][i,indexo].real>0:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #                else:
    #                    if loc=='pi':
    #                        
    #                        if np.round(cmath.phase(eigveckdw[int((dk-1)/2)][i,indexo]),4)>=0 and np.round(cmath.phase(eigveckdw[int((dk-1)/2)][i,indexo]),4)<=cmath.pi:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #                    else:
    #                        if eigveckdw[0][i,indexo].real>0:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                        else:
    #                            plt.scatter(pltcx[i],pltcy[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #    else:
    #        for i in range(len(pltcxmono)):
    #            if spin=='up':
    #                if np.round(eigveckup[0][i,indexo].real,4)>0:
    #                    plt.scatter(pltcxmono[i],pltcymono[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                else:
    #                    plt.scatter(pltcxmono[i],pltcymono[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #            else:
    #                if np.round(eigveckdw[0][i,indexo].real,4)>0 :
    #                    plt.scatter(pltcxmono[i],pltcymono[i],s=orden[i]/max(orden)*4000,c='red',alpha=0.3)    
    #                else:
    #                    plt.scatter(pltcxmono[i],pltcymono[i],s=orden[i]/max(orden)*4000,c='blue',alpha=0.3) 
    #    
    #                
    # 
    #    plt.axis('off')   
    #    if inde!=0:
    #        if inde<0:
    #            plt.savefig('{}/E={}_s={}_U={}/{}{}{}(k={})_E={}_s={}_{}.jpg'.format(name,E,s,UC,name,typ,inde,loc,str(E),str(s),spin)) 
    #        else:
    #            plt.savefig('{}/E={}_s={}_U={}/{}{}+{}(k={})_E={}_s={}_{}.jpg'.format(name,E,s,UC,name,typ,inde,loc,str(E),str(s),spin)) 
    #    else:
    #        plt.savefig('{}/E={}_s={}_U={}/{}{}_E={}_s={}_{}.jpg'.format(name,E,s,UC,name,typ,str(E),str(s),spin)) 
    #    plt.show()
    
    endtime=datetime.datetime.now()
    print('end time : {}'.format(endtime))
    print('duration : {}'.format(endtime-startime))
    #for i in range(5):
    #    plotmo('lumo',i,'pi','up')
    #    plotmo('lumo',i,'0','up')
    #for i in range(len(coor1)):
    #    for k in range(len(coor1)):
    #        if coor1[k][0]=='C*':
    #            if abs(eigveckup[0][k,i])>0.25:
    #                print(k)                
    #            else:
    #                pass
    #        else:
    #            pass
    #plt.figure()
    #plt.ylabel('energy (eV)')
    #plt.xlabel('index of orbital')
    #for i in range(len(eigekup[0])):    
    #    if abs(eigveckup[0][6,i])>0.25 or abs(eigveckup[0][115,i])>0.25 or abs(eigveckup[0][198,i])>0.25 or abs(eigveckup[0][6,i])>0.25 :
    #        plt.plot(i,eigekup[0][i],'o',color='black',markersize=0.7)
    #    else:
    #        plt.plot(i,eigekup[0][i],'o',color='black',markersize=0.7)
    #for i in range(len(coormax)):
    #    for j in range(1,4):
    #        coormax[i][j]=str(coormax[i][j])
    #for i in range(len(coormax)):
    #    w='   '.join(coormax[i])
    #    commax1.write(w+'\n') 
    #commax1.close()
    #densityrecord.close()
    #K=2*cmath.pi/(3*1.40328)
    #Kup,Kdw=hmatrix(np.array([K,K/(3**0.5),0]))
    #Ke,Kv=np.linalg.eig(Kup)
    #Kesoc.append(abs(Ke[0]-Ke[1]))
    #Kesocx.append(0.005*soci)
        
    #from sklearn.linear_model import LinearRegression
    #fitmodel=LinearRegression()
    #fitmodel.fit(Kesocx,Kesoc)
    #plt.figure()
    
    #plt.scatter(Kesocx,Kesoc,s=0.1)
#######wannier function########
    #den=[]
    #for i in range(len(coor1)):
    #    den.append(0)
    #for i in range(len(coor1)):
    #    for j in range(int(len(coor1)/2)-1):
    #        den[i]=den[i]+abs(eigveckup[0][0][i,j])**2
    #unita=[]
    #unitb=[]
    #for i in range(30):
    #    unita.append(0)
    #    unitb.append(0)
    #for i in range(30):
    #    for j in range(len(coor1)):
    #        if coor1[j][1]>-0.1+min(pltcx)+4.26960*i and coor1[j][1]<-0.5+min(pltcx)+4.26960*(i+1):
    #            if coor1[j] in lib_b:
    #                unitb[i]=unitb[i]+den[j]
    #            else:
    #                unita[i]=unita[i]+den[j]
    #unitta=[]
    #for i in unita:
    #    unitta.append(i-3.483)
    #unittb=[]
    #for i in unitb:
    #    unittb.append(i-3.483)
    #plt.bar(list(range(1,31)),unitta)
    #plt.bar(list(range(1,31)),unittb,color='red')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            