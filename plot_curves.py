# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 14:47:09 2016
plot results from  c++

@author: YPC
"""

import csv
import numpy as np
import matplotlib.pyplot as plt

a=[]

a=[]
nn = "SRDD"
with open("%s.txt"%nn) as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
         # format    0 size=kappa , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time,
		#           5 CG_pre count 6 CG_pre time
		#           7 LU time, 8 LU delx, 9 QR time, 10 QR delx
		#						11 full QR time 12 full QR delx
		#						13 Jacobi count 14 Jacobi Time
		#						15 SOR1 "-" 16
		#						17 SOR1.5 "-" 18
		#						19 SOR0.5 "-" 20  

a=np.array(a,dtype=np.float)
a[a==0]=np.NaN
#a = a[:-2,:]       

fig = plt.figure()
#plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')
plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')


#plt.loglog(a[:,0], a[:,4],'yx-', label = 'CG')
plt.loglog(a[:,0], a[:,16],'x-', label = 'GS')


plt.loglog(a[:,0], a[:,14],'gx-', label = 'Jacobi')



plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')


plt.loglog(a[:,0], a[:,2],'bx-', label = 'SD')

plt.loglog(a[:,0], a[:,18],'x-', label = 'SOR1.5')
plt.loglog(a[:-(sum(np.isnan(a[:,20]))),0], a[:-(sum(np.isnan(a[:,20]))),20],'x-', label = 'SOR0.5')


# make some fits            
y = np.log(a[:,14]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='g',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
             
y = np.log(a[:,7]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='r',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

y = np.log(a[:,2]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))


plt.xlim(xmax=max(a[:,0]+200))
plt.title(r'Time of solving SRDD, $\kappa=2$',fontsize=15)
plt.xlabel('Size $N$',fontsize=15)
plt.ylabel('Time in s',fontsize=15)
#plt.legend(loc='best', ncol=2,fontsize=11)
plt.gca().legend(ncol=2)#(loc='center left', ncol=2, bbox_to_anchor=(1, 0.5))
fig.savefig('%s.png'%nn,bbox_inches='tight')

#%%

a=[]
with open("../output/same_size.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
     # format    0 kappa , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
     #           5 CG_pre count 6 CG_pre time 
     #           7 LU time, 8 LU delx, 9 QR time, 10 QR delx  11 full QR time 12 full QR delx       
             
a=np.array(a,dtype=np.float)
fig = plt.figure()

####
#plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')
plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')
plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')
plt.loglog(a[:,0], a[:,4],'yx-', label = 'CG')
plt.loglog(a[:,0], a[:,2],'gx-', label = 'SD')
plt.loglog(a[:,0], a[:,6],'bx-', label = 'CG-pre')


# make some fits
#y = np.log(a[:,11]) 
#x= np.log(a[:,0])
#m,b = np.polyfit(x, y, 1)
#y10= np.exp(m*np.log(a[0,0]) +b)     
#y100= np.exp(m*np.log(a[-(sum(np.isnan(a[:,11]))),0] )+b)
#
#plt.loglog([a[0,0], a[-(sum(np.isnan(a[:,11]))),0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
#               markersize=10,label='Lin. fit m=%s'%round(m,2))

y = np.log(a[:,9]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-(sum(np.isnan(a[:,9]))),0] )+b)

plt.loglog([a[0,0], a[-(sum(np.isnan(a[:,9]))),0]],[y10,y100] ,'--', marker='.',color='c',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
               
               
y = np.log(a[:,7]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='r',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))




y = np.log(a[:,4]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='y',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))


y = np.log(a[:,2]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='g',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

y = np.log(a[:,6]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

####


plt.title(r'Time of solving for different $\kappa$ and $N=100$',fontsize=15)
plt.ylabel('Time in s',fontsize=15)
plt.xlabel(r'$\kappa$',fontsize=15)
plt.legend(loc='best', ncol=2)
fig.savefig('Kappa_time.pdf')

#%%
import csv
import numpy as np
import matplotlib.pyplot as plt

for nn in ["k_2_normal","k_2","k_10"]:
    print "%s.txt"%nn
    a=[]
    with open("%s.txt"%nn) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter=",")
        for line in tsvreader:
            a.append( line)
     # format    0 size , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
     #           5 CG_pre count 6 CG_pre time 
     #           7 LU time, 8 LU delx, 9 QR time, 10 QR delx  11 full QR time 12 full QR delx       
    if nn =="k_10":
        kapp = 10
    else:
        kapp=2
    a=np.array(a,dtype=np.float)
    a[a==0]=np.NaN
    #a = a[:-2,:]       
    
    fig = plt.figure()
    plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')
    plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')
    plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')
    plt.loglog(a[:,0], a[:,4],'yx-', label = 'CG')
    plt.loglog(a[:,0], a[:,2],'gx-', label = 'SD')
    plt.loglog(a[:,0], a[:,6],'bx-', label = 'CG-pre')
    
    
    
    y = np.log(a[:-(sum(np.isnan(a[:,11]))),11]) 
    x= np.log(a[:-(sum(np.isnan(a[:,11]))),0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)     
    y100= np.exp(m*np.log(a[-(sum(np.isnan(a[:,11]))),0] )+b)
    
    plt.loglog([a[0,0], a[-(sum(np.isnan(a[:,11]))),0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
    
    y = np.log(a[:-(sum(np.isnan(a[:,9]))),9]) 
    x= np.log(a[:-(sum(np.isnan(a[:,9]))),0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)
    y100= np.exp(m*np.log(a[-(sum(np.isnan(a[:,9]))),0] )+b)
    
    plt.loglog([a[0,0], a[-(sum(np.isnan(a[:,9]))),0]],[y10,y100] ,'--', marker='.',color='c',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
                   
                   
    y = np.log(a[:,7]) 
    x= np.log(a[:,0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)
    y100= np.exp(m*np.log(a[-1,0] )+b)
    
    plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='r',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
    
    
    
    
    y = np.log(a[:,4]) 
    x= np.log(a[:,0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)
    y100= np.exp(m*np.log(a[-1,0] )+b)
    
    plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='y',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
    
    
    y = np.log(a[:,2]) 
    x= np.log(a[:,0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)
    y100= np.exp(m*np.log(a[-1,0] )+b)
    
    plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='g',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
    
    y = np.log(a[:,6]) 
    x= np.log(a[:,0])
    m,b = np.polyfit(x, y, 1)
    y10= np.exp(m*np.log(a[0,0]) +b)
    y100= np.exp(m*np.log(a[-1,0] )+b)
    
    plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=2,
                   markersize=10,label='Lin. fit m=%s'%round(m,2))
    
    
    plt.xlim(xmax=max(a[:,0]+100))
    plt.title(r'Time of solving for different sizes, $\kappa =%s$'%kapp,fontsize=15)
    plt.xlabel('Size N',fontsize=15)
    plt.ylabel('Time in seconds',fontsize=15)
    plt.legend(loc='best', ncol=2,fontsize=11)
    fig.savefig('%s.pdf'%nn)
