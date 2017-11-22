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

with open("Matrixclass/SRDD2.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
 # format    0 size , 1 JacCount, 2 JacTime, 3 SOR1Count, 4 SOR1Time
 #           5 SOR1.5Count, 6 SOR1.5 Time, 7 SOR05Count, 8 SOR05Time
 #           9 SD Count, 10 SD Time, 11 CG Count, 12 CG time, 
 #           13 LU time, 14 LU delx, 15 QR time, 16 QR delx           
a=np.array(a,dtype=np.float)
#a = a[:-2,:]       
#%%
fig = plt.figure()
plt.loglog(a[:,0], a[:,2],'x-', label = 'Jacobi')
plt.loglog(a[:,0], a[:,8],'x-', label = 'SOR 0.5')
plt.loglog(a[:,0], a[:,10],'mx-', label = 'SD')


plt.loglog(a[:,0], a[:,12],'yx-', label = 'CG')
plt.loglog(a[:,0], a[:,13],'kx-', label = 'LU')

plt.loglog(a[:,0], a[:,4],'gx-', label = 'SOR 1')
plt.loglog(a[:,0], a[:,6],'x-', label = 'SOR 1.5')


plt.loglog(a[:,0], a[:,15],'x-', label = 'QR')

y = np.log(a[:,10]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
y = np.log(a[:,12]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='y',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
               

               

y = np.log(a[:,13]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='k',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
y = np.log(a[:,4]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='g',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

plt.title(r'Time of solving for different sizes',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.xlim(xmax=max(a[:,0]))
plt.ylabel('Time in seconds',fontsize=15)
plt.legend(loc='best',ncol=2,fontsize=10)
fig.savefig('Tex/figs/SRDD_time.png', dpi=500)


#%%
fig = plt.figure()
plt.loglog(a[:,0], a[:,1],'x-', label = 'Jacobi')
plt.loglog(a[:,0], a[:,3],'x-', label = 'SOR 1')
plt.loglog(a[:,0], a[:,5],'x-', label = 'SOR 1.5')
plt.loglog(a[:,0], a[:,7],'x-', label = 'SOR 0.5')

plt.loglog(a[:,0], a[:,9],'x-', label = 'SD')
plt.loglog(a[:,0], a[:,11],'x-', label = 'CG')
plt.xlim(xmax=max(a[:,0]))
plt.title(r'Number of iterations for different sizes',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.ylabel('Iterations',fontsize=15)
plt.legend(loc='best',ncol=2,fontsize=10)
fig.savefig('Tex/figs/SRDD_iters.png', dpi=500)
#%%

y = np.log(a[:,16]) 
x= np.log(a[:,0])
    
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

fig3 = plt.figure()
plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='g',linewidth=3,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
plt.loglog(a[:,0], a[:,14],'x-', label = 'LU')
plt.loglog(a[:,0], a[:,16],'x-', label = 'QR')
plt.xlim(xmax=max(a[:,0]))
plt.title('Accuracy for different sizes',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.ylabel('Error',fontsize=15)
plt.legend(loc='best', fontsize=15, ncol=2)
fig3.savefig('Tex/figs/Tex/figs/SRDD_Accuracy.png', dpi=500)

#%%
a=[]

with open("Matrixclass/difkappa.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
 # format    0 kappa , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
 #           5 LU time, 6 LU delx, 7 QR time, 8 QR delx           
a=np.array(a,dtype=np.float)
y = np.log(a[:65,2]) 
x= np.log(a[:65,0])
    
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[65,0] )+b)
fig = plt.figure(2)

plt.loglog(a[:65,0], a[:65,2],'x-', label = 'SD')
plt.loglog(a[:65,0], a[:65,4],'x-', label = 'CG')
plt.loglog(a[:65,0], a[:65,5],'x-', label = 'LU')
plt.loglog(a[:65,0], a[:65,7],'x-', label = 'QR')
plt.loglog([a[0,0], a[65,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=3,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

plt.title(r'Time of solving for different $\kappa$ and $N=50$',fontsize=15)
plt.ylabel('Time in s',fontsize=15)
plt.xlabel(r'$\kappa$',fontsize=15)
plt.legend(loc='best', ncol=2)
fig.savefig('Tex/figs/Kappa_time.png', dpi=500)

y = np.log(a[:65,6]) 
x= np.log(a[:65,0])
    
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[65,0] )+b)

fig3 = plt.figure()
plt.loglog([a[0,0], a[65,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=3,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
plt.loglog(a[:65,0], a[:65,6],'x-', label = 'LU')
plt.loglog(a[:65,0], a[:65,8],'x-', label = 'QR')
plt.xlim(xmax=max(a[:65,0]))
plt.title(r'Accuracy for different $\kappa$',fontsize=15)
plt.xlabel(r'$\kappa$', fontsize=15)
plt.ylabel('Error',fontsize=15)
plt.legend(loc='best', fontsize=15, ncol=2)
fig3.savefig('Tex/figs/Kappa_Accuracy.png', dpi=500)

#%%
import os
os.chdir('/Users/yuki/Dropbox/!Dphil/5_ES/cpp_matrix')

a=[]

with open("same_kappa_new.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
 # format    0 size , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
 #           5 CG_pre count 6 CG_pre time 
 #           7 LU time, 8 LU delx, 9 QR time, 10 QR delx  11 full QR time 12 full QR delx       
 
a=np.array(a,dtype=np.float)
#a = a[:-2,:]       

fig = plt.figure()
plt.loglog(a[:,0], a[:,2],'x-', label = 'SD')
plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')


plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')

plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')
plt.loglog(a[:,0], a[:,6],'bx-', label = 'CG-pre')
plt.loglog(a[:,0], a[:,4],'x-', label = 'CG')



y = np.log(a[:,11]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)     
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
               
               
y = np.log(a[:,7]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='r',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))


y = np.log(a[:,9]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='c',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

y = np.log(a[:,6]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))




plt.xlim(xmax=max(a[:,0]))
plt.title(r'Time of solving for different sizes, $\kappa =2$',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.ylabel('Time in seconds',fontsize=15)
plt.legend(loc='best', ncol=2,fontsize=11)
fig.savefig('sameKappa_time_new.pdf', dpi=500)

#%%

fig3 = plt.figure()

y = np.log(a[:,8]) 
x= np.log(a[:,0])
    
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=3,
               markersize=10,label='Lin. fit m=%s'%round(m,2))
plt.loglog(a[:,0], a[:,8],'x-', label = 'LU')
plt.loglog(a[:,0], a[:,10],'x-', label = 'QR')
plt.loglog(a[:,0], a[:,12],'x-', label = 'QR (expl. Q)')
plt.xlim(xmax=max(a[:,0]))
plt.title(r'Accuracy for different sizes, $\kappa=2$',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.ylabel('Error',fontsize=15)
plt.legend(loc='best', fontsize=12, ncol=2)
fig3.savefig('Tex/figs/sameKappa_Accuracy_new.png', dpi=500)
