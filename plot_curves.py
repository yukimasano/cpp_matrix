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
nn = "x"
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
plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')
plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')


plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')

plt.loglog(a[:,0], a[:,4],'yx-', label = 'CG')

plt.loglog(a[:,0], a[:,2],'gx-', label = 'SD')

plt.loglog(a[:,0], a[:,6],'bx-', label = 'CG-pre')

plt.loglog(a[:,0], a[:,14],'x-', label = 'Jacobi')

plt.loglog(a[:,0], a[:,16],'x-', label = 'GS')
plt.loglog(a[:,0], a[:,18],'x-', label = 'SOR1.5')
plt.loglog(a[:,0], a[:,20],'x-', label = 'SOR0.5')





y = np.log(a[:,11]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)     
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

y = np.log(a[:,9]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)
y100= np.exp(m*np.log(a[-1,0] )+b)

plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='c',linewidth=2,
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

plot1 = plt.loglog([a[0,0], a[-1,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))


plt.xlim(xmax=max(a[:,0]+10))
plt.title(r'Time of solving for different sizes, $\kappa =10$',fontsize=15)
plt.xlabel('Size N',fontsize=15)
plt.ylabel('Time in seconds',fontsize=15)
#plt.legend(loc='best', ncol=2,fontsize=11)
plt.gca().legend(loc='center left', ncol=2, bbox_to_anchor=(1, 0.5))
fig.savefig('%s.pdf'%nn,bbox_inches='tight')


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
with open("same_size.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
     # format    0 kappa , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
     #           5 CG_pre count 6 CG_pre time 
     #           7 LU time, 8 LU delx, 9 QR time, 10 QR delx  11 full QR time 12 full QR delx       
             
a=np.array(a,dtype=np.float)

####
plt.loglog(a[:,0], a[:,11],'mx-', label = 'QR (expl. Q)')
plt.loglog(a[:,0], a[:,9],'cx-', label = 'QR')


plt.loglog(a[:,0], a[:,7],'rx-', label = 'LU')

plt.loglog(a[:,0], a[:,4],'yx-', label = 'CG')

plt.loglog(a[:,0], a[:,2],'gx-', label = 'SD')

plt.loglog(a[:,0], a[:,6],'bx-', label = 'CG-pre')



y = np.log(a[:,11]) 
x= np.log(a[:,0])
m,b = np.polyfit(x, y, 1)
y10= np.exp(m*np.log(a[0,0]) +b)     
y100= np.exp(m*np.log(a[-(sum(np.isnan(a[:,11]))),0] )+b)

plt.loglog([a[0,0], a[-(sum(np.isnan(a[:,11]))),0]],[y10,y100] ,'--', marker='.',color='m',linewidth=2,
               markersize=10,label='Lin. fit m=%s'%round(m,2))

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

#y = np.log(a[:65,6]) 
#x= np.log(a[:65,0])
#    
#m,b = np.polyfit(x, y, 1)
#y10= np.exp(m*np.log(a[0,0]) +b)
#y100= np.exp(m*np.log(a[65,0] )+b)
#
#fig3 = plt.figure()
#plt.loglog([a[0,0], a[65,0]],[y10,y100] ,'--', marker='.',color='b',linewidth=3,
#               markersize=10,label='Lin. fit m=%s'%round(m,2))
#plt.loglog(a[:65,0], a[:65,6],'x-', label = 'LU')
#plt.loglog(a[:65,0], a[:65,8],'x-', label = 'QR')
#plt.xlim(xmax=max(a[:65,0]))
#plt.title(r'Accuracy for different $\kappa$',fontsize=15)
#plt.xlabel(r'$\kappa$', fontsize=15)
#plt.ylabel('Error',fontsize=15)
#plt.legend(loc='best', fontsize=15, ncol=2)
#fig3.savefig('Tex/figs/Kappa_Accuracy.png', dpi=500)

#%%
import os
import csv
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/Users/yuki/Dropbox/!Dphil/5_ES/cpp_matrix')


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
    plt.title(r'Time of solving for different sizes, $\kappa =10$',fontsize=15)
    plt.xlabel('Size N',fontsize=15)
    plt.ylabel('Time in seconds',fontsize=15)
    plt.legend(loc='best', ncol=2,fontsize=11)
    fig.savefig('%s.pdf'%nn)
#%%
a=[]
with open("same_kappa2_new.txt") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=",")
    for line in tsvreader:
        a.append( line)
 # format    0 size , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time, 
 #           5 CG_pre count 6 CG_pre time 
 #           7 LU time, 8 LU delx, 9 QR time, 10 QR delx  11 full QR time 12 full QR delx       
 
a=np.array(a,dtype=np.float)
#a = a[:-2,:]       

fig = plt.figure()
plt.loglog(a[:,0], a[:,1],'x-', label = 'SD')
plt.loglog(a[:,0], a[:,5],'bx-', label = 'CG-pre')
plt.loglog(a[:,0], a[:,3],'x-', label = 'CG')








plt.xlim(xmax=max(a[:,0]))
plt.title(r'Time of solving for different sizes, $\kappa =10$',fontsize=15)
plt.xlabel('Number of iters',fontsize=15)
plt.ylabel('Time in seconds',fontsize=15)
plt.legend(loc='best', ncol=2,fontsize=11)
fig.savefig('sameKappa10_count_new.pdf', dpi=500)
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
