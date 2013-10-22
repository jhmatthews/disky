#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python

########################################################
#
#				JM/NSH 10/03/13: ZEPHYR.PY
#
#	models SV wind. written by NSH with some additions
#
#
########################################################
import csv, sys, os, array, warnings
import matplotlib.pyplot as plt
import numpy as np
import sv_sub as sv
import disk_sub as disk
from collections import namedtuple
import astro_help as ah
import geo
import astroconst as consts

RADIAN=57.2959
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
RHO2NH=4.217851e23
PI=3.1415 
HEV=4.13620e-15
H=6.6262e-27
WIEN=5.879e10
YR=3.1556925e7



nogrid_read=True

if (len(sys.argv)>1):
	fname=(sys.argv[1])
	em_fname=fname
	for i in range(len(sys.argv)):
		if sys.argv[i]=="-grid":
			zscale=float(sys.argv[i+1])
			nz=int(sys.argv[i+2])
			nogrid_read=False
			em_fname=fname+'_scale1e'+str(np.log10(zscale))
else:
	fname="test"

fname_to_open=fname+".pf"
specfname=fname+".spec"

print'ZEPHYR: models SV93 wind for parameter file', fname_to_open

sv.input_params(fname_to_open)



#################################################
# print out some simple luminosities

ledd=1.3e38*geo.mstar/MSOL
print "The eddington luminosity for this black hole is ",ledd

medd=ledd/(0.1*C*C)
print "Which implies an accretion rate of ",medd/(MSOL/YR)

edd_frac=(geo.disk_mdot)/medd
print "This geometry is accreting at an eddington fraction of ",edd_frac 

omega=np.cos(geo.sv_thetamin)-np.cos(geo.sv_thetamax)
print "With these angles, the wind will subtend a fraction of ",omega," possible viewing angles"

fev=13.62/HEV
c4ev=64.4/HEV
f2500=C/(2500*1e-8)
f2=C/(800*1e-8)
f1=C/(2000*1e-8)
c4a=C/(1600*1e-8)
c4b=C/(1500*1e-8)
df=1e14

# frequency boundaries
print 'frequency boundaries:',fev,f1,f2

# now do some disk stuff
print geo.mstar, geo.disk_mdot, geo.rstar
tref=disk.tdisk(geo.mstar, geo.disk_mdot, geo.rstar)
rtemp=np.linspace(geo.rstar,geo.disk_radmax,num=1000)
print 'disk tref=',tref

#f,b=disk.spec_disk (1e14,1e17,geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax)
spec_freq,disk_spec=disk.spec_disk2 (1e14,1e18,geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax)
#ax=plt.loglog(f,b)

tdisk=[]
rdisk=[]
spectrum=[]
fmax=[]
lum_ioniz=[]
lum_uv=[]
L_2500=0.0
n_ioniz=[]
cum_ioniz=[]
tot_ioniz=0.0
tot_c4=0.0
tot_nioniz=0.0
peak_nioniz=0.0
for j in range(len(rtemp)-1):
	rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
	r=rdisk[j]/geo.rstar
	area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
	t=(disk.teff(tref,r))
	tdisk.append(t)
	fmax.append(t*WIEN)
	temp_ioniz=0.0
	temp_uv=0.0
	temp_nioniz=0.0
	tempc4=0.0
	i=0
	bnu=1e-98
	while (bnu > 1e-99):
		freq=1e15+i*df
		bnu=ah.planck_nu(t,freq)
		if (freq > fev):
			temp_ioniz=temp_ioniz+bnu
			temp_nioniz=temp_nioniz+(bnu/(H*freq))
		if (freq > f1 and freq < f2):
			temp_uv=temp_uv+bnu
		i=i+1
	L_2500=L_2500+2*PI*ah.planck_nu(t,f2500)*area
	temp_ioniz =PI*temp_ioniz*df*area*2
	temp_uv = PI*temp_uv *df*area*2
	temp_nioniz =PI*temp_nioniz*df*area*2
	tot_ioniz=tot_ioniz+temp_ioniz
	tot_nioniz=tot_nioniz+temp_nioniz
	if j>0:
		cum_ioniz.append(cum_ioniz[j-1]+temp_ioniz)
	else:
		cum_ioniz.append(temp_ioniz)
	if temp_nioniz > peak_nioniz:
		peak_nioniz = tot_nioniz
		r_peak=rdisk[j]
	lum_ioniz.append(temp_ioniz)
	lum_uv.append(temp_uv)
	n_ioniz.append(temp_nioniz)

for i in range(len(cum_ioniz)):
	cum_ioniz[i]=cum_ioniz[i]/tot_ioniz

print "The total ionizing luminosity is ",tot_ioniz," which represents ",tot_nioniz," ionizing photons"
print "The peak occurs at a radius of ",r_peak,"cm from the origin"

if geo.QSO==1:
	const_agn=geo.lum_agn / ((( (2.42e18**( geo.alpha_agn + 1.))) -  (4.84e17**( geo.alpha_agn + 1.0))) /(geo.alpha_agn + 1.0))
	QSO_ioniz=const_agn*((((1e20** (geo.alpha_agn + 1.))) -  (fev** (geo.alpha_agn + 1.))) /(geo.alpha_agn + 1.0))
	QSO_nioniz=(const_agn/H)*((((1e20** ((geo.alpha_agn-1.) + 1.))) -  (fev** ((geo.alpha_agn-1.) + 1.))) /((geo.alpha_agn-1.) + 1.0))
	L_2kev=const_agn*4.84e17**geo.alpha_agn
	Xray_Lum=const_agn*((( (1e20**( geo.alpha_agn + 1.))) -  (1e14**( geo.alpha_agn + 1.0))) /(geo.alpha_agn + 1.0))
	agn_spec=np.zeros(len(spec_freq))
	tot_spec=np.zeros(len(spec_freq))
	for i in range (len(spec_freq)):
		agn_spec[i]=const_agn*spec_freq[i]**geo.alpha_agn
		tot_spec[i]=agn_spec[i]+disk_spec[i]
else:
	QSO_ioniz=0.0
	QSO_nioniz=0.0
	L_2kev=0.0

print "Disk_ioniz=",tot_ioniz
print "QSO_ioniz=",QSO_ioniz
print "QSO_nioniz=",QSO_nioniz
print "F_2500=",f2500
print geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax
print "L_2500=",L_2500,disk.lnu_disk(C/(2500*1e-8),geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax)
print "L_1549=",disk.lnu_disk(C/(1549*1e-8),geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax)
print "L_2kev=",L_2kev
print "L_bol=",ledd*edd_frac
print "Xray_lum=",Xray_Lum
print "ALPHA_OX=",-0.3838*np.log10(L_2500/L_2kev)
f_3560=disk.llamb_disk(3560.0*1e-8,geo.mstar,geo.disk_mdot,geo.rstar,geo.disk_radmax)*1e-8
print "L_3560",f_3560
f_3560=f_3560/(4*consts.PI*consts.PC*consts.PC*100.0) 
print "F_3560 (10pc)=",f_3560
print "Mu=",2.5*np.log10(419.6e-11/f_3560)
print "Mbol=",4.75-2.5*np.log10((ledd*edd_frac)/3.82e33)


print 'Plotting various spectra of e.g. disk, censrc...'
fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Spectrum ("+fname+")")
ax.set_xlabel("Frequency(Hz)")
ax.set_ylabel("L(nu) (erg s-1 )")
ymax=10**float(int(np.log10(np.max(disk_spec)))+1)
ymin=ymax/1e20

#ax.loglog(spec_freq,disk_spec,label="Disk")
if geo.QSO==1:
#	ax.loglog(spec_freq,agn_spec,label="AGN")
	ax.loglog(spec_freq,tot_spec*spec_freq/(4.0*np.pi),label="total")
	ax.loglog(spec_freq,tot_spec*spec_freq,label="total")
	ax.loglog([f2500,4.84e17],[L_2500,L_2kev],'x-',label="Alpha OX")
	ymax=10**float(int(np.log10(np.max(tot_spec)))+1)
	ymin=10**float(int(np.log10(np.min(agn_spec))))
#ax.set_ylim([ymin,ymax])
ax.set_ylim([1e40,1e45])
ax.legend()
fig.savefig(fname+"spectrum.png")
	

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Temperature vs radius ("+fname+")")
ax.set_xlabel("Radius (cm)")
ax.set_ylabel("Temperature (K)")
ax.semilogx(rdisk,tdisk)
ax.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=10)
ax.axvline(x=geo.rstar)
ax.axvline(x=geo.disk_radmax)
fig.savefig(fname+"temp.png")


fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Ionising luminosity vs radius ("+fname+")")
ax.set_xlabel("Radius (cm)")
ax.set_ylabel("Luminosity 13.6eV to infinity (erg s-1)")
ax.text(0.02,0.95,"Central mass="+str(geo.mstar/MSOL)+" solar masses",transform=ax.transAxes)
ax.text(0.02,0.9,"Accretion rate="+str(geo.disk_mdot/(MSOL/YR))+" solar masses/yr",transform=ax.transAxes)
ax.text(0.02,0.85,"This is "+str(edd_frac)+" of the eddington limit",transform=ax.transAxes)
ax.text(0.02,0.75,"Total ionizing luminosity from disk="+str(tot_ioniz),transform=ax.transAxes)
if geo.QSO == 0:
	ax.text(0.02,0.8,"There is no central X-ray source",transform=ax.transAxes)
else:
	ax.text(0.02,0.8,"The Central source adds "+str(QSO_ioniz),transform=ax.transAxes)
ax.text(0.02,0.70,"Number of ionizing photons="+str(tot_nioniz),transform=ax.transAxes)
ax.text(0.02,0.65,"Mdot Wind="+str(geo.wind_mdot)+" solar masses/yr",transform=ax.transAxes)
plt.semilogx(rdisk,lum_ioniz)
plt.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=10)
ax.axvline(x=geo.rstar)
ax.axvline(x=geo.disk_radmax)
ax1=ax.twinx()
ax1.semilogx(rdisk,cum_ioniz)
ax1.set_ylabel("Cumulative ionizing luminosity")
fig.savefig(fname+"lum_ioniz.png")

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Number of ionising photons vs radius ("+fname+")")
ax.set_xlabel("Radius (cm)")
ax.set_ylabel("number of photons with energy over 13.6eV (s-1)")
plt.semilogx(rdisk,lum_ioniz)
plt.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=10)
ax.axvline(x=geo.rstar)
ax.axvline(x=geo.disk_radmax)
fig.savefig(fname+"n_ioniz.png")


fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Luminosity 800-2000Ang vs radius ("+fname+")")
ax.set_xlabel("Radius (cm)")
ax.set_ylabel("Luminosity 800-2000Ang (erg s-1)")
plt.semilogx(rdisk,lum_uv)
ax.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=10)
ax.axvline(x=geo.rstar)
ax.axvline(x=geo.disk_radmax)
fig.savefig(fname+"lum_UV.png")

################################################################


#define grid using input file parameters and user specified zscale, otherwise set to rstar
if nogrid_read:
	zscale=geo.rstar/10.0
	nz=100
xscale=geo.sv_rmin
nx=100
# populates x and z arrays with grid
x,z=geo.get_grid(zscale, xscale, nz, nx, geo.wind_rmax, 1)
xcen,zcen=geo.get_grid(zscale, xscale, nz, nx, geo.wind_rmax, 0)
print 'You are defining a', nx, 'by', nz, 'grid.'
print 'zscale:', zscale, 'xscale:', xscale, 'rmax:', geo.wind_rmax


# NSH VERSION
#x=np.logspace(np.log10(geo.rstar/10.0),np.log10(geo.wind_rmax),num=500)
#z=np.logspace(np.log10(geo.rstar/10.0),np.log10(geo.wind_rmax),num=500)

lx=[]
lz=[]

for i in range(len(x)):
	lx.append(np.log10(x[i]))
for j in range(len(z)):
	lz.append(np.log10(z[j]))


#get density, IP and nh in cells, also emission measure
fig1=plt.figure()
ax1=fig1.add_subplot(211)
ax2=fig1.add_subplot(212)
density=np.empty([len(x),len(z)])
nh=np.empty([len(x),len(z)])
IP=np.empty([len(x),len(z)])
lower_boundx, lower_boundz=0.0,0.0
CEM, EM=np.empty([len(x),len(z)]), np.empty([len(x),len(z)]), 
VEM, VCEM=np.zeros([len(z)]), np.zeros([len(z)])
em_sum=0.0
for  j in range(len(z)):
	if j>0: lower_boundz=z[j-1]
	for i in range(len(x)):
		if i>0: lower_boundx=x[i-1]
		rho_min = geo.sv_rmin + z[j] * np.tan (geo.sv_thetamin)
		rho_max = geo.sv_rmax + z[j] * np.tan (geo.sv_thetamax)
		if (x[i]<=rho_max and x[i]>=rho_min):
			xtest=[x[i],0.0,z[j]]
			density[i,j]=(sv.sv_rho(xtest))
			n_h=density[i,j]*RHO2NH
			nh[i,j]=np.log10(n_h)
			IP[i,j]=np.log10(tot_nioniz/((4*PI*(x[i]*x[i]+z[j]*z[j])*((10.0**nh[i,j]))*C)))
		else:
			density[i,j]=-99
			nh[i,j]=-99
			IP[i,j]=-99
			n_h=00
		volume=(z[j]-lower_boundz)*(x[i]-lower_boundx)
		#note not multiplying by volume
		EM[i,j]=(n_h**2)
		em_sum=em_sum+EM[i,j]
		CEM[i,j]=em_sum
		VEM[j]=VEM[j]+EM[i,j]
	VCEM[j]=em_sum
	


#JM: this is the plot I added
#EM and CEM are emission measures in the actual cells
#VEM is vertical emission measure, summed over the cells at a given z 
#(i.e. the emission measure summed over the row)
#VCEM is cumulative VEM

ax1.loglog(z, VCEM)	
ax2.loglog(z, VEM)			
ax1.set_xlabel("Vertical distance (cm)")
ax2.set_xlabel("Vertical distance (cm)")
ax1.set_ylabel("CEM (n_H^2, cumulative)")
ax2.set_ylabel("EM (n_H^2)")
fig1.savefig(em_fname+"_cellemiss.png")

em1=np.transpose(EM) #transpose the EM array for plotting
fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("emission measure by cell ("+fname+")")
#CP=plt.contourf(x,z,nh1,[0,1,2,3,4,5,6,7,8,9,10,11,12])
#print em1
CP=plt.contourf(lx,lz,np.log10(em1),[5,10,15,17,18,19,20,21,22,23,25])
plt.colorbar(CP)
fig.savefig(em_fname+"_emiss_contour.png")



#done plotting emission measure
print 'Sum of emission measure is:', em_sum

#prepare streamlines
stream_roots=np.logspace(np.log10(1.01*geo.sv_rmin),np.log10(0.99*geo.sv_rmax),num=11)
stream_angles=[]
for root in stream_roots:
	stream_angles.append(sv.sv_theta_wind (root))
stream_dist=np.logspace(10,20,num=1001)


# Plot up velocity, rho, emission measure and IP
fig1=plt.figure()
fig2=plt.figure()
fig3=plt.figure()
fig4=plt.figure()
ax1=fig1.add_subplot(111)
ax2=fig2.add_subplot(111)
ax3=fig3.add_subplot(111)
ax4=fig4.add_subplot(111)
ax1.set_title("Streamline velocity ("+fname+")")
ax2.set_title("Streamline density ("+fname+")")
ax3.set_title("Streamline IP ("+fname+")")
ax4.set_title("Streamline CEM (cumulative emission measure) ("+fname+")")
print 'Plotting streamline velocities, IP, emission measure and density...'
for i in range(len(stream_roots)):
	print stream_roots[i],np.degrees(stream_angles[i])
	stream_vel=[]
	stream_rho=[]
	stream_IP=[]
	stream_CEM=[]
	szarr,sxarr=[],[]
	CEM=0
	for j in range(len(stream_dist)):
		sx=stream_roots[i]+stream_dist[j]*np.sin(stream_angles[i])
		sz=stream_dist[j]*np.cos(stream_angles[i])
		szarr.append(sz)
		sxarr.append(sx)

	for j in range(len(stream_dist)):
		sx=stream_roots[i]+stream_dist[j]*np.sin(stream_angles[i])
		sz=stream_dist[j]*np.cos(stream_angles[i])
		xtest=[sx,0.0,sz]
		stream_rho.append(RHO2NH*sv.sv_rho(xtest))
		speed,temp=sv.sv_velocity(xtest)
		stream_vel.append(np.sqrt(temp[0]*temp[0]+temp[2]*temp[2]))
		stream_IP.append(tot_nioniz/((4*PI*(sx*sx+sz*sz)*((stream_rho[j]))*C)))
		if j<(len(stream_dist)-1): vz=szarr[j+1] - szarr[j]
		if i<(len(stream_dist)-1): vx=sxarr[i+1] - sxarr[i]
		if i==0 and j==0: print 'first cell dimensions check:', vz, vx
		volume=vx * vz
		CEM=CEM + (volume*(stream_rho[j] **2))
		stream_CEM.append(CEM)

	if i==0:
		ax1.loglog(szarr,stream_vel,label="Inner streamline")
		ax2.loglog(szarr,stream_rho,label="Inner streamline")
		ax3.loglog(szarr,stream_IP,label="Inner streamline")
	elif i==len(stream_roots)-1:
		ax1.loglog(szarr,stream_vel,label="Outer streamline")
		ax2.loglog(szarr,stream_rho,label="Outer streamline")
		ax3.loglog(szarr,stream_IP,label="Outer streamline")
	else:
		ax1.loglog(szarr,stream_vel)
		ax2.loglog(szarr,stream_rho)
		ax3.loglog(szarr,stream_IP)
		ax4.loglog(szarr,stream_CEM)


ax1.legend(loc='upper left')
ax1.text(0.03,0.8,'Acceleration length='+str(geo.sv_r_scale)+'cm',transform=ax1.transAxes)
ax1.text(0.03,0.76,'Acceleration exponent='+str(geo.sv_alpha),transform=ax1.transAxes)
ax1.text(0.03,0.72,'V infinity ='+str(geo.sv_v_infinity)+' Vesc',transform=ax1.transAxes)
ax1.text(0.03,0.68,'Streamline bunching ='+str(geo.sv_gamma),transform=ax1.transAxes)
ax1.text(0.03,0.64,'Inner streamline root ='+str(geo.sv_rmin),transform=ax1.transAxes)
ax1.text(0.03,0.60,'Inner streamline root ='+str(geo.sv_rmax),transform=ax1.transAxes)
ax2.legend(loc='lower left')
ax2.text(0.03,0.5,'Accretion rate='+str(geo.disk_mdot/(MSOL/YR))+'Msun/yr',transform=ax2.transAxes)
ax2.text(0.03,0.46,'Mdot R exponent='+str(geo.sv_lambda),transform=ax2.transAxes)
ax2.text(0.03,0.42,'Streamline bunching ='+str(geo.sv_gamma),transform=ax2.transAxes)
ax2.text(0.03,0.64,'Inner streamline root ='+str(geo.sv_rmin),transform=ax2.transAxes)
ax2.text(0.03,0.60,'Inner streamline root ='+str(geo.sv_rmax),transform=ax2.transAxes)
ax3.legend(loc='upper left')
ax1.set_xlabel("Vertical distance along streamline(cm)")
ax2.set_xlabel("Vertical distance along streamline(cm)")
ax3.set_xlabel("Vertical distance along streamline(cm)")
ax4.set_xlabel("Vertical distance along streamline(cm)")
ax1.set_ylabel("Poloidal velocity(cms-1)")
ax2.set_ylabel("Hydrogen number density(cm-3)")
ax3.set_ylabel("Ionization parameter")
ax4.set_ylabel("CEM (rho^2 x volume)")
fig1.savefig(fname+"streamline_vel.png")
fig2.savefig(fname+"streamline_rho.png")
fig3.savefig(fname+"streamline_IP.png")
fig4.savefig(fname+"streamline_CEM.png")

################################################################


sightline_array=[71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0]
total_column=[]
for sightline in sightline_array:
	sightline=sightline/RADIAN
	test_sightline=(PI/2.0)-sightline
	launch=1e15
	test_thetamin=(PI/2.0)-geo.sv_thetamin
	test_thetamax=(PI/2.0)-geo.sv_thetamax
	if (test_sightline < test_thetamin):
		z1=(geo.sv_rmin-launch)/(1/np.tan(test_sightline)-1/np.tan(test_thetamin))
		x1=launch+z1/np.tan(test_sightline)
		if (test_sightline < test_thetamax):
			z2=(geo.sv_rmax-launch)/(1/np.tan(test_sightline)-1/np.tan(test_thetamax))
			x2=launch+z2/np.tan(test_sightline)
			if (z2 > geo.wind_rmax):
				z2=geo.wind_rmax
			if (x2 > geo.wind_rmax):
				x2=geo.wind_rmax		
		else:
 			x2=geo.wind_rmax
			z2=geo.wind_rmax
	print 'At sightline ',np.degrees(sightline),' we go from',x1,z1,' to ',x2,z2
	xsl=[]
	zsl=[]
	lxsl=[]
	lzsl=[]
	lsl=[]
	rhosl=[]
	col=0.0
	column=[]
	dl=np.sqrt(((z2-z1)*(z2-z1))+((x2-x1)*(x2-x1)))
	dl=dl/10000.0
	for i in range(10000):
		xsl.append(x1+i*dl*np.cos(test_sightline))
		lxsl.append(np.log10(xsl[i]))
		zsl.append(z1+i*dl*np.sin(test_sightline))
		lzsl.append(np.log10(zsl[i]))
		lsl.append(np.sqrt(xsl[i]*xsl[i]+zsl[i]*zsl[i]))
		rhosl.append((sv.sv_rho([xsl[i],0.0,zsl[i]])*RHO2NH))
		col=col+(rhosl[i]*dl)
		column.append(col)
	print "Cumulative hydrogen density along sightline ",np.degrees(sightline)," to origin is ",column[-1]
	total_column.append(column[-1])

tot_col=np.array(total_column)
tau=tot_col*consts.THOMPSON


out=open("tau_data.dat",'w')
for i in range(len(tau)):
	out.write(str(sightline_array[i])+' '+str(tau[i])+'\n')

out.close()


fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(sightline_array,tau)
ax.set_ylabel("Electron scattering optical depth through wind to origin")
ax.set_xlabel("Sightline")
fig.savefig(fname[:-3]+"_tau_sightline.eps",bbox_inches='tight')
fig.savefig(fname[:-3]+"_tau_sightline.png",bbox_inches='tight')


'''
fig=plt.figure()
ax1=fig.add_subplot(111)
ax.loglog(lsl,rhosl)
ax.set_xlabel("Distance along photon path (cm)")
ax.set_ylabel("Hydrogen density")
ax1=ax.twinx()
ax1.loglog(lsl,column)
ax1.set_ylabel("Cumlative neutral hydrogen column")
fig.savefig(fname+"passage.png")
'''


nh1=np.transpose(nh)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Hydrogen density ("+fname+")")
CP=plt.contourf(x,z,nh1,[0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=2)
plt.plot([geo.rstar,geo.disk_radmax],[0,0],linewidth=10)
plt.plot(xsl,zsl)
plt.colorbar(CP)
fig.savefig(fname+"linear.png")

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Hydrogen density ("+fname+")")
ax.axis([0,1e16,0,1e16])
CP=plt.contourf(x,z,nh1,[0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.plot([geo.sv_rmin,geo.sv_rmax],[0,0],linewidth=2)
plt.plot([geo.rstar,geo.disk_radmax],[0,0],linewidth=10)
plt.plot(xsl,zsl,color='r')
plt.plot([launch,xsl[0]],[0,zsl[0]],color='r')
plt.colorbar(CP)
fig.savefig(fname+"linear_zoom.png")

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Hydrogen density ("+fname+")")
CP=plt.contourf(lx,lz,nh1,[0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.plot([np.log10(geo.sv_rmin),np.log10(geo.sv_rmax)],[lz[0],lz[0]],linewidth=2)
plt.plot([np.log10(geo.rstar),np.log10(geo.disk_radmax)],[lz[0],lz[0]],linewidth=10)
plt.plot(lxsl,lzsl)
plt.colorbar(CP)
fig.savefig(fname+"log.png")

IP1=np.transpose(IP)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_title("Hydrogen density ("+fname+")")
CP=plt.contourf(lx,lz,IP1)
plt.plot([np.log10(geo.sv_rmin),np.log10(geo.sv_rmax)],[lz[0],lz[0]],linewidth=2)
plt.plot([np.log10(geo.rstar),np.log10(geo.disk_radmax)],[lz[0],lz[0]],linewidth=10)
plt.colorbar(CP)
fig.savefig(fname+"IP.png")
