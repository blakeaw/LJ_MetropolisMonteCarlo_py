##########################
# NPT Metropolis Monte-Carlo simulation code
# simulates homogenous system of Lennard-Jones (LJ) particles 
# interacting via the traditional 12-6 LJ pairwise additive potential
#
#
#  
# Written by Blake Wilson
# email: blake.wilson@utdallas.edu
#
# Please e-mail me if you have comments, suggestions, or to report errors/bugs.
#
#########################


# import some libraries
import os
import math
import random


#Define some classes

#Define an atom class
class Atom:
	x = y = z = 0.0
	eps = sig = 0.0
	mass = 0.0


# Define a configuration class - configuration is a collection of
# atoms
class Configuration:
	def __init__(self, na):
		self.natoms = na
		self.COM = [0.0, 0.0, 0.0]
		self.atom = []
		for i in xrange(0, na):
			self.atom.append(Atom())
	
	def CalcCOM(self):
		M = 0.0
		sumx = sumy = sumz = 0.0
		for i in xrange(0, self.natoms):
			m = self.atom[i].mass
			sumx += self.atom[i].x * m
			sumy += self.atom[i].y * m
			sumz += self.atom[i].z * m
			M += m

		Cx = sumx/M
		Cy = sumy/M
		Cz = sumz/M
		self.COM[0]=Cx
		self.COM[1]=Cy
		self.COM[2]=Cz
		return
		
	def RadGyr(self):
		sumgyr = 0.0
		M = 0.0
		self.CalcCOM()
		comx = self.COM[0]
		comy = self.COM[1]
		comz = self.COM[2]
		for i in xrange(0, self.natoms):
			M += self.atom[i].mass
			mc = self.atom[i].mass
			rx = self.atom[i].x
			ry = self.atom[i].y
			rz = self.atom[i].z
			sgx = (rx - comx)**2.0
			sgy = (ry - comy)**2.0
			sgz = (rz - comz)**2.0

			sumgyr += mc*(sgx+sgy+sgz)

		Rsq = sumgyr / M
		R = math.sqrt(Rsq)
		return R

	def ScaleCoordinates(self, factor):
		for i in xrange(self.natoms):
			self.atom[i].x*=factor
			self.atom[i].y*=factor
			self.atom[i].z*=factor

		return

#Running Statistics
class RunningStats:
	n = 0
	Mnold = Mnnew = Snold = Snnew = 0.0
	def Push(self, val):
		self.n += 1
		if self.n == 1:
			self.Mnold = val
			self.Snold = 0.0
		else:
			self.Mnnew = self.Mnold + (val - self.Mnold)/((float)(self.n));
			self.Snnew = self.Snold + (val - self.Mnold)*(val-self.Mnnew);
			self.Mnold = self.Mnnew;
			self.Snold = self.Snnew;
		return
	def Mean(self):
		if self.n > 0:
			return self.Mnnew
		else:
			return 0.0
	def Variance(self):
		if self.n > 1:
			vary = self.Snnew/(float(self.n)-1.0)
			return vary
		else:
			return 0.0
	def Deviation(self):
		dev = sqrt(self.Variance())
		return dev
	def Reset():
		self.n = 0
		return

# A boundary condition class 
class HardSphericalWall:
	def __init__(self, cent, rad):
		self.center = cent
		self.radius = rad

	def CheckIfInside(self, x, y, z):
		dx = x - self.center[0]
		dy = y - self.center[1]
		dz = z - self.center[2]
		r = math.sqrt(dx*dx + dy*dy + dz*dz)
		if r < self.radius:
			return 1
		else:
			return 0
	def SetRad(self, l):
		self.radius = l
		return
	def ResetSizebyVolume(self, Vol):
		r = (3.0*Vol/(4.0*math.pi))**(1.0/3.0)
		self.radius = r;
	def PrintDetails(self):
		print "Hard Spherical Wall Boundary"
		print "Spherical Boundary Radius: ", self.radius, " center: ", self.center
		print "Volume: ", self.ComputeVolume()
		return

# Another boundary condition class 
class HardCubicBox:
	def __init__(self, cent, l):
		self.center = cent
		self.bx = l
		self.by = l
		self.bz = l
		#lower box bounds
		self.bxl = self.center[0]-self.bx/2.0
		self.byl = self.center[1]-self.by/2.0
		self.bzl = self.center[2]-self.bz/2.0
		#upper box bounds
		self.bxh = self.center[0]+self.bx/2.0
		self.byh = self.center[1]+self.by/2.0
		self.bzh = self.center[2]+self.bz/2.0
		return

	def CheckIfInside(self, x, y, z):
		dx = x - self.center[0]
		dy = y - self.center[1]
		dz = z - self.center[2]
		#x bound
		if dx<self.bxl or dx>self.bxh :
			return 0
		#y bound
		elif dy<self.byl or dy>self.byh :
			return 0
		#z bound
		elif dz<self.bzl or dz>self.bzh :
			return 0
		else:
			return 1

	def ComputeVolume(self):
		return self.bx*self.by*self.bz

	def SetBox(self, l):
		self.bx = l
		self.by = l
		self.bz = l
		#lower box bounds
		self.bxl = self.center[0]-self.bx/2.0
		self.byl = self.center[1]-self.by/2.0
		self.bzl = self.center[2]-self.bz/2.0
		#upper box bounds
		self.bxh = self.center[0]+self.bx/2.0
		self.byh = self.center[1]+self.by/2.0
		self.bzh = self.center[2]+self.bz/2.0
		return
	def ResetSizebyVolume(self, Vol):
		l = Vol**(1.0/3.0)
		self.bx = l
		self.by = l
		self.bz = l
		#lower box bounds
		self.bxl = self.center[0]-self.bx/2.0
		self.byl = self.center[1]-self.by/2.0
		self.bzl = self.center[2]-self.bz/2.0
		#upper box bounds
		self.bxh = self.center[0]+self.bx/2.0
		self.byh = self.center[1]+self.by/2.0
		self.bzh = self.center[2]+self.bz/2.0

	def PrintDetails(self):
		print "Hard Cubic Box Boundary"
		print "Cubic Box Size: ", self.bx, " center: ", self.center
		print "Volume: ", self.ComputeVolume()
		return
	
#Define Functions

#Pair Energy between two atoms using LJ12-6
def LJPairEnergy(atom1, atom2):
	eps = atom1.eps
	sig = atom1.sig
	rcs = (3.0*sig)**2.0	

	dx = atom1.x - atom2.x
	dy = atom1.y - atom2.y
	dz = atom1.z - atom2.z
	
	rs = dx**2.0 + dy**2.0 + dz**2.0
	
	if rs<rcs:
		E = 4.0*eps*( (sig/rs)**6.0 - (sig/rs)**3.0 )
	else:
		E = 0.0

	return E	
#Pair Energy between an atom and old coordinates of other using LJ12-6
def LJPairEnergyO(atom, xo, yo, zo):
	eps = atom.eps
	sig = atom.sig
	rcs = (3.0*sig)**2	

	dx = atom.x - xo
	dy = atom.y - yo
	dz = atom.z - zo
	rs = dx**2.0 + dy**2.0 + dz**2.0
	if rs<rcs:
		E = 4.0*eps*( (sig/rs)**6.0 - (sig/rs)**3.0 )
	else:
		E = 0.0

	return E	
#Computes the total pair energy of configuration object using LJ12-6
def TotalPairEnergy(c):
	N = c.natoms
	Et = 0.0
	for i in xrange(0, N-1):
		for j in xrange(i+1, N):
			
			Et += LJPairEnergy(c.atom[i], c.atom[j])
	return Et
#computes the difference in energy from moving an atom
def DeltaE(c, aran, xo, yo, zo):
	Et = 0.0
	Eto = 0.0
	for i in xrange(0, c.natoms):
		if i != aran:
			Et += LJPairEnergy(c.atom[aran], c.atom[i])
			Eto += LJPairEnergyO(c.atom[i], xo, yo, zo)

	dE = Et - Eto
	return dE
#computes the probability or probability ratio using Boltzmann weight
def ConfProb(dE, Beta):
	#print "value: ", dE, " beta: ", Beta
	prob = math.exp(-Beta*dE)
	return prob




#------Params to set-------------------------

#Number of atoms
N = 17
# lj epsilon 
epsilon = 1.0
#lj sigma
sigma = 1.0
#Mass
mass = 1.0
#Number of trial moves
MCit = 4000000
#Number of MCit to devote to initial equilibration steps
Eqit = 2000000
# frequency of trial moves to do volume deform
vit = 2*N
#Step interval to collect data after equilibration
Dit = 2*vit
#Step interval to output running averages to screen
Scit = 100000
#Boundary Conditions -- using hard cubic box so need a center and box length
center = 0.0, 0.0, 0.0
box = 8.0
#Temperature
Temp = 2.0
#Pressure
Press = 1.0


#Boltzmann Constant--
# make sure energy portion of units match the epsilon energy units 
# for reduced units use 1.0
#kb = 0.0019872041 
kb=1.0
# if using real units this is the corresponding reduced temperature
rt  = kb*Temp/epsilon
#Trial move step size - factor and then actual
factor = 0.15
#This is approximate formula- acceptance ratio should be 
# between 0.3 to 0.5 , so factor may need to be adjusted 
step = factor*box*rt
vstep = factor*rt*box**3

# --------------------------------------------------------


#seed random number generator - a=None uses default clock time
random.seed(a=None)

#Initialize the boundary
boundary = HardCubicBox(center, box)

print " "
print "Hello! I will run a Monte Carlo Simulation of a homogenous LJ system!"
print " "
print "Parameters have been set to: "
print "Number of atoms: ", N, " epsilon: ", epsilon, " sigma: ", sigma, " mass: ", mass
print "The Temperature: ", Temp, " Number of Trial Moves: ", MCit, "Number of Equilibration Steps: ", Eqit
print "The interval to collect data after equilibration: ", Dit, " trial move step size: ", step
print "The Pressure: ", Press
#print "Spherical Boundary Radius: ", radius, " center: ", center
boundary.PrintDetails()
print "reduced temperature equivalent: ", rt
print " "

#Initialize config with N atoms
c = Configuration(N)
#set the epsilon, sigma, and mass parameters -- only works for homogenous system
print "Assigning atom parameters..."
print " "
for index in xrange(0, N):
	c.atom[index].eps = epsilon
	c.atom[index].sig = sigma
	c.atom[index].mass = mass
	



#Set some initial coordinates - 
# assign randomly inside a cube 
#box = radius*1.40
print "Assigning Initial Coordinates..."
print " "
Energy = 100.0
while Energy>10.0:
	for index in xrange(0, N):
		c.atom[index].x = center[0]+ box * (random.random()-0.5)
		c.atom[index].y = center[1]+ box *  (random.random()-0.5)
		c.atom[index].z = center[2]+ box *  (random.random()-0.5)
	Energy = TotalPairEnergy(c)
	#print "Energy of initial configuration: ", Energy, " xval: ", x
#Compute the Energy of the Initial configuration
#Energy = TotalPairEnergy(c)
print "Energy of initial configuration: ", Energy
#Compute volume
Volume = boundary.ComputeVolume()
print "Volume of initial configuration: ", Volume
#Enthalpy
H = Energy + Press*Volume
print "Enthalpy of initial configuration: ", H
print " "
#Get ready for sampling

#Initialize all of the needed tracking of observables
#Average energy
Estat = RunningStats()
#Average enthalpy
Hstat = RunningStats()
#Average Enthalpy Squared
Hsqstat = RunningStats()
#heat capacity
Cpstat = RunningStats()
#radius of gyration
Rgstat = RunningStats()
#Volume
Vstat = RunningStats()

#Beta 
Beta = 1.0 / (kb * Temp)

#trial move acceptance counters
acceptcount = 0
vacceptcount = 0
xmovecount = 0
vmovecount =0

print "Beginning the sampling loop. Beta is: ", Beta, ". Please wait while I work..."
print " "
factor = 1.0
#start the sampling loop
for i in xrange(0, MCit):
	#save energy before trial move
	Enb = Energy
	Vnb = boundary.ComputeVolume()
	Hnb = Energy + Press*Vnb
	#Pick a random atom - by index
	Aran = int (random.random() * N)
	#Save its coordinates
	xold = c.atom[Aran].x
	yold = c.atom[Aran].y
	zold = c.atom[Aran].z
	
	btag = 1
	vtag = 0
	xtag = 0
	vprob=1.0
	# Volume deformation trial move
	if i%vit==0:
		vmovecount+=1
		#lnVnew = math.log(Vnb) + (random.random()-0.5)*vstep
		Vnew = Vnb + (random.random()-0.5)*vstep
		# prevent Volume from going below zero
		while Vnew<0.0:
			Vnew = Vnb + (random.random()-0.5)*vstep
		Vscale = Vnew/Vnb
		ranprob = random.random()
		#vprob = 1.0
		#if Vscale<1.0:
		dN = float(N)
		vprob = Vscale**dN
		
		#reset new test boundary size 
		boundary.ResetSizebyVolume(Vnew)
		factor = Vscale**(1.0/3.0)
		c.ScaleCoordinates(factor)
		Energy = TotalPairEnergy(c)
		vtag=1

	else:
		xmovecount+=1
		#Perform trial move - random translation in 3d
		xnew = c.atom[Aran].x = xold +  step * (random.random()-0.5)
		ynew = c.atom[Aran].y = yold +  step * (random.random()-0.5)
		znew = c.atom[Aran].z = zold +  step * (random.random()-0.5)
		#ar = c.atom[Aran]
		#Apply Acceptance Criteria
	
		#First check boundary condition -- 1 = pass, 0 = fail
		btag = boundary.CheckIfInside(xnew, ynew, znew)
		xtag = 1
		
	
		#if boundary condition passed compute change in energy and test acceptance
		if btag == 1:
			dE = DeltaE(c, Aran, xold, yold, zold)
			Energy += dE
			
	etag = 0
	ttag = 0
	V = boundary.ComputeVolume()
	H = Energy + Press*V
	dH = H - Hnb
	prob = ConfProb(dH, Beta)
	ranprob = random.random()
	if prob*vprob > ranprob:
		etag = 1
	# if failed trial move
	if btag == 0 or etag == 0:
		#was volume deformation
		if vtag==1:
			Unscale =1.0/factor
			c.ScaleCoordinates(Unscale)
			Energy = Enb
			H = Hnb
			V = Vnb
			boundary.ResetSizebyVolume(V);
		elif xtag==1:
			# was coordinate deform
			c.atom[Aran].x = xold
			c.atom[Aran].y = yold
			c.atom[Aran].z = zold
			Energy = Enb
			H = Hnb
			V = Vnb
	
	else:
		#trial move was accepted
		if vtag==1:
			vacceptcount += 1
		elif xtag==1:
			acceptcount += 1

	#collect data 
	if i> Eqit and i % Dit == 0:
		Estat.Push(Energy)
		Hstat.Push(H)
		Hsqstat.Push(H * H)	
		Rg = c.RadGyr()
		Rgstat.Push(Rg)
		Cp = (Hsqstat.Mean() - (Hstat.Mean()*Hstat.Mean()))*(Beta/Temp)
		Cpstat.Push(Cp)
		Vstat.Push(V)

	if i % Scit == 0:
		print "## Iteration ", i, " current averages: "
		if i<Eqit:
			print "Equilibrating..."
		else:
			print "Energy: ", Estat.Mean(), "Enthalpy: ", Hstat.Mean(), " Heat Capacity: ", Cpstat.Mean(), " Enthalpy Squared: ", Hsqstat.Mean()
			print "Radius of Gyration: ", Rgstat.Mean(), " Volume ", Vstat.Mean()
		print " "

# End of sampling loop
xar = float(acceptcount)/float(xmovecount)
var = float(vacceptcount)/float(vmovecount)
ccv = (Hsqstat.Mean() - (Hstat.Mean()*Hstat.Mean()))*(Beta/Temp)
print "Simulation Complete!"
print "Trial Move Acceptance ratios: "
print "- Coordinate Deform: ", xar
print "- Volume Deform: ", var
print "Here are some final average values-"
print "Energy: ", Estat.Mean()
print "Enthalpy ", Hstat.Mean()
print "Enthalpy Squared: ", Hsqstat.Mean() 
print "Heat Capacity: ", Cpstat.Mean()
print "Radius of Gyration: ", Rgstat.Mean()
print "Volume: ", Vstat.Mean()
print "Final Instantaneous Heat Capacity: ", ccv
print " "
print " Thank You and Have a Nice Day!"
print " "
print "... End of Line ..."



