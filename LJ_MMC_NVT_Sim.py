##########################
# NVT Metropolis Monte-Carlo simulation code
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
		
		return 0.0
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
#Computes the total pair energy
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
	prob = math.exp(-Beta*dE)
	return prob




#------Params to set-------------------------

#Number of atoms
N = 17
# lj epsilon 
epsilon = 0.238
#lj sigma
sigma = 3.405
#Mass
mass = 39.95
#Number of trial moves
MCit = 2000000
#Number of initial equilibration steps
Eqit = 100000
#Step interval to collect data after equilibration
Dit = 10
#Step interval to output running averages to screen
Scit = 100000
#Boundary Conditions -- using hard sphere so need a center and radius
center = 0.0, 0.0, 0.0
radius = 6.6 * sigma
#Temperature
Temp = 74.0
#Boltzmann Constant--
# make sure energy portion of units match the epsilon energy units 
# for reduced units use 1.0
kb = 0.00198 
# if using real units this is the corresponding reduced temperature
rt  = kb*Temp/epsilon
#Trial move step size - factor and then actual
factor = 2.5
#This is approximate formula- acceptance ratio should be 
# between 0.3 to 0.5 , so factor may need to be adjusted 
step = factor*radius*rt


# --------------------------------------------------------


#seed random number generator - a=None uses default clock time
random.seed(a=None)
print " "
print "Hello! I will run a Monte Carlo Simulation of a homogenous LJ system!"
print " "
print "Parameters have been set to: "
print "Number of atoms: ", N, " epsilon: ", epsilon, " sigma: ", sigma, " mass: ", mass
print "The Temperature: ", Temp, " Number of Trial Moves: ", MCit, "Number of Equilibration Steps: ", Eqit
print "The interval to collect data after equilibration: ", Dit, " trial move step size: ", step
print "Spherical Boundary Radius: ", radius, " center: ", center
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
	

#Initialize the boundary
boundary = HardSphericalWall(center, radius)

#Set some initial coordinates - 
# assign randomly inside a cube inscribed in our spherical boundary
box = radius*1.40
print "Assigning Initial Coordinates..."
print " "
for index in xrange(0, N):
	c.atom[index].x = center[0]+ box * (random.random()-0.5)
	c.atom[index].y = center[1]+ box *  (random.random()-0.5)
	c.atom[index].z = center[2]+ box *  (random.random()-0.5)
	
#Compute the Energy of the Initial configuration
Energy = TotalPairEnergy(c)
print "Energy of initial configuration: ", Energy
print " "
#Get ready for sampling

#Initialize all of the needed tracking of observables
#Average energy
Estat = RunningStats()
#Average energy squared
Esqstat = RunningStats()
#heat capacity
Cvstat = RunningStats()
#radius of gyration
Rgstat = RunningStats()

#Beta 
Beta = 1.0 / (kb * Temp)

acceptcount = 0
print "Beginning the sampling loop. Please wait while I work..."
print " "
#start the sampling loop
for i in xrange(0, MCit):
	#save energy before trial move
	Enb = Energy
	#Pick a random atom - by index
	Aran = int (random.random() * N)
	#Save its coordinates
	xold = c.atom[Aran].x
	yold = c.atom[Aran].y
	zold = c.atom[Aran].z
	#Perform trial move - random translation in 3d
	xnew = c.atom[Aran].x = xold +  step * (random.random()-0.5)
	ynew = c.atom[Aran].y = yold +  step * (random.random()-0.5)
	znew = c.atom[Aran].z = zold +  step * (random.random()-0.5)
	ar = c.atom[Aran]
	#Apply Acceptance Criteria
	
	#First check boundary condition -- 1 = pass, 0 = fail
	btag = boundary.CheckIfInside(xnew, ynew, znew)
	
	etag = 0
	
	#if boundary condition passed compute change in energy and test acceptance
	if btag == 1:
		dE = DeltaE(c, Aran, xold, yold, zold)
		Energy += dE
		prob = ConfProb(dE, Beta)
		ranprob = random.random()
		if ranprob < prob:
			etag = 1

	ttag = 0

	if btag == 0 or etag == 0:
		c.atom[Aran].x = xold
		c.atom[Aran].y = yold
		c.atom[Aran].z = zold
		Energy = Enb
	else:
		acceptcount += 1
	#collect data 
	if i> Eqit and i % Dit == 0:
		Estat.Push(Energy)
		Esqstat.Push(Energy * Energy)	
		Rg = c.RadGyr()
		Rgstat.Push(Rg)
		Cv = (Esqstat.Mean() - (Estat.Mean()*Estat.Mean()))*(Beta/Temp)
		Cvstat.Push(Cv)


	if i % Scit == 0:
		print "## Iteration ", i, " current averages: "
		print "Energy: ", Estat.Mean(), " Heat Capacity: ", Cvstat.Mean()
		print "Radius of Gyration: ", Rgstat.Mean()
		print " "

# End of sampling loop
ar = float(acceptcount)/float(MCit)
ccv = (Esqstat.Mean() - (Estat.Mean()*Estat.Mean()))*(Beta/Temp)
print "Simulation Complete!"
print "Trial Move Acceptance ratio: ", ar
print "Here are some final average values-"
print "Energy: ", Estat.Mean()
print "Energy Squared: ", Esqstat.Mean() 
print "Heat Capacity: ", Cvstat.Mean()
print "Radius of Gyration: ", Rgstat.Mean()
#print "Final Instantaneous Heat Capacity: ", ccv
print " "
print " Thank You and Have a Nice Day!"
print " "
print "... End of Line ..."



