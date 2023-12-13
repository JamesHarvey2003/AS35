#!/usr/bin/env python3

# University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Isochrone Fitter

# Import Python Libraries
import math
import numpy as np
import matplotlib
matplotlib.use('GTK3Agg')	# Added to force use of X display for forwarding ; FC 2020/04/28
import matplotlib.pyplot as plt
import curses
import time
import sys
import getopt


class CMDdata() : 	# CMD data class
	
	def __init__(self,name) : 
		self.name=name
		self.v=[]
		self.bv=[]		
	
	def loaddata(self,datafile) : 
		try : 
			f = open(datafile,"r")
			self.name = datafile
			for line in f : 
				if(line[0]=='#') : 	# Comment lines
					continue
				data = line.split()
				# Modified to read in RA, DEC, B, V format produced for AS32 scripts.
				self.add_point(data[2],data[3])
		except : 
			print ("Error loading data CMD from file " + datafile)
		
	def add_point(self,bmag,vmag) : 
		self.v.append(float(vmag))
		self.bv.append(float(bmag)-float(vmag))
	
	def plot(self) : 
		plt.clf()
		plt.plot(self.bv, self.v,'.',label=self.name)
		plt.xlim(-2.0,3.0)
		plt.gca().invert_yaxis()
		plt.xlabel('B-V colour')
		plt.ylabel('V magnitude')
		
class CMDmodel(CMDdata) : 	# CMD model class
	
	def __init__(self,age,Z) :
		self.name="MODEL"
		self.age=age
		self.Z=Z
		self.v=[]
		self.bv=[]
		self.R=3.0	#	Extinction coefficient in E(B-V) to A(V)
	
	def add_point(self,bmag,vmag) :
#		print bmag,vmag,float(bmag)-float(vmag)
		self.v.append(float(vmag))
		self.bv.append(float(bmag)-float(vmag))	

	def plot(self,dist=0.0,ext=0.0) :
		vcor=dist + (ext * self.R)
		plt.plot(np.array(self.bv) + ext,np.array(self.v) + vcor,'-',lw=2)
	
		
def loadmodels(modelfile="PadovaCMD.dat") :

	CMDs = []
	
	try :
		f = open(modelfile,"r")
		for line in f : 
#			print line
			if(line[0]=='#') :      # Comment lines
				continue	
			data = line.split()
			Z = data[0]
			age = data[1]
			v = data[10]
			b = data[9]
			
			fnd=0
			for c in CMDs : 
				if(c.age == age and c.Z == Z) :
					c.add_point(b,v)
					fnd=1

			if(not fnd) :
				CMDs.append(CMDmodel(age,Z))
#				print "added new CMD"
				CMDs[-1].add_point(b,v)
#				print "added new point"
		
		f.close()
		
	except : 
			print ("Error loading model CMDs from file " + modelfile) 
			return 
			
	return CMDs

def loaddata(datafile='CMD.dat') : 
	''' 
	Loads a data file with format B-V, V
	''' 
	cmd = CMDdata(datafile)
	
	try : 
		f = open(datafile,"r")
		for line in f : 
			if(line[0]=='#') : 	# Comment lines
				continue
			data = line.split()
			# Modified to read in RA, DEC, B, V format produced for AS32 scripts.
			cmd.add_point(data[2],data[3])
	except : 
		print ("Error loading data CMD from file " + datafile)
		return
		
	return cmd
		
def chisq(cmd,model,dist,ext) :
	
	chisq=0.0
	n=0
	
	for ov,obv in np.array(zip(cmd.v,cmd.bv)) : 
		for mv,mbv in np.array(zip(model.v + dist + model.R*ext,model.bv + ext)) :		
			chisq += ((ov-mv)**2 + (mv-mbv)**2)
			n+=1
	
	return chisq/float(n)
	

# def autofit (cmd,models,modid,dist,ext) : 
	# '''
	# Autofit the best choice model around the current preference
	
	# cmd = data to fit
	# models = array of models
	# modid = mid-point model (proxy for age)
	# dist = mid-point distance
	# ext = mid-point extinction
	# '''
	
	# distrange=3.0
	# extrange=0.2
	# modrange=6
	
	# res=[]
	
	
	# for d in np.arange(dist-distrange/2.0,dist+distrange/2.0,distrange/5.0) : 
		# for e in np.arange(max(0.0,ext-extrange/2.0),ext+extrange/2.0,extrange/5.0) :
			# for m in np.arange(max(0,modid-modrange/2),min(len(models)-1,modid+modrange/2),1,dtype='int') : 
				# chi2 = chisq(cmd,models[m],d,e)
				# print d,e,m,chi2
				# res.append([models[m].age,d,e,chi2])
	# return res
	
def isofit (cmd=None,models=None) : 
	
	if models is None :
		print ('loading default model files')
		models = loadmodels()
	
	if cmd is None: 
		print ('loading default CMD data')
		cmd = loaddata()
		
	
	#plot the data file
	#get a fixed fig reference
	cmd.plot()			# Plot the CMD file
	fig = plt.gcf()		# Get the figure for use later
	ax = plt.gca()		# Get the axes for use later

	ax.set_autoscale_on(False)
	plt.ion()
	plt.show()
	fig.canvas.draw()
	plt.pause(0.001)
	
	# Fire into curses. Yeah, 1980s here we come :)
	
	stdscr = curses.initscr()
	curses.cbreak()
	curses.noecho()
	stdscr.keypad(1)
	stdscr.nodelay(1)
#	stdscr.timeout(0)
	
	stdscr.clear()
	stdscr.addstr(0,10,"ISOFIT")
	stdscr.addstr(3,10,"'q' to quit")
	stdscr.addstr(4,10,"'p' to print. Output file will be " + cmd.name + "_HHMMSS.pdf")
	stdscr.addstr(6,10,"'o'/'y' to increase/decrease model age")
	stdscr.addstr(7,10,"'+'/'-' to increase/decrease model distance")
	stdscr.addstr(8,10,"'9'/'0' to increase/decrease model reddening")

	stdscr.refresh()
	
	key = ''
	dist=10.0
	ext=0.1
	modid=int(len(models)/2)
	
	age=models[modid].age
	
	#go into loop with keyboard
	#   get key
	#   change parameters
	#   remove old plot
	#   update new plot
	#   quit on 'q'
	# 
	
	mplot=None
	i=0
	pulse='.+0+'
	

	lab = '10^{:.3s}'.format(models[modid].age) + "yr " + '{:04.1f}'.format(dist) + "mag " + '{:05.3f}'.format(ext) + " E(B-V)"
	mplot = ax.plot(np.array(models[modid].bv) + ext,np.array(models[modid].v) + dist + models[modid].R * ext,color='r',lw=3,label=lab)
	ax.legend(loc='best')
	fig.canvas.draw()
	plt.pause(0.001)
	
	while True :
		key = stdscr.getch()
		plt.pause(0.005)
		
		if key == -1 :
#			stdscr.addstr(10,15,pulse[i])
#			i = i+1
#			if i==len(pulse)  :
#				i=0
#			plt.pause(0.1)
#			stdscr.refresh()
			continue

		# Wipe the previous print statement...
		stdscr.addstr(12,10,"                                                      ")

		# Process the input			
		if key == ord('q') : 
			break;		
		if key == ord('+') or key==ord('=') :
			dist=min(15.0,dist+0.1)
		elif key == ord('-') or key==ord('_') : 
			dist=max(3.0,dist-0.1)
		elif key == ord('y') : 
			modid=max(0,modid-1)
		elif key == ord('o') : 
			modid=min(len(models)-1,modid+1)
		elif key == ord('9') : 
			ext=min(3.0,ext+0.025)
		elif key == ord('0') : 
			ext=max(0.0,ext-0.025)
		elif key == ord('p') : 
			ltime = time.localtime()
			fname = cmd.name + "_" + '{:02d}'.format(ltime.tm_hour) \
								   + '{:02d}'.format(ltime.tm_min) \
								   + '{:02d}'.format(ltime.tm_sec) \
								   + ".pdf"
			stdscr.addstr(12,10,"Saved current figure to " + fname)
			fig.savefig(fname)

		stdscr.addstr(10,10,"                                                           ")
		stdscr.addstr(10,10,"Age="+models[modid].age+" Distance="+str(dist)+" Extinction="+str(ext))
		stdscr.refresh()
					
		if mplot is not None : 
			mplot.pop(0).remove()
		
		lab = '10^{:.3s}'.format(models[modid].age) + "yr " + '{:04.1f}'.format(dist) + "mag " + '{:05.3f}'.format(ext) + " E(B-V)"
		
		mplot = ax.plot(np.array(models[modid].bv) + ext,np.array(models[modid].v) + dist + models[modid].R * ext,color='r',lw=3,label=lab)
		ax.legend(loc='best')
		fig.canvas.draw()
		plt.pause(0.001)
		
		curses.flushinp()	# Kill other characters so we don't get a lag...
				
	curses.endwin()
	
	return 
	
def version() :
	return "2020/04/28"

def main(argv=None) :
	if argv is None :
		argv = sys.argv
	
	cmdfile=None
	modelfile=None
	
	print (len(argv), argv)

	if(len(argv)>1) : 
		cmdfile=argv[1]
	if(len(argv)>2) : 
		modelfile=argv[2]

	cmd = loaddata(cmdfile)
	models = loadmodels(modelfile)
	
	isofit(cmd,models)
	
	return

if __name__ == "__main__" : 
	main()
