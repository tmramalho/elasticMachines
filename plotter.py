'''
Created on Feb 6, 2013

@author: tiago
'''
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.delaunay.triangulate as triang
import os
import itertools

class Plotter(object):
	'''
	classdocs
	'''


	def __init__(self, d, pt, path=""):
		'''
		Constructor
		'''
		if path == "" or not os.path.isdir(path):
			self.path = os.getcwd() + "/data/"
			if not os.path.isdir(self.path):
				os.mkdir(self.path)
		else:
			self.path = path
		self.i = 0
		#plt.ion()
		self.d = d
		self.dh = d/2
		self.md = d*2
		self.figsize = (8/(10*d),6)
		self.minxlim = -0.001
		self.maxxlim = 1.001
		self.minylim = -0.001
		self.maxylim = 0.501
		self.pressThresh = pt
	
	def calcDelaunay(self):
		#delaunay
		self.tt = triang.Triangulation(self.xp,self.yp)
		self.edges = [e for e in self.tt.edge_db if self.dist(e,self.xp,self.yp) < self.md]
		self.numEdges = len(self.edges)
	
	def setLims(self, maxX, maxY):
		self.maxxlim = maxX+0.001
		self.maxylim = maxY+0.001
		
	def plotDelaunay(self):
		plt.clf()
		plt.figure(figsize=self.figsize)
		for t in self.tt.triangle_nodes:
			# t[0], t[1], t[2] are the points indexes of the triangle
			t_i = [t[0], t[1], t[2], t[0]]
			plt.plot(self.xp[t_i],self.yp[t_i])

		plt.plot(self.xp,self.yp,'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "step"+str(self.i)+".png")
		
	def next(self):
		self.i += 1
		
	def setupCellSystem(self, cs):
		self.xp = np.copy(cs.x[:,3])
		self.yp = np.copy(cs.y[:,3])
		self.V = np.copy(cs.V)
		self.state = np.copy(cs.state)
		self.nCells = cs.size
		self.setup()
		
	def setupFromArrays(self, x, y, s, v):
		self.xp = x
		self.yp = y
		self.state = s
		self.V = v
		self.nCells = len(v)
		self.setup()
		
	def setup(self):
		self.calcDelaunay()
		dic = self.mapPointsToPos(self.xp, self.yp)
		self.extedges, op = self.findExternalEdges(self.tt.triangle_nodes, self.xp, self.yp)
		fx, fy = self.createAuxPoints(self.extedges, op, self.xp, self.yp)
		self.tri = triang.Triangulation(fx,fy)
		self.vertices, self.fstate, self.fpress = self.createVertexArrayFromTriangulation(self.tri, dic)
	
	def ccworder(self, A):
		A= A- np.mean(A, 1)[:, None]
		return np.argsort(np.arctan2(A[1, :], A[0, :]))
		
	def plotCellSystem(self):
		plt.clf()
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		#xp, yp is same as c.x, c.y
		for i in xrange(self.nCells):
			ptest = self.V[i] < self.pressThresh
			if self.state[i] == 1 and ptest:
				co = 'g'
			elif self.state[i] == 1:
				co = 'y'
			elif ptest:
				co = 'b'
			else:
				co = 'r'
			circle = plt.Circle((self.xp[i], self.yp[i]), radius=self.dh, color=co, alpha = 0.5)
			ax.add_patch(circle)
			ax.text(self.xp[i], self.yp[i], str(i))
			ax.text(self.xp[i], self.yp[i]-0.002, "%.2e" % self.V[i], fontsize=9, color='m')
			
		for e in self.edges:
			plt.plot(self.xp[e],self.yp[e], color = 'b')
			
		plt.plot(self.xp,self.yp,'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "system"+str(self.i).zfill(3)+".png")
		
	def drawCellStates(self, debug=False):
		plt.clf()
		coll = PolyCollection(self.vertices, facecolors=self.fstate, edgecolors='k', alpha = 0.5)
		
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		ax.add_collection(coll)
		plt.plot(self.xp, self.yp, 'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "state"+str(self.i).zfill(3)+".png")
		
		if debug:
			for e in self.tri.edge_db:
				plt.plot(self.tri.x[e],self.tri.y[e], color = 'b')
			for e in self.extedges:
				ei = np.array(e)
				plt.plot(self.xp[ei],self.yp[ei], color = 'r')
			plt.xlim(1.1*self.minxlim, 1.1*self.maxxlim)
			plt.ylim(1.1*self.minylim, 1.1*self.maxylim)
			plt.savefig(self.path + "vorDebug"+str(self.i)+".png")
		
	def drawCellPressure(self):
		plt.clf()
		pcols = []
		for fp in self.fpress:
			if fp < self.pressThresh:
				pcols.append('b')
			else:
				pcols.append('c')
		coll = PolyCollection(self.vertices, facecolors=pcols, edgecolors='k', alpha = 0.5)
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		ax.add_collection(coll)
		plt.plot(self.xp, self.yp, 'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "pstate"+str(self.i).zfill(3)+".png")
		
	def drawPressureValues(self):
		plt.clf()
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		#jet = plt.get_cmap('jet') 
		cNorm  = colors.Normalize(vmin=0, vmax=self.pressThresh*2)
		coll = PolyCollection(self.vertices, array=np.array(self.fpress), cmap=cm.rainbow, norm=cNorm, edgecolors='k')
		ax.add_collection(coll)
		fig.colorbar(coll, ax=ax)
		plt.plot(self.xp, self.yp, 'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "press"+str(self.i).zfill(3)+".png")
		
	def createVertexArrayFromTriangulation(self, tri, dic):
		auxVerts = [[] for _ in xrange(len(tri.x))] #array with number of points
		for cc, t in zip(tri.circumcenters, tri.triangle_nodes):
			'''add this circumcenter as a vertex of each of
			the points which make up the triangle'''
			auxVerts[t[0]].append(cc)
			auxVerts[t[1]].append(cc)
			auxVerts[t[2]].append(cc)
		vertices = []
		fstate = []
		fpress = []
		for i,v in enumerate(auxVerts):
			try:
				cidx = dic[(tri.x[i], tri.y[i])]
			except KeyError:
				#belongs to fake points
				continue
			#sort vertices and append to array
			v = np.array(v)
			vidx = self.ccworder(v.T)
			v = v[vidx]
			vertices.append(v)
			#create array with state colors
			if self.state[cidx] == 1:
				fstate.append('g')
			else:
				fstate.append('r')
			#array with pressure values
			fpress.append(self.V[cidx])
		return vertices, fstate, fpress
		
	def mapPointsToPos(self, xp, yp):
		dic = {}
		for i in xrange(len(xp)):
			dic[(xp[i],yp[i])] = i
		return dic
	
	def dist(self, e, xp, yp):
		dx = xp[e[0]] - xp[e[1]]
		dy = yp[e[0]] - yp[e[1]]
		return np.sqrt(dx*dx+dy*dy)
	
	def findExternalEdges(self, tri, xp, yp):
		edges = set()
		opposites = dict()
		for t in tri:
			if np.any(t < 0):
				continue
			pts = list(itertools.combinations(t, 2))
			badTri = False
			#check if this triangle has a removed edge
			#if yes ignore it
			for p in pts:
				d = self.dist(p, xp, yp)
				if d > self.md:
					badTri = True
					break
			if badTri:
				continue
			for p in pts:
				ps = tuple(sorted(p))
				op = [x for x in t if x not in ps]
				if ps in edges:
					edges.remove(ps)
					del opposites[ps]
				else:
					edges.add(ps)
					opposites[ps] = op[0]
		return edges, opposites
	
	def createAuxPoints(self, ee, op, xp, yp):
		'''
		This method by Miguel Miranda 2013
		'''
		fx = np.copy(xp)
		fy = np.copy(yp)
		for e in ee:
			#index of third point in the triangle
			#which contains the external edge
			opi = op[e]
			dx = xp[e[0]] - xp[e[1]]
			dy = yp[e[0]] - yp[e[1]]
			daux = dx*dx+dy*dy
			dw = np.round(np.sqrt(daux-daux/4),4)
			if abs(dy) < 1e-8:
				wn = -1 if yp[opi] > yp[e[0]] else 1
				xn = xp[e[0]] - dx/2
				yn = yp[e[0]] + wn * dw
			elif abs(dx) < 1e-8:
				wn = -1 if xp[opi] > xp[e[0]] else 1
				xn = xp[e[0]] + wn * dw
				yn = yp[e[0]] - dy/2
			else:
				rx = dy * 0.7
				ry = -dx * 0.7
				
				ox = xp[opi] - xp[e[0]]
				oy = yp[opi] - yp[e[0]]

				zxc = (ox+rx)*(ox+rx)+(oy+ry)*(oy+ry)
				asd = (ox-rx)*(ox-rx)+(oy-ry)*(oy-ry)
				
				wn = 1 if zxc < asd else -1
								
				xn = (xp[e[0]] - dx/2) + wn * rx
				yn = (yp[e[0]] - dy/2) + wn * ry
				
			fx = np.append(fx, xn)
			fy = np.append(fy, yn)
		#maxs
		bottom = np.min(fy)
		top = np.max(fy)
		left = np.min(fx)
		right = np.max(fx)
		#fix bottom corners
		fx = np.append(fx, left)
		fy = np.append(fy, bottom)
		fx = np.append(fx, right)
		fy = np.append(fy, bottom)
		#fix top corners
		fx = np.append(fx, left)
		fy = np.append(fy, top)
		fx = np.append(fx, right)
		fy = np.append(fy, top)
		return fx,fy
	