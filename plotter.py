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


	def __init__(self, d, path=""):
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
		self.minxlim = -0.01
		self.maxxlim = 1.01
		self.minylim = -0.01
		self.maxylim = 0.5
	
	def setLims(self, maxX, maxY):
		self.maxxlim = maxX+0.01
		self.maxylim = 3*maxY+0.01
		
	def plotDelaunay(self, xp, yp, triangles):
		for t in triangles:
			# t[0], t[1], t[2] are the points indexes of the triangle
			t_i = [t[0], t[1], t[2], t[0]]
			plt.plot(xp[t_i],yp[t_i])

		plt.plot(xp,yp,'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		#plt.draw()
		plt.savefig(self.path + "step"+str(self.i)+".png")
		plt.clf()
		
	def next(self):
		self.i += 1
	
	def plotVoronoi(self, cs, xp, yp, edges, circumcenters, neighbors, triangles):
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		vertices = [[] for _ in xrange(len(xp))] #array with number of points
		for cc, t in zip(circumcenters, triangles):
			'''add this circumcenter as a vertex of each of
			the points which make up the triangle'''
			vertices[t[0]].append(cc)
			vertices[t[1]].append(cc)
			vertices[t[2]].append(cc)
		for i,v in enumerate(vertices):
			v = np.array(v)
			idx = self.ccworder(v.T)
			v = v[idx]
			vertices[i] = v
		fc = []
		for i in xrange(cs.size):
			if cs.state[i] == 1:
				fc.append('g')
			else:
				fc.append('r')
		coll = PolyCollection(vertices, facecolors=fc, edgecolors='k', alpha = 0.5)
		ax.add_collection(coll)
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "vor"+str(self.i).zfill(3)+".png")
		plt.clf()
	
	def ccworder(self, A):
		A= A- np.mean(A, 1)[:, None]
		return np.argsort(np.arctan2(A[1, :], A[0, :]))
		
	def plotCellStates(self, cs, xp, yp, edges, circumcenters, neighbors):
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		#xp, yp is same as c.x, c.y
		for i in xrange(cs.size):
			if cs.state[i] == 1:
				co = 'g'
			else:
				co = 'r'
			circle = plt.Circle((cs.x[i, 3], cs.y[i, 3]), radius=self.dh, color=co, alpha = 0.5)
			ax.add_patch(circle)
			ax.text(cs.x[i, 3]+0.01, cs.y[i, 3]+0.01, str(i))
			
#		for i in xrange(len(neighbors)):
#			j = neighbors[i][0]
#			k = neighbors[i][1]
#			l = neighbors[i][2]
#			if j > -1:
#				plt.plot([circumcenters[i][0], circumcenters[j][0]], [circumcenters[i][1], circumcenters[j][1]], color='k')
#			if k > -1:
#				plt.plot([circumcenters[i][0], circumcenters[k][0]], [circumcenters[i][1], circumcenters[k][1]], color='k')
#			if l > -1:
#				plt.plot([circumcenters[i][0], circumcenters[l][0]], [circumcenters[i][1], circumcenters[l][1]], color='k')
		
		for e in edges:
			plt.plot(xp[e],yp[e], color = 'b')
			
		plt.plot(xp,yp,'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "state"+str(self.i).zfill(3)+".png")
		plt.clf()
		
	def drawCells(self, cs, triangles, xp, yp):
		dic = self.mapPointsToPos(xp, yp)
		#fx, fy = self.createFakePoints(xp, yp)
		ee, op = self.findExternalEdges(triangles, xp, yp)
		fx, fy = self.createAuxPoints(ee, op, xp, yp)
		tri = triang.Triangulation(fx,fy)
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
			if cs.state[cidx] == 1:
				fstate.append('g')
			else:
				fstate.append('r')
			#array with pressure values
			fpress.append(cs.V[cidx])
		
		coll = PolyCollection(vertices, facecolors=fstate, edgecolors='k', alpha = 0.5)
		
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		ax.add_collection(coll)
		plt.plot(xp, yp, 'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "vor"+str(self.i).zfill(3)+".png")
		
		for e in tri.edge_db:
			plt.plot(tri.x[e],tri.y[e], color = 'b')
		for e in ee:
			ei = np.array(e)
			plt.plot(xp[ei],yp[ei], color = 'r')
		plt.xlim(self.minxlim - 0.2, self.maxxlim + 0.2)
		plt.ylim(self.minylim - 0.2, self.maxylim + 0.2)
		plt.savefig(self.path + "vorDebug"+str(self.i)+".png")
		plt.clf()
		
		fig = plt.figure(figsize=self.figsize)
		ax = fig.add_subplot(1, 1, 1)
		#jet = plt.get_cmap('jet') 
		cNorm  = colors.Normalize(vmin=0, vmax=0.01)
		coll = PolyCollection(vertices, array=np.array(fpress), cmap=cm.rainbow, norm=cNorm, edgecolors='k')
		ax.add_collection(coll)
		fig.colorbar(coll, ax=ax)
		plt.plot(xp, yp, 'o')
		plt.xlim(self.minxlim, self.maxxlim)
		plt.ylim(self.minylim, self.maxylim)
		plt.savefig(self.path + "press"+str(self.i).zfill(3)+".png")
		plt.clf()
		
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
		fx = np.copy(xp)
		fy = np.copy(yp)
		for e in ee:
			#index of third point in the triangle
			#which contains the external edge
			opi = op[e]
			dx = xp[e[0]] - xp[e[1]]
			dy = yp[e[0]] - yp[e[1]]
			if abs(dy) < 1e-8:
				xn = xp[opi]
				yn = 2*yp[e[0]] - yp[opi]
			elif abs(dx) < 1e-8:
				xn = 2*xp[e[0]] - xp[opi]
				yn = yp[opi]
			else:
				#y = ax + b defines the external edge
				a = dy/dx
				b = yp[e[1]] - xp[e[1]]*a
				#reflect third point along external edge
				d = (xp[opi] + (yp[opi] - b)*a)/(1 + a*a)
				xn = 2*d - xp[opi]
				yn = 2*d*a - yp[opi] + 2*b
			fx = np.append(fx, xn)
			fy = np.append(fy, yn)
		#fix bottom corners
		fx = np.append(fx, -self.dh)
		fy = np.append(fy, -self.dh)
		fx = np.append(fx, 1+self.dh)
		fy = np.append(fy, -self.dh)
		#fix top corners
		top = np.max(yp)
		fx = np.append(fx, -self.dh)
		fy = np.append(fy, top+self.dh)
		fx = np.append(fx, 1+self.dh)
		fy = np.append(fy, top+self.dh)
		return fx,fy
		
	def createFakePoints(self, xp, yp):
		self.xmax = np.max(xp)
		self.xmin = np.min(xp)
		self.ymax = np.max(yp)
		self.ymin = np.min(yp)
		right = self.xmax + self.d
		left = self.xmin - self.d
		top = self.ymax + self.d
		bottom = self.ymin - self.d
		xran = np.concatenate((np.arange(left,right,self.d), [right]))
		yran = np.concatenate((np.arange(bottom,top,self.d), [top]))
		#right
		fx = np.append(xp, xran)
		fy = np.append(yp, np.ones(len(xran))*top)
		#bottom
		fx = np.append(fx, xran)
		fy = np.append(fy, np.ones(len(xran))*bottom)
		#left
		fx = np.append(fx, np.ones(len(yran))*left)
		fy = np.append(fy, yran)
		#right
		fx = np.append(fx, np.ones(len(yran))*right)
		fy = np.append(fy, yran)
		return fx, fy