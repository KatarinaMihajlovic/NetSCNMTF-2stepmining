# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
#import seaborn as sns


#Load network
def Load_Network(fname):
	net = nx.read_edgelist(fname)
	nodes = [n for n in net.nodes()]
	nodeset = set(nodes)
	nb_nodes = len(nodes)
	nodes2ind = {}
	for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n
	return net, nodes, nodes2ind



#Load a symmetric x.x matrix and create the ordering of x-nodes
def Load_X_to_X(fname, esc):
	edges=[]
	x_nodes = set()
	ifile=open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split(esc)
		if len(lspt)>1:
			x1 = lspt[0]
			x2 = lspt[1]
			edges.append([x1,x2])
			x_nodes.add(x1)
			x_nodes.add(x2)
	
	x_nodes = list(x_nodes)
	
	nb_xnodes = len(x_nodes)
	
	print(" - %s: %i nodes, %i interactions"%(fname, nb_xnodes, len(edges)))
	
	x_nodes2ind = {}
	for n in range(nb_xnodes):
		x_nodes2ind[x_nodes[n]]=n
	
	X2X = np.zeros((nb_xnodes, nb_xnodes))
	for x1,x2 in edges:
		indx1 = x_nodes2ind[x1]
		indx2 = x_nodes2ind[x2]
		X2X[indx1][indx2] = 1.
		X2X[indx2][indx1] = 1.
	
	return X2X, x_nodes, x_nodes2ind



#Load a symmetric a.a matrix, where the ordering of a-nodes is given
def Load_A_to_A(fname, a_nodes2ind, esc):
	edges=[]
	ifile=open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split(esc)
		if len(lspt)>1:
			a1 = lspt[0]
			a2 = lspt[1]
			if a1 in a_nodes2ind and a2 in a_nodes2ind:
				edges.append([a1,a2])
	
	nb_anodes = len(a_nodes2ind)
	
	print(" - %s: %i nodes, %i interactions"%(fname, nb_anodes, len(edges)))
	
	A2A = np.zeros((nb_anodes, nb_anodes))
	for a1,a2 in edges:
		inda1 = a_nodes2ind[a1]
		inda2 = a_nodes2ind[a2]
		A2A[inda1][inda2] = 1.
		A2A[inda2][inda1] = 1.
	
	return A2A




#Load a x.a matrix, where the ordering of a-nodes is given, and create the ordering of x-nodes
def Load_X_to_A(fname, a_nodes2ind, esc):
	edges=[]
	a_nodes = set()
	x_nodes = set()
	ifile=open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split(esc)
		if len(lspt)>1:
			x = lspt[0]
			a = lspt[1]
			if a in a_nodes2ind:
				edges.append([x,a])
				x_nodes.add(x)
				a_nodes.add(a)
	
	x_nodes = list(x_nodes)
	
	nb_xnodes = len(x_nodes)
	nb_anodes = len(a_nodes2ind)
	
	print(" - %s: %i - %i nodes, %i interactions"%(fname, nb_xnodes, nb_anodes, len(edges)))
	
	x_nodes2ind = {}
	for n in range(nb_xnodes):
		x_nodes2ind[x_nodes[n]]=n
	
	X2A = np.zeros((nb_xnodes, nb_anodes))
	for x,a in edges:
		indx = x_nodes2ind[x]
		inda = a_nodes2ind[a]
		X2A[indx][inda] = 1.
	
	return X2A, x_nodes, x_nodes2ind




def Load_GO(fname, allowed_nodes):
	node2go=[]
	nodes = set()
	gos = set()
	ifile=open(fname)
	for line in ifile.readlines():
		lspt = line.strip().split(';')
		if len(lspt)>1:
			g = lspt[0]
			t = lspt[1]
			if g in allowed_nodes:
				node2go.append([g,t])
				nodes.add(g)
				gos.add(t)
	
	nodes = list(nodes)
	gos = list(gos)
	
	nb_nodes = len(nodes)
	nodes2ind = {}
	for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n
		
	nb_gos = len(gos)
	gos2ind = {}
	for n in range(nb_gos):
		gos2ind[gos[n]]=n
	
	gene2go = np.zeros((nb_nodes, nb_gos))
	for g,t in node2go:
		indg = nodes2ind[g]
		indt = gos2ind[t]
		gene2go[indg][indt] = 1.
	
	return gene2go, nodes, nodes2ind, gos, gos2ind




#Computing laplacian matrix, L
def Make_Laplacian(net, net_nodes, net_n2i):
	print("Computing laplacian matrix")
	nb_nodes = len(net_nodes)
	A = np.zeros((nb_nodes,nb_nodes))
	D = np.zeros((nb_nodes,nb_nodes))
	for n1 in range(nb_nodes):
		D[n1][n1] = float(net.degree(net_nodes[n1]))
	
	for e in net.edges():
		n1=net_n2i[ e[0] ]
		n2=net_n2i[ e[1] ]
		A[n1][n2] = 1.
		A[n2][n1] = 1.
	
	L = D-A
	return L



#Computing Adjacency matrix, A
def Make_Adj(net, net_nodes, net_n2i, MOL_or_EXP):
	if MOL_or_EXP == 'MOL':
		nb_nodes = len(net_nodes)
		A = np.zeros((nb_nodes,nb_nodes))
		
		for e in net.edges():
			if e[0] in net_n2i and e[1] in net_n2i:
				n1=net_n2i[ e[0] ]
				n2=net_n2i[ e[1] ]
				A[n1][n2] = 1.
				A[n2][n1] = 1.
		return A

	elif MOL_or_EXP == 'EXP':
		gene_nodes = net_nodes
		SC_nodes = list(net)
		i2SC = dict(list(enumerate(list(net))))
		A = np.zeros((len(gene_nodes),len(SC_nodes)))
		
		for gene in range(len(gene_nodes)):
			for SC in range(len(SC_nodes)):
				A[gene][SC] = net[SC_nodes[SC]][gene_nodes[gene]]
		return A, i2SC




#Compute PPMI matrix from adjacency matrix A
def Make_PPMI(A, context_window=10, Bin=False):
	print("Computing PPMI matrix")
	degrees = np.sum(A,axis = 0)
	volume_of_graph = sum(degrees)
	
	#fix for disconnected nodes
	degree_corrected = np.array( [float(i) for i in degrees] )
	for i in range(len(degree_corrected)):
		if degree_corrected[i] == 0.:
			degree_corrected[i] = 1.
	
	diag_degrees_inv = np.diag([1./float(i) for i in degree_corrected])
	
	tpm = np.matmul(diag_degrees_inv,A)
	power_tpm = np.zeros([len(A),len(A)]) + tpm
	pmi = np.zeros([len(A),len(A)]) + power_tpm
	
	for r in range(2,context_window+1):
		#print 'Iteration: ',r
		power_tpm = np.matmul(power_tpm, tpm)
		pmi += power_tpm
	
	pmi = pmi/float(context_window)
	pmi = volume_of_graph*np.matmul(pmi,diag_degrees_inv)
	pmi_matrix = np.log(pmi, out=np.zeros_like(pmi), where=(pmi!=0))
	
	for i in pmi_matrix:
		i[i<0]=0.
	
	if Bin==True:
		for i in pmi_matrix:
			i[i>0.]=1. 
	
	return pmi_matrix


def Sub_Matrix(P, nodes, new_node2ind):
	NP = np.zeros((len(new_node2ind), P.shape[1]))
	for i in range(len(nodes)):
		if nodes[i] in new_node2ind:
			ni = new_node2ind[nodes[i]]
			for j in range(P.shape[1]):
				NP[ni][j] = P[i][j]
	return NP



def Load_SemSims(fname, go2ind):
	SS = np.zeros((len(go2ind), len(go2ind)))
	ifile = open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		if len(lspt)> 2:
			go1 = int(lspt[0])
			go2 = int(lspt[1])
			ss12 = float(lspt[2])
			gname1 = "GO:%07i"%(go1)
			gname2 = "GO:%07i"%(go2)
			id1 = go2ind[gname1]
			id2 = go2ind[gname2]
			SS[id1][id2] = ss12
			SS[id2][id1] = ss12
	ifile.close()
	return SS

#Compute the PPMI matrix of a bipartite adjacency matrix between x and y nodes.
def Make_PPMI_Bipartite(Mat, x2ind, y2ind, context_window=10, Bin=False):
	nb_x, nb_y = Mat.shape
	Full_Adj = np.zeros( (nb_x+nb_y, nb_x+nb_y) )
	for i in range(nb_x):
		#Full_Adj[i][i] = 1.
		for j in range(nb_y):
			Full_Adj[i][nb_x+j] = Mat[i][j]
			Full_Adj[nb_x+j][i] = Mat[i][j]
	
	#for j in range(nb_y):
	#	Full_Adj[nb_x+j][nb_x+j] = 1.
	
	Full_PPMI = Make_PPMI(Full_Adj, context_window, Bin)
	Sub_PPMI = np.zeros((nb_x, nb_y))
	for i in range(nb_x):
		for j in range(nb_y):
			Sub_PPMI[i][j] = Full_PPMI[i][nb_x+j]
	return Sub_PPMI



def Save_Square_Matrix(matrix, rownames, fname):
	ofile = open(fname, 'w')
	n,m = matrix.shape
	for i in range(n):
		ni = rownames[i]
		for j in range(i+1,n):
			if matrix[i][j] > 0.:
				nj = rownames[j]
				ofile.write("%s\t%s\t%s\n"%(ni, nj, str(matrix[i][j])))
	ofile.close()



def Save_Rectangular_Matrix(matrix, rownames, colnames, fname):
	ofile = open(fname, 'w')
	n,m = matrix.shape
	for i in range(n):
		ni = rownames[i]
		for j in range(m):
			if matrix[i][j] > 0.:
				nj = colnames[j]
				ofile.write("%s\t%s\t%s\n"%(ni, nj, str(matrix[i][j])))
	ofile.close()


def Load_Square_Matrix(node2ind, fname, Bin=False):
	n = len(node2ind)
	M = np.zeros((n,n))
	ifile = open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		if len(lspt)> 2:
			n1 = node2ind[lspt[0]]
			n2 = node2ind[lspt[1]]
			val = float(lspt[2])
			if Bin==True:
				val = 1.
			M[n1][n2] = val
			M[n2][n1] = val
	ifile.close()
	return M


def Load_Rectangular_Matrix(row2ind, col2ind, fname, Bin=False):
	n = len(row2ind)
	m = len(col2ind)
	M = np.zeros((n,m))
	ifile = open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		if len(lspt)> 2:
			n1 = row2ind[lspt[0]]
			n2 = col2ind[lspt[1]]
			val = float(lspt[2])
			if Bin==True:
				val = 1.
			M[n1][n2] = val
	ifile.close()
	return M


def Plot_Mat(matrix, labels, fname):
	plt.matshow(matrix,cmap="bwr")
	
	
	#ax = sns.heatmap(matrix, linewidth=0., xticklabels=20, yticklabels=20)
	
	#plt.ylabel("Patients", fontsize=20)
	#plt.xlabel("Clusters", fontsize=20)
	plt.xticks(labels[0], labels[1])
	plt.yticks(labels[0], labels[1])
	plt.colorbar()
	#plt.tick_params(axis='both', which='major', labelsize=18)
	plt.savefig(fname, dpi=600)
	plt.clf()


def Save_Clustering(Clusters, fname):
	ofile = open(fname, 'w')
	for i in range(len(Clusters)):
		for j in range(len(Clusters[i])):
			item = Clusters[i][j]
			ofile.write("%s\t%i\n"%(item, i))
	ofile.close()



def Save_Matrix_Factor(M, fname, rownames, colnames):
	n,m = M.shape
	ofile = open(fname, 'w')
	#column header
	for i in range(m):
		ofile.write("\t%s"%(colnames[i]))
	ofile.write("\n")
	#rows
	for i in range(n):
		ofile.write("%s"%(rownames[i]))
		for j in range(m):
			ofile.write("\t%s"%(str(M[i][j])))
		ofile.write("\n")
	ofile.close()

def Save_Matrix_Factor_no_headers(M, fname, rownames, colnames):
	n,m = M.shape
	ofile = open(fname, 'w')
	#column header
	# for i in range(m):
	# 	ofile.write("\t%s"%(colnames[i]))
	# ofile.write("\n")
	#rows
	for i in range(n):
		# ofile.write("%s"%(rownames[i]))
		for j in range(m):
			if j > 0: ofile.write("\t")
			ofile.write("%s"%(str(M[i][j])))
		ofile.write("\n")
	ofile.close()
	
	
def Load_Matrix_Factor(fname):
	ifile = open(fname, 'r')
	#column header
	line = ifile.readline()
	colnames = line.strip().split('\t')
	rownames = []
	matrix = []
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		rownames.append(lspt[0])
		matrix.append( [float(i) for i in lspt[1:]])
	ifile.close()
	return matrix, rownames, colnames


