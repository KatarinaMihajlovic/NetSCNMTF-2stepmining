# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import numpy as np
import scipy.sparse.linalg as ssl
from sklearn.decomposition import TruncatedSVD
import seaborn as sn

epsilon = 1e-5

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()



def Get_Soft_Clusters_SM(M, nodes, thr=0.1):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	
	for i in range(n):
		norm_row = softmax(M[i])
		print(max(norm_row))
		for j in range(k):
			if norm_row[j] >= thr:
				Clusters[j].append(nodes[i])
	return Clusters




def Get_Soft_Clusters(M, nodes, thr=0.1):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	
	for i in range(n):
		norm_row = sum(M[i])
		for j in range(k):
			if M[i][j]/norm_row >= thr:
				Clusters[j].append(nodes[i])
	return Clusters


def Get_Hard_Clusters(M, nodes):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	for i in range(n):
		idx = np.argmax(M[i])
		Clusters[idx].append(nodes[i])
	return Clusters



def Get_Hard_Clusters_prime(M, nodes, selected_nodes):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	for i in range(n):
		if nodes[i] in selected_nodes:
			idx = np.argmax(M[i])
			Clusters[idx].append(nodes[i])
	return Clusters



	

class PSBLike_Drug_NMTF:
	#
	#  Min  || DTI - G.S.D^T ||^2_F + gamma D^T.Lchem.D
	#  subject to D >=0
	
	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		
	
	def Score(self, DTI, G, S, D, LChem, gamma, norm_DTI):
		DT = np.transpose(D)
		Err = DTI - np.matmul(G,  np.matmul(S, DT))
		norm_Err = np.linalg.norm(Err, ord='fro')**2
		rel = norm_Err / norm_DTI
		Err_Cst = gamma*np.trace(np.matmul(DT, np.matmul(LChem, D)))
				
		return norm_Err+Err_Cst, rel
	
	
	
	def Solve_MUR(self, DTI, Lchem, G, k1, k2, gamma, init="rand"):
		n,m = np.shape(DTI)
		
		
		norm_DTI = np.linalg.norm(DTI, ord='fro')
		
		if init == "rand":
			print(" * Using random initialization")
			D = np.random.rand(m, k2)
		elif init == 'SVD':
			print(" * Using SVD based initialization")
			
			D = np.zeros((m, k2)) + epsilon
			
			J,K,L = np.linalg.svd(DTI)
			for i in range(m):
				for j in range(k2):
					if L[j][i] > 0.:
						D[i][j] = L[j][i]
			
		
		else :
			print("Unknown initializer: %s"%(init))
			exit(0)
		
		# Computing init S
		GT = np.transpose(G)
		GTG = np.matmul( GT, G)
		GTG_inv = np.linalg.inv(GTG)
		GTG_inv_GT = np.matmul(GTG_inv, GT)
		
		DT = np.transpose(D)
		DTD = np.matmul( DT, D)
		DTD_inv = np.linalg.inv(DTD)
		
		S = np.matmul(GTG_inv_GT, np.matmul(DTI, np.matmul(D, DTD_inv)))
		OBJ, REL = self.Score(DTI, G, S, D, Lchem, gamma, norm_DTI)
		print(" - Init:\t OBJ:%.4f\t REL:%.4f"%(OBJ, REL))
		
		
		#Begining M.U.R.
		DTIT = np.transpose(DTI)
		
		for it in range(1,self.max_iter+1):
			
			# Update rule for D
			DTIT_G_S = np.matmul(DTIT, np.matmul(G, S))
			DTIT_G_S_pos = (np.absolute(DTIT_G_S) + DTIT_G_S)/2.0
			DTIT_G_S_neg = (np.absolute(DTIT_G_S) - DTIT_G_S)/2.0    
			
			ST_GTG_S = np.matmul(np.transpose(S), np.matmul(GTG, S))
			D_ST_GTG_S_pos = np.matmul(D, ((np.absolute(ST_GTG_S) + ST_GTG_S)/2.0))
			D_ST_GTG_S_neg = np.matmul(D, ((np.absolute(ST_GTG_S) - ST_GTG_S)/2.0))
			
			DNum = DTIT_G_S_pos + D_ST_GTG_S_neg
			DDenom = DTIT_G_S_neg + D_ST_GTG_S_pos
			
			# Contrain on G2, contribution from G2T.L2.G2
			DConst_pos = np.matmul(((np.absolute(Lchem) + Lchem)/2.0), D)
			DConst_neg = np.matmul(((np.absolute(Lchem) - Lchem)/2.0), D)     
			
			
			#Check that we are note dividing by zero  
			Dupdate = np.sqrt(np.divide(DNum + gamma*DConst_neg, DDenom + gamma*DConst_pos + epsilon))
			D = np.multiply(D, Dupdate) + epsilon
			
			# Computing S using closed formula
			DT = np.transpose(D)
			DTD = np.matmul( DT, D)
			DTD_inv = np.linalg.inv(DTD)
			
			S = np.matmul(GTG_inv_GT, np.matmul(DTI, np.matmul(D, DTD_inv)))
			
			if (it%self.verbose == 0) or (it==1):
				OBJ, REL = self.Score(DTI, G, S, D, Lchem, gamma, norm_DTI)
				print(" - It %i:\t OBJ:%.4f\t REL:%.4f"%(it, OBJ, REL))
		return D, S



class PD_SSNMTF:
	#
	#  Mol[i] = G.S[i].G^T
	#  subject to G >=0

	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		

	def Score(self, MOLs, EXP, G1, G2, S_Mol, S_EXP, norm_Mol, norm_EXP):
		Square_Error = 0.    
		norm_Err_MOLs = []
		rel_MOLs = []
		G1T = np.transpose(G1)
		for i in range(len(MOLs)):
			if norm_Mol[i] != 0:
				Err_MOLs = MOLs[i] - np.matmul(G1,  np.matmul(S_Mol[i], G1T))
				norm_Err_MOLs.append( np.linalg.norm(Err_MOLs, ord='fro'))
				rel_MOLs.append( norm_Err_MOLs[i] / norm_Mol[i] )
				Square_Error += norm_Err_MOLs[i]**2

		rel_EXP = []
		G2T = np.transpose(G2)
		Err_EXP = EXP - np.matmul(G1,  np.matmul(S_EXP, G2T))
		norm_Err_EXP = np.linalg.norm(Err_EXP, ord='fro')
		rel_EXP.append(norm_Err_EXP / norm_EXP)

		Square_Error += norm_Err_EXP**2
		rel_MOLs_EXP = rel_MOLs + rel_EXP

		return Square_Error, rel_MOLs_EXP

	
		
	
	def Solve_MUR(self, MOLs, EXP, k1, k2, init="rand"):
		import matplotlib.pyplot as plt
		import time

		g,_ = np.shape(MOLs[0])
		c = np.shape(EXP[0])[0]
		nb_mols = len(MOLs)

		norm_MOLs = [np.linalg.norm(MOLs[i], ord='fro') for i in range(nb_mols)]
		norm_EXP = np.linalg.norm(EXP, ord='fro')

		#--------------------------------------------------------
		#                     
		# Initialize G1 and G2
		#
		#-------------------------------------------------------- 

		if init == "rand":
			print(" * Using random initialization")
			G1 = np.random.rand(g, k1)
			G2 = np.random.rand(c, k2)
			#S = [np.random.rand(k_g,k_g) for i in range(nb_mols)]

		elif init == 'SVD':
			print(" * Using SVD based initialization")

			G1 = np.zeros((g, k1)) + epsilon
			G2 = np.zeros((c, k2)) + epsilon
			#S = [np.zeros((k_g,k_g)) + epsilon for i in range(nb_mols)]

			nbf = float(nb_mols)
			SMOL = np.zeros((g,g))
			for i in range(nb_mols):
				SMOL += MOLs[i]


			if np.linalg.norm(SMOL, ord='fro') != 0:
				print(" * -- Eig decomposition on MOL to get G1.Si.G1^T")
				#print(np.linalg.norm(SMOL, ord='fro'), SMOL.ndim, SMOL.shape, type(SMOL))
				M, N, O = ssl.svds(SMOL, k1, solver='lobpcg')
				M = abs(M)
				#print(M)
				for i in range(g):
					for j in range(k1):
						G1[i][j] += M[i][j]

			
			SEXP = EXP
			# A = J * K * L
			if k1 > k2:            
				J, K, L = ssl.svds(SEXP, k1, solver='lobpcg')
			else:
				J, K, L = ssl.svds(SEXP, k2, solver='lobpcg')
			
			J = abs(J)
			L = abs(L)
			#print(J)
			#print(L)
			
			print(" * -- Eig decomposition on EXP to get G1.S_EXP.G2^T")
			for i in range(c):
				for j in range(k2):
					G2[i][j] += L[j][i]

			for i in range(g):
				for j in range(k1):
					G1[i][j] += J[i][j]


		else :
			print("Unknown initializer: %s"%(init))
			exit(0)
		



        #--------------------------------------------------------
        #                     
        # Initialize S_Mol and S_EXP
        #
        #-------------------------------------------------------- 

		#computing Ss for Molecular matrices
		#  S_Mol[i] <- (G1T.G1)^-1 G1T.MOLs[i].G1 (G1T.G1)^-1
		
		S_Mol = []
		G1T = np.transpose(G1)
		G1TG1 = np.matmul( G1T, G1)
		#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
		G1TG1_inv = np.linalg.inv(G1TG1)
		G1_G1TG1_inv = np.matmul(G1, G1TG1_inv)
		G1TG1_inv_G1T = np.matmul(G1TG1_inv, G1T)
		for i in range(nb_mols):
		    S_Mol.append( np.matmul(G1TG1_inv_G1T, np.matmul(MOLs[i], G1_G1TG1_inv)) )


		# computing S for Expression matrices
		#  S_EXP <- (G1T.G1)^-1 G1T.EXP.G2 (G2T.G2)^-1
		
		G2T = np.transpose(G2)
		G2TG2 = np.matmul( G2T, G2)
		#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
		G2TG2_inv = np.linalg.inv(G2TG2)
		G2_G2TG2_inv = np.matmul(G2, G2TG2_inv)

		S_EXP = np.matmul(G1TG1_inv_G1T, np.matmul(EXP, G2_G2TG2_inv))
		

		OBJ, RELs = self.Score(MOLs, EXP, G1, G2, S_Mol, S_EXP, norm_MOLs, norm_EXP)
		print(" - Init:\t OBJ:%.4f\t REL:"%(OBJ), RELs)

		#plot initial OBJ
		its =[]
		OBJs = []
		its.append(0)
		OBJs.append(OBJ)
		plt.xlim(-0.1,200)
		plt.plot(its, OBJs,  marker='o', color = 'b')
		plt.pause(0.05)

		#Begining M.U.R.
		stop_criterion = 1e-3
		Relative_OBJ_error = 1
		
		for it in range(1,self.max_iter+1):
			if Relative_OBJ_error > stop_criterion:
				#--------------------------------------------------------
			    #                     
			    # First, update all Ss
			    #
			    #-------------------------------------------------------- 
				

				# Computing all S[i] using closed formula
				G1T = np.transpose(G1)
				G1TG1 = np.matmul( G1T, G1)
				G1TG1_inv = np.linalg.inv(G1TG1)
				G1_G1TG1_inv = np.matmul(G1, G1TG1_inv)
				G1TG1_inv_G1T = np.matmul(G1TG1_inv, G1T)
				for i in range(nb_mols):
				    S_Mol[i] = np.matmul(G1TG1_inv_G1T, np.matmul(MOLs[i], G1_G1TG1_inv))


				# computing S for Expression matrices
				#  S_EXP <- (G1T.G1)^-1 G1T.EXP.G2 (G2T.G2)^-1
				
				G2T = np.transpose(G2)
				G2TG2 = np.matmul( G2T, G2)
				G2TG2_inv = np.linalg.inv(G2TG2)
				G2_G2TG2_inv = np.matmul(G2, G2TG2_inv)
				G1TG1_inv_G1T = np.matmul(G1TG1_inv, G1T)
				S_EXP = np.matmul(G1TG1_inv_G1T, np.matmul(EXP, G2_G2TG2_inv))
				
			    #--------------------------------------------------------
			    #                     
			    # Second, update G1, which depends on two decompositions
			    #
			    #--------------------------------------------------------  
			    
			    # G1 = G1.sqrt((G1Num1 + G1Num2) / (G1Denom1 + G1Denom2))
			    
			    #Update rule for G1, contribution from Ai(MOLs[i]) = G1.S[i].G1^T
			    # G1 <- G1.sum(([MOLi.G1.SMoli]_pos + G1.[SMoli^T.G1^T.G1.SMoli]_neg) / ([MOLi.G1.SMoli]_neg + G1.[SMoli^T.G1^T.G1.SMoli]_pos ))
			    
				G1Num1 = np.zeros((g,k1))
				G1Denom1 = np.zeros((g,k1))
			    
				for i in range(nb_mols):
					MOLi_G1_SMoli = np.matmul(MOLs[i], np.matmul(G1, S_Mol[i]))
					MOLi_G1_SMoli_pos = (np.absolute(MOLi_G1_SMoli)+MOLi_G1_SMoli)/2.0
					MOLi_G1_SMoli_neg = (np.absolute(MOLi_G1_SMoli)-MOLi_G1_SMoli)/2.0

					SMoli_G1TG1_SMoli = np.matmul(S_Mol[i], np.matmul(G1TG1, S_Mol[i]))
					SMoli_G1TG1_SMoli_pos = (np.absolute(SMoli_G1TG1_SMoli)+SMoli_G1TG1_SMoli)/2
					SMoli_G1TG1_SMoli_neg = (np.absolute(SMoli_G1TG1_SMoli)-SMoli_G1TG1_SMoli)/2

					G1_SMoli_G1TG1_SMoli_pos = np.matmul(G1,SMoli_G1TG1_SMoli_pos)
					G1_SMoli_G1TG1_SMoli_neg = np.matmul(G1,SMoli_G1TG1_SMoli_neg)

					G1Num1 = G1Num1 + MOLi_G1_SMoli_pos + G1_SMoli_G1TG1_SMoli_neg
					G1Denom1 = G1Denom1 + MOLi_G1_SMoli_neg + G1_SMoli_G1TG1_SMoli_pos
			    
			    #Update rule for G1, contribution from A5(EXP) = G1.S_EXP.G2^T
			    # G1 <- G1 (([EXP.G2.S_EXPT]_pos + G1.[S_EXP.G2TG2.S_EXPT]_neg) / ([EXP.G2.S_EXPT]_neg + G1.[S_EXP.G2TG2.S_EXPT]_pos)) ) 
			    
				EXP_G2_S_EXPT = np.matmul(EXP, np.matmul(G2, S_EXP.T))
				EXP_G2_S_EXPT_pos = (np.absolute(EXP_G2_S_EXPT) + EXP_G2_S_EXPT)/2.0
				EXP_G2_S_EXPT_neg = (np.absolute(EXP_G2_S_EXPT) - EXP_G2_S_EXPT)/2.0    

				S_EXP_G2TG2_S_EXPT = np.matmul(S_EXP, np.matmul(G2.T, np.matmul(G2, S_EXP.T)))
				G1_S_EXP_G2TG2_S_EXPT_pos = np.matmul(G1, ((np.absolute(S_EXP_G2TG2_S_EXPT) + S_EXP_G2TG2_S_EXPT)/2.0))
				G1_S_EXP_G2TG2_S_EXPT_neg = np.matmul(G1, ((np.absolute(S_EXP_G2TG2_S_EXPT) - S_EXP_G2TG2_S_EXPT)/2.0))

				G1Num2 = EXP_G2_S_EXPT_pos + G1_S_EXP_G2TG2_S_EXPT_neg
				G1Denom2 = EXP_G2_S_EXPT_neg + G1_S_EXP_G2TG2_S_EXPT_pos  

				#Check that we are not dividing by zero
				G1update = np.sqrt(np.divide(G1Num1 + G1Num2 + epsilon , G1Denom1 + G1Denom2 + epsilon))
				G1 = np.multiply(G1, G1update)  + epsilon    

				#--------------------------------------------------------
				#                     
				# Now update G2
				#
				#--------------------------------------------------------    

				# Update rule for G2, contribution from A5(EXP) = G1.S_EXP.G2^T
				# G2 <- G2 (([EXPT.G1.S_EXP]_pos + G2.[S_EXP^T.G1^T.G1.S_EXP]_neg) / ([EXPT.G1.S_EXP]_neg + G2.[S_EXP^T.G1^T.G1.S_EXP]_pos))
				EXPT_G1_S_EXP = np.matmul(EXP.T, np.matmul(G1, S_EXP))
				EXPT_G1_S_EXP_pos = (np.absolute(EXPT_G1_S_EXP) + EXPT_G1_S_EXP)/2.0
				EXPT_G1_S_EXP_neg = (np.absolute(EXPT_G1_S_EXP) - EXPT_G1_S_EXP)/2.0    

				S_EXPT_G1TG1_S_EXP = np.matmul(S_EXP.T, np.matmul(G1.T, np.matmul(G1, S_EXP)))
				G2_S_EXPT_G1TG1_S_EXP_pos = np.matmul(G2, ((np.absolute(S_EXPT_G1TG1_S_EXP) + S_EXPT_G1TG1_S_EXP)/2.0))
				G2_S_EXPT_G1TG1_S_EXP_neg = np.matmul(G2, ((np.absolute(S_EXPT_G1TG1_S_EXP) - S_EXPT_G1TG1_S_EXP)/2.0))

				G2Num = EXPT_G1_S_EXP_pos + G2_S_EXPT_G1TG1_S_EXP_neg
				G2Denom = EXPT_G1_S_EXP_neg + G2_S_EXPT_G1TG1_S_EXP_pos

				#Check that we are not dividing by zero  
				G2update = np.sqrt(np.divide(G2Num + epsilon , G2Denom + epsilon))
				G2 = np.multiply(G2, G2update) + epsilon


				if (it%self.verbose == 0) or (it==1):
					OBJ_old = OBJ
					OBJ, RELs = self.Score(MOLs, EXP, G1, G2, S_Mol, S_EXP, norm_MOLs, norm_EXP)
					print(" - It %i:\t OBJ:%.4f\t REL:"%(it, OBJ), RELs)

					its.append(it)
					OBJs.append(OBJ)
					plt.plot(its, OBJs,  marker='o', color = 'b')
					plt.pause(0.05)
					Relative_OBJ_error = np.absolute((OBJ-OBJ_old)/OBJ)
					print(Relative_OBJ_error)
					
		plt.close()
		return G1, G2, S_Mol, S_EXP
		


class Partial_NMF:
	
	# Compute Partial Non-negative Matrix Factorization (P_NMF) of a rectangular matrix n*m A
	# A = PN^T 
	#
	#           n×m                 n×k            k×m           
	#     ┌───────────────┐        ┌───┐     ┌──────────────┐                        
	#     │               │        │   │  x  │      N'      │                 
	#     │               │        │   │     └──────────────┘
	#     │      A        │    ≈   │ P │ 
	#     │               │        │   │
	#     │               │        │   │
	#     │               │        │   │
	#     └───────────────┘        └───┘
	#
	#  With: A and P given by the user, N >= 0
	
	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		
	
	def Score(self, A, P, N, norm_A):
		NT = np.transpose(N)
		
		Err1 = A - np.matmul(P, NT)
		norm_err1 = np.linalg.norm(Err1, ord='fro')
		
		rel_err1 = norm_err1 / norm_A
		
		return norm_err1, rel_err1
	
	
	def Solve_MUR(self, A, P, k):
		print("Starting NMTF")
		n,m=np.shape(A)
		norm_A = np.linalg.norm(A, ord='fro')
		
		print(" * Using random initialization")
		N = np.random.rand(m,k)
		
		OBJ, REL = self.Score(A, P, N, norm_A)
		print(" - Init:\t OBJ:%.4f\t REL:%.4f"%(OBJ, REL))
		
		
		#Begining M.U.R.
		AT = np.transpose(A)
		AT_P = np.matmul(AT, P)
		PT_P = np.matmul(np.transpose(P), P)
		
		for it in range(1,self.max_iter+1):
			
			#update rule for N
			NT = np.transpose(N)
			
			N_PT_P = np.matmul(N, PT_P)
			N_mult = np.divide(AT_P, N_PT_P + epsilon)
			
			# Applying M.U.R.
			N = np.multiply(N, N_mult) + epsilon
			
			if (it%self.verbose == 0) or (it==1):
				OBJ, REL = self.Score(A, P, N, norm_A)
				print(" - It %i:\t OBJ:%.4f\t REL:%.4f"%(it, OBJ, REL))
		return N






class ISF_PONMF:
	#
	#  EXP = P.S1.B^T
	#  Mol[i] = G[i].S[i].B^T
	#  subject to P, S, B, G >=0, and B^T.B = I
	
	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		
	
	def Score(self, EXP, MOLs, P, B, G, S1, S, norm_EXP, norm_Mol):
		BT = np.transpose(B)
		
		Square_Error = 0.
		
		Err_EXP = EXP - np.matmul(P,  np.matmul(S1, BT))
		norm_Err_EXP = np.linalg.norm(Err_EXP, ord='fro')
		rel_EXP = norm_Err_EXP / norm_EXP
		Square_Error += norm_Err_EXP
		
		norm_Err_MOLs = []
		rel_MOLs = []
		for i in range(len(MOLs)):
			Err_MOLs = MOLs[i] - np.matmul(G[i],  np.matmul(S[i], BT))
			norm_Err_MOLs.append( np.linalg.norm(Err_MOLs, ord='fro'))
			rel_MOLs.append( norm_Err_MOLs[i] / norm_Mol[i] )
			Square_Error += norm_Err_MOLs[i]
				
		return Square_Error, rel_EXP, rel_MOLs
	
	
	
	def Solve_MUR(self, EXP, MOLs, k_p, k_g, init="rand"):
		n_p,n_g = np.shape(EXP)
		nb_mols = len(MOLs)
		
		norm_EXP = np.linalg.norm(EXP, ord='fro')
		norm_MOLs = [np.linalg.norm(MOLs[i], ord='fro') for i in range(nb_mols)]
		
		if init == "rand":
			print(" * Using random initialization")
			P = np.random.rand(n_p, k_p)
			B = np.random.rand(n_g, k_g)
			G = [np.random.rand(n_g, k_g) for i in range(nb_mols)]
			
			S1 = np.random.rand(k_p,k_g)
			S = [np.random.rand(k_g,k_g) for i in range(nb_mols)]
			
		elif init == 'SVD':
			print(" * Using SVD based initialization")
			P = np.zeros((n_p, k_p)) + epsilon
			B = np.zeros((n_g, k_g)) + epsilon
			G = [np.zeros((n_g, k_g)) + epsilon for i in range(nb_mols)]
			
			S1 = np.zeros((k_p,k_g)) + epsilon
			S = [np.zeros((k_g,k_g)) + epsilon for i in range(nb_mols)]
			
			nbf = 1.+float(nb_mols)
			
			print(" * -- SVD on EXP to get P.S1.B^T")
			J,K,L = np.linalg.svd(EXP)
			for i in range(n_p):
				for j in range(min(n_p, k_p)):
					if J[i][j] > 0.:
						P[i][j] = J[i][j]
			for i in range(min(min(k_p, n_p), k_g)):
				S1[i][i] = K[i]
			for i in range(n_g):
				for j in range(k_g):
					if L[j][i] > 0.:
						B[i][j] = L[j][i]/nbf
			
			for xx in range(nb_mols):
				print(" * -- Eig decomposition on MOL[%i] to get Gi.Si.B^T"%(xx+1))
				M, N = ssl.eigsh(MOLs[xx], k_g)
				for i in range(n_g):
					for j in range(k_g):
						if N[i][j] > 0.:
							G[xx][i][j] = N[i][j]
				for i in range(k_g):
					S[xx][i][i] = abs(M[i])
				for i in range(n_g):
					for j in range(k_g):
						if N[i][j] > 0.:
							B[i][j] += N[i][j]/nbf
			
			
		
		else :
			print("Unknown initializer: %s"%(init))
			exit(0)
		OBJ, REL1, RELs = self.Score(EXP, MOLs, P, B, G, S1, S, norm_EXP, norm_MOLs)
		print(" - Init:\t OBJ:%.4f\t REL:"%(OBJ), REL1, RELs)
		
		
		#Begining M.U.R.
		
		EXPT = np.transpose(EXP)
		MolT = [np.transpose(Mol) for Mol in MOLs]
		
		for it in range(1,self.max_iter+1):
			
			### First, update B, which depends on multiple decompositions
			
			# Update rule for B, orthonormal, contribution from EXP = P.S1.B^T
			BT = np.transpose(B)
			Numb = np.matmul(np.matmul(EXPT, P), S1)
			Denb = np.matmul(B, np.matmul(BT, Numb))
			# Update rule for B, orthonormal, contribution from all MOL[i] = Gi.Si.B^T
			for i in range(nb_mols):
				Num3 = np.matmul(np.matmul(MolT[i], G[i]), S[i])
				Numb += Num3
				Denb += np.matmul(B, np.matmul(BT, Num3))
			#Composed update rule
			B_mult = np.sqrt( np.divide(Numb, Denb + epsilon) )
			# Applying M.U.R. on B
			B = np.multiply(B, B_mult) + epsilon
			
			
			### Second, update all matrix factors that depend only on one decomposition, i.e., P, G, S1, S2, S3, and S4
			
			# Update rule for P, not orthonormal, from EXP=P.S1.B^T
			BT = np.transpose(B)
			S1T = np.transpose(S1)
			B_S1T = np.matmul(B, S1T)
			EXP_B_S1T = np.matmul(EXP, B_S1T)
			P_S1_BT = np.matmul(P, np.matmul(S1, BT))
			P_S1_BT_B_S1T = np.matmul(P_S1_BT, B_S1T)
			P_mult = np.sqrt( np.divide(EXP_B_S1T, P_S1_BT_B_S1T + epsilon) )
			# Applying M.U.R. on P
			P = np.multiply(P, P_mult) + epsilon
			
			#update rule for S1, from EXP=P.S1.B^T
			PT = np.transpose(P)
			PT_EXP_B = np.matmul(PT, np.matmul(EXP, B))
			PT_P = np.matmul(PT, P)
			BT_B = np.matmul(BT, B)
			PT_P_S1_BT_B = np.matmul( PT_P, np.matmul(S1, BT_B))
			S1_mult = np.sqrt( np.divide(PT_EXP_B,PT_P_S1_BT_B + epsilon))
			# Applying M.U.R. on S1
			S1 = np.multiply(S1, S1_mult) + epsilon
					
			# Update rule for all G[i], not orthonormal, from MOL[i] = Gi.Si.B^T
			for i in range(nb_mols):
				SiT = np.transpose(S[i])
				B_SiT = np.matmul(B, SiT)
				Mol_B_SiT = np.matmul(MOLs[i], B_SiT)
				G_Si_BT = np.matmul(G[i], np.matmul(S[i], BT))
				G_Si_BT_B_SiT = np.matmul(G_Si_BT, B_SiT)
				G_mult = np.sqrt( np.divide(Mol_B_SiT, G_Si_BT_B_SiT + epsilon) )
				# Applying M.U.R. on G
				G[i] = np.multiply(G[i], G_mult) + epsilon
			
			
			#update rule for all Si, from MOL[i] = Gi.Si.B^T
			for i in range(nb_mols):
				GT = np.transpose(G[i])
				GT_Mol_B = np.matmul(GT, np.matmul(MOLs[i], B))
				GT_G = np.matmul(GT, G[i])
				#BT_B = np.matmul(BT, B)
				GT_G_Si_BT_B = np.matmul( GT_G, np.matmul(S[i], BT_B))
				Si_mult = np.sqrt( np.divide(GT_Mol_B,GT_G_Si_BT_B + epsilon))
				# Applying M.U.R. on S3
				S[i] = np.multiply(S[i], Si_mult) + epsilon
			
			
			if (it%self.verbose == 0) or (it==1):
				OBJ, REL1, RELs = self.Score(EXP, MOLs, P, B, G, S1, S, norm_EXP, norm_MOLs)
				print(" - Init:\t OBJ:%.4f\t REL:"%(OBJ), REL1, RELs)
		return P, B, G, S1, S



#test
#n_g = 2000
#n_p = 500
#
#EXP = np.random.rand(n_p,n_g)
#Mol1 = np.random.rand(n_g,n_g)
#Mol1 = (Mol1+np.transpose(Mol1))/2.
#
#Mol2 = np.random.rand(n_g,n_g)
#Mol2 = (Mol2+np.transpose(Mol2))/2.
#
#Mol3 = np.random.rand(n_g,n_g)
#Mol3 = (Mol3+np.transpose(Mol3))/2.
#
#
#Solver = ISF_PONMF(max_iter=500, verbose = 10)
#P, B, G, S1, S = Solver.Solve_MUR(EXP, [Mol1, Mol2, Mol3], 50, 50, init='rand')



