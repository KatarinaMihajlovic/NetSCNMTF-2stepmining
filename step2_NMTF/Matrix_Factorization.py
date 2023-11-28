# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""

import numpy as np
import scipy.sparse.linalg as ssl

import warnings
warnings.simplefilter("ignore", UserWarning)



epsilon = 1e-5


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

			SMOL = np.zeros((g,g)) + epsilon
			for i in range(nb_mols):
				SMOL += MOLs[i]
			print(np.count_nonzero(SMOL), SMOL.shape, SMOL.dtype)

# 			for x in SMOL:
# 				for i in x:
# 					if i != 0.0:
# 						if i != 1.0:
# 							print(i)


			if np.linalg.norm(SMOL, ord='fro') != 0:
				print(" * -- Eig decomposition on MOL to get G1.Si.G1^T")
				#print(np.linalg.norm(SMOL, ord='fro'), SMOL.ndim, SMOL.shape, type(SMOL))
				try:
					M, N, O = ssl.svds(SMOL, k1, solver='lobpcg')
					M = abs(M)
					#print(M)
					for i in range(g):
						for j in range(k1):
							G1[i][j] += M[i][j]
				except ValueError:
					pass

			
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
		
		fig = plt.figure(figsize=(16, 9))
		ax = fig.add_subplot(111)
		ax.set_xlabel('Iterations', fontsize = 24)
		ax.set_ylabel('Value of Objective Function', fontsize = 24)
		ax.tick_params(axis='y', which='major', labelsize=20)
		ax.tick_params(axis='x', which='major', labelsize=20)
		#ax.set_xlim(left=-0.2, right=None) 
		ax.plot(its, OBJs,  marker='o', color = 'b')
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
					ax.plot(its, OBJs,  marker='o', markersize=15, color = 'b')
					plt.pause(0.05)
					Relative_OBJ_error = np.absolute((OBJ-OBJ_old)/OBJ)
					print(Relative_OBJ_error)
		

		#plt.close()
		return G1, G2, S_Mol, S_EXP, fig
		

class NMTF_basic:
	#
	#  Mol[i] = G.S[i].G^T
	#  subject to G >=0

	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		

	def Score(self, EXP, G1, G2, S_EXP, norm_EXP):
		Square_Error = 0.    

		rel_EXP = []
		G2T = np.transpose(G2)
		Err_EXP = EXP - np.matmul(G1,  np.matmul(S_EXP, G2T))
		norm_Err_EXP = np.linalg.norm(Err_EXP, ord='fro')
		rel_EXP.append(norm_Err_EXP / norm_EXP)

		Square_Error += norm_Err_EXP**2

		return Square_Error, rel_EXP
		
	
	def Solve_MUR(self, EXP, k1, k2, init="rand"):
		import matplotlib.pyplot as plt

		g, c = np.shape(EXP)
		print(g,c)
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
			
			SEXP = EXP
			# A = J * K * L
			if k1 > k2:            
				J, K, L = ssl.svds(SEXP, k1, solver='lobpcg')
				#J, K, L = sl.svd(SEXP, k1)

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
        # Initialize S_EXP
        #
        #-------------------------------------------------------- 

		#computing Ss for Molecular matrices
		#  S_Mol[i] <- (G1T.G1)^-1 G1T.MOLs[i].G1 (G1T.G1)^-1

		# computing S for Expression matrices
		#  S_EXP <- (G1T.G1)^-1 G1T.EXP.G2 (G2T.G2)^-1
		G1T = np.transpose(G1)
		G1TG1 = np.matmul( G1T, G1)
		#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
		G1TG1_inv = np.linalg.inv(G1TG1)
		G1_G1TG1_inv = np.matmul(G1, G1TG1_inv)
		G1TG1_inv_G1T = np.matmul(G1TG1_inv, G1T)
        
		G2T = np.transpose(G2)
		G2TG2 = np.matmul( G2T, G2)
		#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
		G2TG2_inv = np.linalg.inv(G2TG2)
		G2_G2TG2_inv = np.matmul(G2, G2TG2_inv)

		S_EXP = np.matmul(G1TG1_inv_G1T, np.matmul(EXP, G2_G2TG2_inv))
		

		OBJ, RELs = self.Score(EXP, G1, G2, S_EXP, norm_EXP)
		print(" - Init:\t OBJ:%.4f\t REL:"%(OBJ), RELs)

		#plot initial OBJ
		its =[]
		OBJs = []
		its.append(0)
		OBJs.append(OBJ)
		
		fig = plt.figure(figsize=(16, 9))
		ax = fig.add_subplot(111)
		ax.set_xlabel('Iterations', fontsize = 24)
		ax.set_ylabel('Value of Objective Function', fontsize = 24)
		ax.tick_params(axis='y', which='major', labelsize=20)
		ax.tick_params(axis='x', which='major', labelsize=20)
		#ax.set_xlim(left=-0.2, right=None) 
		ax.plot(its, OBJs,  marker='o', color = 'b')
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
			    

			    #Update rule for G1, contribution from A5(EXP) = G1.S_EXP.G2^T
			    # G1 <- G1 (([EXP.G2.S_EXPT]_pos + G1.[S_EXP.G2TG2.S_EXPT]_neg) / ([EXP.G2.S_EXPT]_neg + G1.[S_EXP.G2TG2.S_EXPT]_pos)) ) 
			    
				EXP_G2_S_EXPT = np.matmul(EXP, np.matmul(G2, S_EXP.T))
				EXP_G2_S_EXPT_pos = (np.absolute(EXP_G2_S_EXPT) + EXP_G2_S_EXPT)/2.0
				EXP_G2_S_EXPT_neg = (np.absolute(EXP_G2_S_EXPT) - EXP_G2_S_EXPT)/2.0    

				S_EXP_G2TG2_S_EXPT = np.matmul(S_EXP, np.matmul(G2.T, np.matmul(G2, S_EXP.T)))
				G1_S_EXP_G2TG2_S_EXPT_pos = np.matmul(G1, ((np.absolute(S_EXP_G2TG2_S_EXPT) + S_EXP_G2TG2_S_EXPT)/2.0))
				G1_S_EXP_G2TG2_S_EXPT_neg = np.matmul(G1, ((np.absolute(S_EXP_G2TG2_S_EXPT) - S_EXP_G2TG2_S_EXPT)/2.0))

				G1Num = EXP_G2_S_EXPT_pos + G1_S_EXP_G2TG2_S_EXPT_neg
				G1Denom = EXP_G2_S_EXPT_neg + G1_S_EXP_G2TG2_S_EXPT_pos  

				#Check that we are not dividing by zero
				G1update = np.sqrt(np.divide(G1Num + epsilon , G1Denom + epsilon))
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
					OBJ, RELs = self.Score(EXP, G1, G2, S_EXP, norm_EXP)
					print(" - It %i:\t OBJ:%.4f\t REL:"%(it, OBJ), RELs)

					its.append(it)
					OBJs.append(OBJ)
					ax.plot(its, OBJs,  marker='o', markersize=15, color = 'b')
					plt.pause(0.05)
					Relative_OBJ_error = np.absolute((OBJ-OBJ_old)/OBJ)
					print(Relative_OBJ_error)
		

		#plt.close()
		return G1, G2, S_EXP, fig
		


