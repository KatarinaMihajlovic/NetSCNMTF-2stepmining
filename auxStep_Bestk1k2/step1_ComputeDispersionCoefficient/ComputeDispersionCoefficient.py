import pandas as pd
import numpy as np
import os, pickle,sys
from sklearn.cluster import KMeans
from sklearn_extra.cluster import KMedoids


def kMeans(M, n_clusters):
    M = M.tolist()
   
    kMeans = KMeans(n_clusters=n_clusters).fit(M)
    KMeans_labels = list(kMeans.labels_)

    return KMeans_labels

def kMedoids(M, n_clusters):
    kmedoids = KMedoids(n_clusters=n_clusters, init = 'random').fit(M)
    kmedoids_labels = list(kmedoids.labels_)
    return kmedoids_labels


def getHardClusters(G):
    """
    Compute hard clustering from a G factor comming from a NMTF factorization
    
    Parameters
    ----------
    G: numpy array
        Contains G factor (cluster indicator) of a NMTF factorization
        
    Return
    ------
    Returns the cluster indicator of each entry in a list.
    
    """

    return [np.argmax(G[i]) for i in range(G.shape[0])]

def createConnectivityMatrix(clusters, k, col):
    """
    Create the connectivity matrix
    
    Parameters
    ----------
    clusters: DataFrame
        Each column contain the cluster indicators for the elements in the index
    k: integer
        Indicates the number of clusters
    col: string
        Indicates the column from which the connectivy matrix must be computed
        
    Return
    ------
    Returns the created connectivity matrix in a DataFrame.
    
    """

    
    connectivityMatrix = pd.DataFrame(np.zeros((clusters.shape[0], clusters.shape[0])), \
                                      index=clusters.index, columns=clusters.index)
    for nC in range(k):
        nodes = clusters[clusters[col] == nC].index
        cMTmp = pd.DataFrame(np.ones((nodes.shape[0], nodes.shape[0])), index=nodes, columns=nodes)
        connectivityMatrix.update(cMTmp)

    return connectivityMatrix.values

def copeheneticCorreletionCoef(matrix):
    """
    Computed the copehenetic correlation coefficient from a matrix
    
    Parameters
    ----------
    matrix: numpy
        contains the matrix with the clusters (consensus matrix)
        
    Return
    ------
    Returns the correlaction coefficient.
    
    """
    
    n = matrix.shape[0]
    coef = 0
    for row in matrix:
        for item in row:
            coef += 4*(item - (1/2))**2
    return 1/n**2 * coef

def computeDispersionCoefficient(clusterElements, cell_cond, k1s = [25,50,75,100,125,150,175,200,250,275], k2s = [10,20,30,40,50,60,70,80,90,100], R=10, in_dir='output/', savePath='output/', verbose=True, cur_dir = os.getcwd()):
    """
    Compute the dispersion coefficient for each parameter from the cluster indicatiors factors of a GNMTF factorization for different parametrizations and runs. It saves the clusters of the different runs in a same file and the consensus matrix computed.
    
    Parameters
    ----------
    clusterElements: list of lists
        Indicates the names of the elements of the different cluster-indicator matrices (genes and SCs)
    k1s: list
        Indicated the values for the parameter k1. By default [25,50,75,100,125,150,175,200,250,275] 
    k2s: list
        Indicated the values for the parameter k2. By default [10,20,30,40,50,60,70,80,90,100]
    R: list
        Indicated the number of runs. By default 10
    in_dir: string
        Indicated the in_dir where the factor matrices are stored.    
    savePath: string
        Indicated the path where the results must be saved.
        
    verbose: boolean
        Indicates whether information of the dataset have to be printed. By default True

        

    """

    Net_comb = 'ALL'
    res= []

    cm = 'kmeans'

    for k1 in k1s:
        for k2 in k2s:
            save_path2 = savePath + '/k1_' + str(k1) + '_k2_' + str(k2)
            try:
                os.makedirs(save_path2)
            except FileExistsError:
                pass
            
            res_file = open(save_path2 + '/DispCoeff_' + cell_cond + '_k1_' + str(k1) + '_k2_' + str(k2) + '_' + cm + '.txt', 'w')
            res_file.write('Net_comb\tk1\tk2\tdispersionCoefG1\tdispersionCoefG2\n')
            if verbose: print('k1=' + str(k1) + ', k2=' + str(k2))

            path_parent = os.path.dirname(os.getcwd())
            os.chdir(path_parent) 

            #Compute the connectivity and consensus matrices from the G factors

            #Initializate the clusters and consensus matrices
            clustersG1 = pd.DataFrame(index=clusterElements[0])
            clustersG2 = pd.DataFrame(index=clusterElements[1])

            consensusMatrixG1 = np.zeros((len(clusterElements[0]), len(clusterElements[0])))
            consensusMatrixG2 = np.zeros((len(clusterElements[1]), len(clusterElements[1])))
            for r in range(R):
                if verbose:  print('Run: ', r)
                vars_dir = in_dir + '/k1_' + str(k1) + '_k2_' + str(k2) + '/' + str(r)  + '/' + Net_comb
                G1 = np.load(vars_dir + '/' + Net_comb + '_G1.npy')
                G2 = np.load(vars_dir + '/' + Net_comb + '_G2.npy')

                G1_NumClusters = np.shape(G1)[1]
                genes = clusterElements[0]
                G2_NumClusters = np.shape(G2)[1]
                SCs = clusterElements[1]
                
                clustersG1.insert(clustersG1.shape[1], str(r), kMeans(G1, G1_NumClusters))
                clustersG2.insert(clustersG2.shape[1], str(r), kMeans(G2, G2_NumClusters))


                if verbose:  print(' CM-G1 Completed')
                consensusMatrixG1 += createConnectivityMatrix(clustersG1, k1, str(r))
                if verbose:  print(' CM-G2 Completed')
                consensusMatrixG2 += createConnectivityMatrix(clustersG2, k2, str(r))

            #Save the clusters
            out_dir = savePath + '/k1_' + str(k1) + '_k2_' + str(k2)
            try:
                os.makedirs(out_dir)
            except FileExistsError:
                pass
            clustersG1.to_csv(out_dir + '/G1_clusters_' + cm + '.csv')
            clustersG2.to_csv(out_dir + '/G2_clusters_' + cm + '.csv')

            #Normalize the consensus matrices
            consensusMatrixG1 = consensusMatrixG1 / R
            consensusMatrixG2 = consensusMatrixG2 / R

            #Saving the consensus Matrices - SAVING IT REQUIRES TOO MUCH HARD DISK MEMORY
            #np.save(out_dir + '/G1_ConsensusMatrix.npy', consensusMatrixG1)
            #np.save(out_dir + '/G2_ConsensusMatrix.npy', consensusMatrixG2)

            #Compute Dispersion Coefficient
            dispersionCoefG1 = copeheneticCorreletionCoef(consensusMatrixG1)
            dispersionCoefG2 = copeheneticCorreletionCoef(consensusMatrixG2)

            res_file.write(Net_comb + '\t' + str(k1) + '\t' + str(k2) + '\t' + str(dispersionCoefG1) + '\t' + str(dispersionCoefG2) + '\n')
            res_file.close()



















