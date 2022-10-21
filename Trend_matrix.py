import numpy as np

def Trend(matrixfile,outputfile,k):
    """
    A function that transform a numeric matrix in a binary matrix.
    Useful to built a trend matrix
    
    If : 
        (*) L2FC[t-1] < L2FC[t] => binMatrix = 1
        (*) L2FC[t-1] > L2FC[t] => binMatrix = -1
    
    PARAMETER
    - matrixfile : (str) path to the .csv file text that contains a numeric matrix
    - outputfile : (str) path of the output file
    - k : (float) k > 0  is threshold of the affine function
        Only the genes that have a highest/lowest coefficient rate are selected.
        It's a way to select the significant variation between 2 timepoints. 
    
    For example a L2FC matrix : 
         (*) T : the number of timepoints
         (*) M : the number of genes
    
    #### Input file structure ##############################
    
    #geneNames #L2FC_t0 #L2FC_t1   ...  #L2FC_t+n 
       #1        ...     ...       ...     ...         
      ....       ...     ...       ...     ...     
      ....       ...     ...       ...     ...    
      ....       ...     ...       ...     ...     
       #M        ...     ...       ...     ...     
      
    ########################################################
    
    RETURN
    - binmatrix : (np.matrix) an np.array binary matrix size(M*T-1)
    
    Write this matrix in a csv file txt 
    """
    matrix = np.loadtxt(matrixfile,delimiter = ',')[:,1:]
    M,T = np.shape(matrix) # M genes and T timepoints
    binmatrix = np.zeros((M,T-1)) # Building the binary matrix

    for gene in range(M):
        for timepoint in range(T-1):
            if (matrix[gene,timepoint+1]-matrix[gene,timepoint])>k:
                binmatrix[gene,timepoint]=1
            elif (matrix[gene,timepoint+1]-matrix[gene,timepoint])<-k:
                binmatrix[gene,timepoint]=-1
            else :
                binmatrix[gene,timepoint]=0
                
    np.savetxt(outputfile, binmatrix,fmt='%i', delimiter=',')
    return(binmatrix)
