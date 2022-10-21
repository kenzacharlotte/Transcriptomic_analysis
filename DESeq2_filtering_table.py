############## STEP I ###############################################
# Filter the significative L2FC in one timepoint
# DESeq2 output

def FilterTable(file, name, padjThreshold = 0.05, L2FCThreshold = 1):
    """
    Function that reads the input file, select and save in a text file the genes with : 
            (1) a significative padj ie. padj < padj threshold
            (2) a signifcative Log2 fold change ie.  |L2FC| > L2FC threshold 
    
    Output file structure : 
    #geneNames #L2FC #padj
    
    NB : 1 file for 1 condition/timepoint
            
    ### WARNING #############################################################
    The input file (from R DESeq2 package) should be structered like this : 
    
    - 1st column = gene name
    - 3nd column = log2 fold change
    - last column = padj
    
    > special case for the header : 
    - 2nd element = Log2FoldChange
    - last element = padj
    #######################################################################
    
    PARAMETERS  
    - file : (str) path to the file
    - name : (str) name of the output file
    - padjThreshold : (int or float) the padj threshold (0.05 by default)
    - L2FCThreshold : (int or float) the L2FC threshold (1 by default)
    
    RETURN 
    - FilteredDict : (dict) keys are the genes name associated to a 2 element list
    This list contain the L2FC and padj-value associated to the gene that are significant
    
    Write a file text containing the filtered values 
    """
    
    # I - Initialisation
    rf = open(file, 'r')
    wf = open('FilteredTable.'+name+'.padj_'+str(padjThreshold)+'.L2FC1_'+str(L2FCThreshold)+".txt",'w')
    cpt = 0
    cptF = 0
    FilteredDict = {}
    
    # II - Filtering
    # (i) - Reading the file
    for line in rf :
        line = line.split() #each line is splitted in a list of 7 element
        # first element of the list (index = 0) : gene name
        # 3d element of the list (index = 2) : L2FC
        # last element of the list (index = -1) : padj
        
        # (ii) - Particular case for the header
        if cpt ==0 : 
            wf.write('GeneName\t'+line[1]+'\t'+line[-1]+'\n')
            FilteredDict['GeneName'] = [line[1],line[-1]]
            
        # (iii) - Filtering based on the padj value and L2FC
        else :
            if (line[-1] != 'NA') : #excluding the NA values 
                if (float(line[-1]) < padjThreshold and (abs(float(line[2])) >= L2FCThreshold)): 
                    #we're taking only the values that are in between the thresholds
                    wf.write(line[0]+'\t'+line[2]+'\t'+line[-1]+'\n') #saving it in a file
                    FilteredDict[line[0]] = [line[2],line[-1]] #saving it in a dictionary
                    cptF+=1 #counting gene number
        cpt+=1 
    wf.write('\nNumber of genes filtered in : '+str(cptF))
    return FilteredDict 


############## STEP II ###############################################
# From a merged table : padj + L2FC accross all the timepoints
# Select only the genes that have a signicative L2FC in AT LEAST one timepoint

def FilterMergedTable(file, name, ListDict):
    """
    Function that read the input file and select the genes that have a significant differential expression
    in AT LEAST ONE time point.
            
    ### WARNING #############################################################
    * 1 file for all the timepoints
    * The input file is a merged file from DeSeq2 R package  
    #########################################################################
    
    PARAMETERS  
    - file : (str) path to the file
    - name : (str) name of the output file
    - ListDict : (list) containing all the dictionnary from the function FilterTable(file, name, padjThreshold, L2FCThreshold)
    
    RETURN 
    - FilteredGenes : (list) containing the names of all the genes that have been selected
    Write a file text containing the filtered genes & values accross all the timepoints

    ### Output file structure ################################################
    
    #geneNames #L2FC_t0 #padj_t0 #L2FC_t1 #padj_t1 ... #L2FC_t+n #padj_t+n
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      
    #########################################################################


    """
    FilteredGenes = []
    nbCheck = 0
    # 0 - Extracting all the dictionnary keys 
    # To extract all the genes that have a significant differential expression in at least 1 timepoiny
    for FilteredDic in ListDict : 
        for gene in FilteredDic.keys() : 
            if gene not in FilteredGenes :
                FilteredGenes.append(gene)

    # I - Initialisation
    rf = open(file, 'r')
    wf = open('AllSignificativeGenesTable.'+name+".txt",'w')
    cpt = 0
    
    # II - Filtering
    
    # (i) - Reading the file
    for line in rf :
        lineS = line.split() #each line is splitted in a list of (n*2)+1 element
        # first element of the list (index = 0) : gene name
        # n the number of timepoints, each timepoint has 2 values : L2FC and padj
        
        # (ii) - Particular case for the header
        if cpt ==0 : 
            wf.write(line+'\n')

        # (iii) - Filtering based on the padj value and L2FC
        else :
            if (lineS[0] in FilteredGenes) : #writing only the genes that has been filtered-in in the new file
            #we're taking only the values that are significant in at least one time point
                nbCheck +=1 
                wf.write(line) #saving all the line in a file
       
        cpt+=1 
    print('\nNumber of genes filtered-in accross all the timepoints : '+str(nbCheck)+'\n')
    return FilteredGenes

############## STEP III ###############################################
# From a file with all the filtered genes and their L2FC+padj values
# Creates a matrix containing only the L2FC ready for K-means clustering algorithm

def KmeanMatrix(file, ofile):
    """
    Function that from a file creates a file to use for K-mean clusterisation. 
    
    PARAMETERS
    - file : (str) path to the file
    - ofile : (str) name of the output file 
    
    ### Intput file structure ################################################
    #geneNames #L2FC_t0 #padj_t0 #L2FC_t1 #padj_t1 ... #L2FC_t+n #padj_t+n
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
      ....       ...     ...       ...     ...     ...     ...      ...
    ##########################################################################    

    """
    # I - Initialisation
    rf = open(file, 'r')
    wf = open('L2FC_Matrix.'+ofile+".txt",'w')
    cpt = 0
    NA_nb = 0
    G_nb = 0
    
    # II - Filtering-out the padj values
    # (A) - Reading the file
    for line in rf :
        lineS = line.split() 
        
        # (A)(i) - Particular case for the header
        if cpt == 0 :
            header = "#"+lineS[0]
            for i in range(1,len(lineS),2) : 
                header += "," + lineS[i]
            header += '\n'
            wf.write(header)
            
        # (A)(ii) - Saving the L2FC accross all the timepoints in a file
        else :
            if "NA" not in lineS : 
                if len(lineS)!=0 :
                    row = lineS[0]
                    for i in range(1,len(lineS),2) : 
                        row += ","+ lineS[i]
                    row += '\n'
                    wf.write(row)
                    G_nb +=1
            if ("NA" in lineS) or (len(lineS)==0) : 
                NA_nb +=1
        cpt+=1
    wf.close
    rf.close
    
    print("Some informations :\n(*)",NA_nb," NA values\n(*)",G_nb," genes ready to use for k-means clusterisating\n")
    return 
