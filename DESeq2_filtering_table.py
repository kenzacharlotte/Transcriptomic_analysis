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
