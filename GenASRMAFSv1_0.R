#myAln <- "~/c04_demo/alignment_c04.fasta"
#myMCC <- "~/c04_demo/tree.newick.txt"
#myASR <- "~/c04_demo/seq.joint.txt" #fastml with -qf option
#myMAF <- "~/c04_demo/anno_c04.tsv.txt"
#normName <- "c04_N6"
#patientID <- "c04"



##Function: gen_ASR_MAF
###Purpose: Return to the user a MAF with the mutations from a MAF mapped to their first occurance, branchwise.
###Mandatory Input
#### myAln, the fasta alignment used to make the tree
#### myMCC, the newick tree with internal nodes labeled. FASTML does this by default (tree.newick.txt)
#### myASR, the ancestral state reconstruction file (v1 supports only fastml)
#### normName, the name of the normal sample to in the alignment file to check against the MAF reference.
###Optional Input
#### patientID, if included, a column with the patientID (not to be confused with sampleID!) will be included
###### as a column. Helpful if going to combine many MAFS later. 
#### iq_or_fast, is the ancestral state iqtree ("iq") or fastml ("fast"). Need to work out a bug with iqtree so
###### right now the only option that will work is "fast" unless also providing a internally labelled newick tree
###### corresponding your iqtree input (in which case will work)
#### keepExtraMAF, if TRUE will save all the non-essential MAF columns from the input MAF and include them in
###### the returned MAF data structure
#### outFile, if a character string, will save the MAF to file in addition to returing the data frame.
###Depends
#### ape, seqinr
###Output
#### a data frame and, optionally, a maf file which include the maf elements + new annotations
#### Origin_Node, the first ancestral node the variant appears in. Not that helpful
#### Branch_Number, the number of the branch where the variant first appears. is consistent with a dfs
#### Binding_Nodes, the 2 nodes that bind the branch.
###Other Notes
####ASSUMES THE ALIGNMENT AND MAF ARE IN THE SAME ORDER

gen_ASR_MAF <- function(myAln,myMCC,myASR,myMAF,normName,patientID=NULL,iq_or_fast = NULL,keepExtraMAF = F,outFile=NULL) {
  #load requires  
  library("ape")
  library("seqinr")
  
  #Used to print if there are any terminal leaf ambiguities and, if so, how many there are!
  leafAmbigFlag<-F
  leafAmbigCount<-0
    
  #read in the original alignment, remove any blank lines
  origAln <- readLines(myAln)
  origAln <- origAln[nchar(origAln) > 1]
    
  ###Prep the ASR###
  #if the input ASR is iq tree  
  if (iq_or_fast == "iq") {
    asrAln <- origAln #start with the original alignment
    tmpASR <- read.table(myASR, header = TRUE) #read in the ASR
    uniNodes <- unique(tmpASR$Node) #get all internal nodes
    for (i in uniNodes) { #for each internal node
      asrAln <- c(asrAln, paste0(">", i)) #add name to alignment fasta
      asrAln <-
        c(asrAln, paste0(tmpASR$State[which(tmpASR$Node == i)], collapse = "")) #add seq
    }#if it is fastml
  } else if (iq_or_fast == "fast") {
    asrAln <- readLines(myASR) #easy, they already have it in the format we want!
  } else{
    stop("must select either iqtree (iq) or fastml (fast)") #catch if user doesn't specify which ASR
  }
  #we write the plaintext to a file, as seqinr doesn't work with internal con
  #to avoid potential collisions from users running code in parallel in same working dir
  #generate a random number and prepend it to the filenames
  rNameNum<-as.character(sample(1:1000000000,1)) #technically collisions still possible but 1 in 1 billion chance
  writeLines(asrAln, con = paste0("./",rNameNum,"-asrAln.fasta"))
  #read Ancestral Alignment as fasta with seqinr
  asrAln <- read.alignment(paste0("./",rNameNum,"-asrAln.fasta"), format = "fasta")
  #tidy up by removing the tmp file
  removedASRFile<-file.remove(paste0("./",rNameNum,"-asrAln.fasta"))
  
  ###Read in the tree
  thisTree <- read.tree(myMCC)
  ###ASSUMES THE ALIGNMENT AND MAF ARE IN THE SAME ORDER###
  
  ###Read in the MAF (which is just a tsv)
  origMAF <- read.csv(myMAF, sep = "\t")
  
  ###Check MAF for having the minimal essential MAF elements###
  ###define those elements
  mafCheck <-
    c(
      "Sample_ID",
      "Chromosome",
      "Start_Position",
      "Reference_Allele",
      "Tumor_Seq_Allele2"
    )
  #flag if there is an issue in the inital passover
  issueFlag <- F
  for (i in mafCheck) {
    if (!i %in% colnames(origMAF)) {#if the column names don't match
      if (i != "Sample_ID") {issueFlag <- T; print(paste0(i, " not in origMaf!"))}
    }
  }
  #Alerted user, now try and fix with common alternatives (from vcf files)
  if (issueFlag) {
    print("There was an issue in the MAF column names")
    print("Expected column names are")
    print(mafCheck)
    print("Attempting to autofix MAF column names")
    colnames(origMAF) <-
      gsub(pattern = "ref", "Reference_Allele", colnames(origMAF))
    colnames(origMAF) <-
      gsub(pattern = "alt", "Tumor_Seq_Allele2", colnames(origMAF))
    colnames(origMAF) <-
      gsub(pattern = "chr", "Chromosome", colnames(origMAF))
    colnames(origMAF) <-
      gsub(pattern = "pos", "Start_Position", colnames(origMAF))
    print("Checking autofixed names")
    for (i in mafCheck) {
      if (!i %in% colnames(origMAF)) {
        if (i != "Sample_ID") {
          stop("Mismatched Columns Persist!")
        }
      }
    }
    print("Autofixing of MAF column names appears successful!")
  }
  
  #instantiate the new MAF to add our processed data to.   
  newMAF <- data.frame(matrix(ncol = 8))
  colnames(newMAF) <-
    c(
      "Sample_ID",
      "Chromosome",
      "Start_Position",
      "Reference_Allele",
      "Tumor_Seq_Allele2",
      "Origin_Node",
      "Branch_Number",
      "Binding_Nodes"
    )
  #get the normal sequence from the ancestral alignment to compare to the MAF reference
  #this isn't strictly necessary but is an easy way to make sure the aln and maf line up!!!!
  normSeq <- asrAln$seq[which(asrAln$nam == normName)][[1]]
  ###This is the main loop of the thing.
  ### because preserving ordering is helpful and important
  ### we go line by line in the original maf
  ### even though computationally, doing the dfs every time
  ### is way more intenstive. Doesn't matter for small trees
  ### for which pretty much every cancer tree will be.
  for (i in 1:nrow(origMAF)) { #for each row in the original MAF
    #get the information from that row (for readability)
    thisRow <- origMAF[i, ]
    thisSampleID <- thisRow$Sample_ID
    thisChr <- thisRow$Chromosome
    thisStart <- thisRow$Start_Position
    thisRef <- thisRow$Reference_Allele
    thisAlt <- thisRow$Tumor_Seq_Allele2
    #start DFS (depth first search)
    #lookupVec will have the actual names corresponding to every numbered node in the correct order
    lookupVec <- c(thisTree$tip.label, thisTree$node.label)
    for (j in 1:nrow(thisTree$edge)) { #for each branch/edge in the tree
      thisN1 <- thisTree$edge[j, 1] #where does the edge start? (Node)
      thisN2 <- thisTree$edge[j, 2] #where does it go? (Node)
      seq1 <- asrAln$seq[which(asrAln$nam == lookupVec[thisN1])][[1]] #get the sequence for that Node
      seq2 <- asrAln$seq[which(asrAln$nam == lookupVec[thisN2])][[1]] #above
      seq1 <- toupper(substring(seq1, i, i)) #we only want the nucleotide corresponding to the 
      seq2 <- toupper(substring(seq2, i, i)) #current row in the MAF file.
      
      #handle terminal ambiguities
      #if your alignment has Ns in it, then the rest of the code would treat them as distinct mutations
      #which we have no evidence for. So, this code ignores them, but alerts the user.
      if (seq2 == "N") {
        if (thisN2 <= thisTree$Nnode) { #this means it is a terminal/leaf node
          if(leafAmbigFlag==F){print("At least 1 Leaf Nucleotide is ambiguous (N), setting to MRCA for purposes of mutation origination");leafAmbigFlag<-T}
          #only print the warning once
          seq2 <- seq1 #a cheeky way to skip the rest of the loop
          leafAmbigCount<-leafAmbigCount+1
        }
      }
      if (seq1 == seq2) { #if they are the same, no mutation/change happened. Move to the next edge/branch!
        next()
      } #not the first time it has changed or no change
      #here is the sanity check to make sure the ref from the aln and maf line up
      if (thisRef  !=  toupper(substring(normSeq, i, i))) {
        stop(
          "Mismatch between MAF and Alignment for Normal")
      }
      #warn the user if the MAF and the sequence don't match. This could be legit
      #think A->T->C
      #so we use the alignment one, but note it for the user.
      if (thisAlt != seq2) {
        print(
          "MAF alt and aln seq don't match. This could be due to 2 sequential mutations or misaligned MAF-to-aln"
        )
        print(paste0("Position: ", as.character(i)))
        print(paste0("MAF Alt: ", thisAlt))
        print(paste0("Aln Alt: ", seq2))
        print(paste0("From: ", lookupVec[thisN1]))
        print(paste0("To: ", lookupVec[thisN2]))
        print("")
        warning("Using alignment sequence state due to MAF alt and aln seq mismatch! Tread with caution")
      }
      #if you made it this far, this is the first time the mutation was seen on the tree
      #add this result to the newMAF data structure, preserving the origMAF ordering.
      newMAF[i, ] <-
        c(paste0("B", j),thisChr,thisStart,thisRef,seq2,lookupVec[thisN2],j,
          paste0(lookupVec[thisN1], "-", lookupVec[thisN2])
        )
    }
  }
  #report back to the user the number of leaf ambiguities
  print(paste0("There were ",leafAmbigCount," leaf ambiguities."))
  #if the want to keep the extra columns, we do so here.
  if(keepExtraMAF){
    #prevent collisions
    toPort<-which(!colnames(origMAF)%in%colnames(newMAF))
    newMAF<-cbind(newMAF,origMAF[,toPort])
  }
  #if they specified a patientID, we include it here
  if(!is.null(patientID)){
    Patient_ID<-rep(x = patientID,nrow(newMAF))
    name(Patient_ID)<-"Patient_ID"
    newMAF<-cbind(Patient_ID,newMAF)
  }
  #if they want the result written to file, we do so here.
  if(!is.null(outFile)){
    write.table(x = newMAF,sep = "\t",row.names = F,col.names = T,file = outFile)
  }
  #return complete answer
  return(newMAF)
}

##Function: renderTree
###Purpose: plot a tree that maps to the ancestral MAF for interpretability.
###Mandatory Input
#### myMCC, same tree as for the ASR MAF mapping function above, newick formatted
###Optional Input
#### outputFile, will save file as png with that name. Renders in RGraphics otherwise.
###Depends
#### ape
###Output
#### a plot, either as a png or active RGraphics
renderTree<-function(myMCC,outputFile=NULL){
  library(ape)
  thisTree<-read.tree(myMCC)
  if(!is.null(outputFile)){png(filename = outputFile,width = 1000,height = 1000)}
  plot(thisTree)
  nodelabels(thisTree$node.label)
  edgelabels(paste0("B",1:(nrow(thisTree$edge))))
  if(!is.null(outputFile)){dev.off()}
}
