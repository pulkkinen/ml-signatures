##This program was used to draw the pyrimidine context (5'and 3' neighbours) of substitutions presented in Table S2B of [Iorio et al., A Landscape of Pharmacogenomic Interactions in Cancer, Cell 166:740--754 (2016)] from the ensembl database.

##In the outmost loop, we go through the list of unique ccds id´s, and the samples with that unique ccds id in an inner loop.
##setwd('/Users/opulkki/Desktop/Iorio_data');
##Run: source("query_per_ccds_vectorized.R")

library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl");

ioriodata<-read.csv(file="IorioTableS2B.csv",header=FALSE); #Loads Iorio et al. Table S2B

data_start <- 20;
data_end <- 727247;
sample<-ioriodata[data_start:data_end,1];#Sample names
type<-ioriodata[data_start:data_end,2];#Cancer types
cDNA<-ioriodata[data_start:data_end,5];#type and location of mutation

#Make a list of unique ccds ids by stripping version indicator. This has been done and saved as ccds0.dat, so commented out 
#ccds<-ioriodata[data_start:data_end,4]; #gene ccds id's
#ccds0 <- c();
#for (i in seq(1, data_end-data_start+1)){
#    ccds_id <- strsplit(as.character(ccds[i]),split="\\."); 
#    ccds0 <- append(ccds0,ccds_id[[1]][1]);
#    }
#write.table(ccds0,'ccds0.dat');
    
rm(ioriodata);

#get ccds id's without version number and order them 
ccds0 <- read.table('ccds0.dat');
ccds0 <- ccds0[,1]
ccds_unique <- unique(as.character(ccds0))[order(unique(as.character((ccds0))))];
ccds_unique <- ccds_unique[1:18854]; #The rest are labeled by transcript seq ENST... Apprarently, no ccds for those exist
N_ccds_unique<- length(ccds_unique);

##Order the samples in numerical/alphabetical order 
ordered_samples<-unique(as.character(sample))[order(unique(as.character((sample))))];
N_samples <- length(ordered_samples); #This is the number of unique samples = 6815
    
substitution_classes_pyr <- c( "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"); #substitution of pyrimidine base
substitution_classes_pur <- c( "G>T", "G>C", "G>A", "A>T", "A>G", "A>C"); #substitution of purine base
    
#substitution class matrix (N_samples x 96 substitutions)
csc_mat <- matrix(0, nrow=N_samples, ncol=96);
    
##Convert bases to numbers A=0,C=1,G=2,T=3
alphabet <- c("A", "C", "G", "T");
base2num <- function(x){match(x, alphabet)-1};

#N_ccds_unique
for (i in seq(1, N_ccds_unique)){

    #if ( substr(as.character(ccds_unique[i]),start=1,stop=4)=="CCDS" ){ #THIS NO LONGER NEEDED BECAUSE ccds_unique contains ccds only

        sample_ind_for_this_ccds <- which(ccds0==ccds_unique[i]); #These are the original row indices of samples with this ccds
        samples_for_this_ccds <- as.character(sample[sample_ind_for_this_ccds]); #The original sample names
        dnaseq <- getSequence(id = ccds_unique[i], type="ccds", seqType="coding", mart = mart);#Use biomaRt to get the gene sequence

        if (length(dnaseq[[1]])>=1){
            dnaseq <- strsplit(dnaseq[[1]],split=""); 
            dnaseq <- dnaseq[[1]]; #Now individual bases can be called as dnaseq[i]
            
	    num_chars <- nchar(as.character(cDNA[sample_ind_for_this_ccds]));
            
            #Make a mask of exon base substitutions: 1 for those indices that describe a substitution in the coding region, zero for others (indels etc.)
            valid_substitutions <- as.numeric(grepl("-",as.character(cDNA[sample_ind_for_this_ccds]))==FALSE)*as.numeric(grepl("\\+",as.character(cDNA[sample_ind_for_this_ccds]))==FALSE)*as.numeric(substr(as.character(cDNA[sample_ind_for_this_ccds]),start=num_chars-1,stop=num_chars-1)==">");
            valid_substitutions = which(valid_substitutions == 1);

            if  ( length(valid_substitutions) >= 1 ){
                subs_classes <- substr(as.character(cDNA[sample_ind_for_this_ccds]),start=num_chars-2,stop=num_chars);#List of substitutions for ccds (e.g. ["G>A", "C>T", "elC])
                subs_classes <- subs_classes[valid_substitutions]; #only valid mutations left, e.g. strip out "elC" above
                wts <- substr(subs_classes,start=1, stop=1); #this are the wildtype bases (for checking that it matches the one in ccds)	
                SC_pur <- match(subs_classes,substitution_classes_pur ); #Make a list of substitutions in purine base. Pyrimidine base subst. are NA
                revcomps <- as.numeric(complete.cases(SC_pur)); #A mask for purine base substitutions, NA for pyrimidine
                SC_pur[is.na(SC_pur)]<-0; #Replace NA's by zeros    
                SC_pyr <- match(subs_classes,substitution_classes_pyr ); #substitutions in pyrimidine base
                SC_pyr[is.na(SC_pyr)]<-0; #Replace NA's, i.e purine base substitutions by zeros
                
                SC <- SC_pur + SC_pyr; #Merge pur and pyr cases. The vector SC contains a 
                revcomps <- revcomps[SC>0];#This cuts out any remining indels etc.
                wts<- wts[SC>0]; #Clean up the list of wildtypes
                
                valid_cdna <- as.character(cDNA[sample_ind_for_this_ccds[SC>0]]);#cDNA value (e.g. c.491G>A, c.237+1G>A) for the samples corresp. to ccds
                valid_cdna <- valid_cdna[valid_substitutions]; #Clean up
                locs <- na.omit(unlist(strsplit(unlist(valid_cdna), "[^0-9]+"))); #mutation locations, Every second element of this list is blank
                locs <- as.numeric(locs[seq(2,length(locs),2)]);#remove blanks			  

                valid_samples_for_this_ccds <- as.character(samples_for_this_ccds[SC>0]); #List of samples
                valid_samples_for_this_ccds <- valid_samples_for_this_ccds[valid_substitutions]; #Cleanup
                  
                #Some ccds sequences have changed, and the designated substitution location can be larger than the seq length.
                #We remove those locations
                valid_samples_for_this_ccds <- valid_samples_for_this_ccds[locs >= 2 & locs <= (length(dnaseq)-1)];
                wts <- wts[locs >= 2 & locs <= (length(dnaseq)-1)];
                revcomps <- revcomps[locs >= 2 & locs<= (length(dnaseq)-1)];
                SC <- SC[locs >=2 & locs<= (length(dnaseq)-1)];
                
                locs <- locs[locs >=2 & locs <= (length(dnaseq)-1)];

                #Now we get the context triplets for each substitution:
                middle_bases <- dnaseq[locs];
                contexts_5prime <- dnaseq[locs-1];
                contexts_3prime <- dnaseq[locs+1];                               

                #If any of the wt bases does not match the designated mutating base, the whole cds is discarded. This is because we do not know 
                #whether the rest of bases produce a match just by chance (1/4 prob for each).
                if ( sum(ifelse(middle_bases == wts,0,1))==0 ){

                    #Map the the current samples to the list of all, ordered samples:
                    row_indices <- match(valid_samples_for_this_ccds,ordered_samples);
                    
                    #This will map the mutation into one of the 96 context substitution classes, respecting the possible reverse complement
                    if ( (length(contexts_5prime) != length(revcomps)) || (length(contexts_3prime) != length(revcomps)) ){
                        print("warning: vector lengths don't match");
                        print(i);
                    }
                    coeff1 <- base2num(contexts_5prime)*(rep(1,length(revcomps))-revcomps) + (rep(3,length(revcomps)) - base2num(contexts_3prime))*revcomps;
                    coeff2 <- base2num(contexts_3prime)*(rep(1,length(revcomps))-revcomps) + (rep(3,length(revcomps)) - base2num(contexts_5prime))*revcomps;
                    
                    column_indices <-16*(SC-1) + 4*coeff1 +coeff2 + 1;
                    
                    row_indices <- row_indices[!is.na(column_indices)];
                    column_indices <- column_indices[!is.na(column_indices)];

                    #The following does not respect the amount of protein alterations due to alternate splicing - just the gene. So +1 is not added multiple times even if the same combination of ccds, sample and substitution occur multiple times
                    indices <- cbind(row_indices,column_indices);
                    csc_mat[indices] <- csc_mat[indices] + 1;
                    
                } else {print("warning: wt does not match middle base. This ccds will be skipped"); print(i);}
            }  
        }
    #}

    if (i%%500 == 0){write.table(csc_mat, paste0("csc_mat", i,".dat")); print(i)}
    
}
    
write.table(csc_mat, "csc-no_salience.dat")










