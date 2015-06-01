
import h5py
import allel
import numpy as np


def open_HDF5(filename):
    h5_file = h5py.File(filename, 'r')
    return h5_file


def extract_genotype_data(h5_file):

    '''Transform the genotype data in h5_file to a GenotypeArray'''
     
  
    #genotypes = h5_file['calldata']['genotype']
    genotypes = h5_file['calldata']['GT']
    return allel.model.GenotypeArray(genotypes)

    


def filter_snps_by_maf(genotype_array, max_freq=1):
   
    '''Eliminates the SNP with a allelic frequencie bigger than max_freq'''
    
    af = genotype_array.count_alleles().to_frequencies()
    maf = np.amax(af, axis=1)
    is_polimorf = maf < max_freq	
    return genotype_array[is_polimorf, :, :]


def calculate_total_alleles(genotype_array):
 
    '''With the genotype array, guess the ploidy and the number of samples to calculate'''
    '''how many alleles per SNP have the array'''

    g_ploidy = genotype_array.ploidy
    g_samples = genotype_array.n_samples
    return g_ploidy*g_samples


def filter_snps_by_missing_calls(genotype_array, num_al_snp, max_missing):
    
    '''Eliminate SNP whith less missing data than max_missing'''
    '''num_al_snp is the total number of alleles per SNP''' 

    al = genotype_array.count_alleles()
    suma = np.sum(al, axis=1)
    min_all = num_al_snp-max_missing
    filt = suma > min_all
      
    return genotype_array[filt, :, :]
    

def filter_snps_by_min_calls(genotype_array, min_calls):
    
    '''Eliminates SNP whith less data than min_calls'''

    al = genotype_array.count_alleles()
    suma = np.sum(al, axis=1)
    filt = suma >= min_calls      
    return genotype_array[filt, :, :]



def filter_gn(gn, num_al_snp):
 
    '''Eliminates the rows with only 2'''
    '''This filter is requiered if you have any row whith all 2, because PCA falls'''
    '''Works with ploidy==2'''
    suma = np.sum(gn, axis=1)
    div =suma/num_al_snp 
    filt = div<1 #If all are 2, div will be 1.
    return gn[filt, :]
   

    
def filter_(genotype_array, filters=True, num_al_snp=None, max_missing=0, min_calls=1):

    '''Uses at least maf filter'''
    '''If filters==2 filts with missing call, 3 with min calls, and 4 with both'''
        
    if filters:
        g_filt = filter_snps_by_maf(genotype_array)
    if filters==2:
       g_filt = filter_snps_by_missing_calls(g_filt, num_al_snp, max_missing) 
    elif filters==3:
       g_filt = filter_snps_by_min_calls(g_filt, min_calls)
    elif filters==4:
       g_filt = filter_snps_by_missing_calls(g_filt, num_al_snp, max_missing)
       g_filt = filter_snps_by_min_calls(g_filt, min_calls)
    
    return g_filt
    




def obtain_pca_points(pca):
    '''Obtain axes X, Y, Z like the first, second and third components'''
      
    X = []
    Y = []
    Z = []
    valores=pca[0] #pca[0] are point data, pca[1] is the object
    for i in range(len(valores)):          
        X.append(pca[0][i][0])
        Y.append(pca[0][i][1])
        Z.append(pca[0][i][2])
        #i is the sample number
        #The third number is the component
           
    
    return X, Y, Z


def main(filename=None, filters=None, max_freq=1):
    
    
   
    h5_file = open_HDF5(filename)
    genotype_array = extract_genotype_data(h5_file)
    num_al_snp = calculate_total_alleles(genotype_array)
    g_filt = filter_(genotype_array, filters=2, num_al_snp=num_al_snp, max_missing=500)
    gn = g_filt.to_n_alt()
    gn = filter_gn(gn, num_al_snp)
    
    #Deleting the variables with genotypes that we don't need, increases speed and
    #the max number of snp that can calculate the pca function

    del(genotype_array)
    del(g_filt)
    
    pca = allel.stats.decomposition.pca(gn) 
    
    X, Y, Z = obtain_pca_points(pca)

    return X, Y, Z


if __name__=='__main__':

    main()    
    

