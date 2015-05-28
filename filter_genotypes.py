import h5py

import allel

import numpy as np
from scipy.stats import describe

from crear_grafico import crear_grafico3D_angulo

def abrir_HDF5(filename):
    
    f = h5py.File(filename, 'r')

    return f



def extraer_datos_genotype(f):

    '''Extre un GenotypeArray de un fichero HDF5 cuyo npyarray '''
    ''' fue creado con la opcion calldata_2d (libreria: vcfnp)'''
    
    genotypes = f['calldata']['genotype']
    #genotypes = f['calldata']['GT']
    g = allel.model.GenotypeArray(genotypes)

    return g



########################## 

   
def cargar_filtrar_calcular(filename, max_freq=1):
    
    f = abrir_HDF5(filename)
    g = extraer_datos_genotype(f) 
    af = g.count_alleles().to_frequencies()
    is_polimorf = af<max_freq
    af_fil = (af[is_polimorf[:, 3]])
    
    return pca

#--------------------------------------------------------------------------

def filter_snps_by_maf(genotypes, max_freq):
   
    '''Elimina los SNP que tengan alguna frecuencia alelica superior a max_freq'''
    
    af = genotypes.count_alleles().to_frequencies()
    maf = np.amax(af, axis=1)
    is_polimorf = maf < max_freq	
    return genotypes[is_polimorf, :, :]

def filter_snps_by_missing_calls_freq(genotypes, num_al_snp, max_missing):
    
    '''Elimina aquellos SNP que posean hasta 'max_missing' de datos daltante''' 

    al = genotypes.count_alleles()
    suma = np.sum(al, axis=1)
    min_all = num_al_snp-max_missing
    filt = suma > min_all
      
    g_filt = genotypes[filt, :, :]
    
    return g_filt
    
  
   
    
def filter_snps_by_min_calls(genotypes, min_calls):

    '''Elimina los SNP que no tengan como minimo el numero de 'min_calls' '''

    al = genotypes.count_alleles()
    suma = np.sum(al, axis=1)
    filt = suma >= min_calls      
    g_filt = genotypes[filt, :, :]
    
    return g_filt



def filter_gn(gn, num_al_snp):

    '''Elimina aquellas filas que solo contengan el valor 2(todos hom alt)'''
    '''Si no se eliminan el PCA no se puede calcular'''
    suma = np.sum(gn, axis=1)
    div =suma/num_al_snp #coincide el numero de alelos total, con el valor para los hom alt
    filt = div<1 #si todos son 2, div sera 1.
    gn = gn[filt, :]
   
    return gn
    
def cargar_filtrar_calcular2(filename, max_freq=1):
    
    
    f = abrir_HDF5(filename)
    
    
    genotype = extraer_datos_genotype(f)
    print genotype 
    g = filter_snps_by_maf(genotype, max_freq=1)
    print g
    g_ploidy = g.ploidy
    g_samples = g.n_samples
    num_al_snp = g_ploidy*g_samples
    #filter_snps_by_miss(g, num_al_snp, max_missing=5)
    
    g = filter_snps_by_missing_calls_freq(g, num_al_snp, max_missing=100)
    
    #g = filter_snps_by_min_calls(g, min_calls=496)
    
    print g.shape
    
    g = allel.model.GenotypeArray(g)
    gn = g.to_n_alt()
    gn = filter_gn(gn, num_al_snp)
   
    

    pca = allel.stats.decomposition.pca(gn) 
    print pca
    return pca

####################################################################################################

def filtrar_gn(gn, modo=1):

    '''Elimina todas las lineas compuestas unicamente por homocigotos de un solo tipo por defecto'''
    '''En modo=2 elimina todas las que tengan mas de 0.95 de hom ref o hom alt. '''
    gn_filt=None
    
    
    for linea in gn:
        if (2 and 0 in linea) or 1 in linea:
            if modo==2:
                num = len(linea)
                suma = np.sum(linea)
                div = suma/num
                if div<0.005 or div>1.995:
                    continue
                
            if gn_filt is None:
                gn_filt=np.array(linea)
            else:
                gn_sig = np.array(linea)
                gn_filt = np.vstack((gn_filt, gn_sig))


    return gn_filt



def calcular_pca_con_hom_het(filename, modo=1):
    
    f = abrir_HDF5(filename)
    g = extraer_datos_genotype(f)
    gn = g.to_n_alt()
    gn_filt = filtrar_gn(gn, modo=modo)
    pca = allel.stats.decomposition.pca(gn_filt, n_components=3)

    return pca




################################################################
 

def separar_puntos_pca(pca):
    '''Obtiene los eje X, Y, Z como los 3 primeros componentes'''
    #Mirar si se pueden coger directamente   
    X = []
    Y = []
    Z = []
    valores=pca[0] #Los datos, 1 seria el objeto
    for i in range(len(valores)):          
        X.append(pca[0][i][0])
        Y.append(pca[0][i][1])
        Z.append(pca[0][i][2])
        #i sera el numero de cada muestra, es decir si hay 153 muestras
        #habra 153 filas con 10 componentes calculados
        #El tercer numero en el componente
        #Ej pca[0][152][1] es el valor del segundo componente de la ultima muestra.   
    
    return X, Y, Z

if __name__=='__main__':
    
    pca = cargar_filtrar_calcular2(filename)
    #pca = calcular_pca_con_hom_het(filename, modo=1)
    X, Y, Z = separar_puntos_pca(pca)
    crear_grafico3D_angulo(X, Y, Z, fhand='out.png')
    
    
