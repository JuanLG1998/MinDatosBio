#####################################
# CONTEO DE GENES DE ARCHIVOS FASTQ #
# Autor: Alan Omar Lozada Sánchez   #
# Date: 08/03/2022                  #
#####################################

# Install packages and dependencies
options(repos = list(CRAN="http://cran.rstudio.com/")) # Repositorio de CRAN

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!requireNamespace("ShortRead", quietly = TRUE)){
  BiocManager::install("ShortRead")
}

if (!requireNamespace("Rsubread", quietly = TRUE)){
  BiocManager::install("Rsubread")
}

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if (!requireNamespace("SeqWins", quietly = TRUE)){
  devtools::install_github("jliu678/SeqWins")
  install.packages("SeqWins")
}

# Load packages
library(SeqWins)

# Load argumentos
args = commandArgs(trailingOnly = TRUE)
data_args = list("fastq_files_path"=args[1],
                 "genomic_fna_path"=args[2],
                 "genomic_gtf_path"=args[3])

# Seleccionar los archivos fastq.gz
f1 = list.files(path = data_args$fastq_files_path, pattern = "fastq.gz")
f2 = list.files(path = data_args$fastq_files_path, pattern = "paired.fastq.gz")

# Filtro de archivos "_trimed.fastq.gz" (ejecuciones anteriores)
temp = c()
for (file in f1){
  if (!(endsWith(file,'_trimed.fastq.gz'))){
    temp = c(temp,file)
  }
}
f1 = temp

if (length(f2) > 0){
  # Muestras PAIRED
  cat('\nExisten muestras PAIRED\n')
  paired_files = c()
  for (file in f2){
    paired_files = c(paired_files,gsub(file,pattern = '_paired.fastq.gz',replacement = '.fastq.gz'))
  }
  sort(f2)
  sort(paired_files)
  single_files = c()
  for (file in f1){
    if (!(file %in% f2 || file %in% paired_files)){
      single_files = c(single_files,file)
    }
  }
  f1 = c(paired_files,single_files)
  
  cat('f1 (single and paired_1): ',length(f1),'\n')
  cat('f2 (paired): ',length(f2),'\n')
  rm(file,temp)
  
  # EJECUCION PARA LAS MUESTRAS PAIRED
  # Contar los genes. Genera archivos _trimed.fastq.gz, bam files, index y rnaFeatureCount.rds
	#Dependiendo de los recursos informaticos de su equipo, los parámetros subReadThreads y shortreadRAM deben ser modificados
	# Donde:
	# 	subReadThreads es el número de hilos que se utilizará para el procesamiento
	#	shortreadRAM es la cantidad de ram que el programa puede utilizar
	# Para ambos casos se recomienda utilizar la mitad de los recursos disponibles (mitad de hilos del procesador y mitad de RAM de la computadora) 
  seqW(fileList1=paired_files,fileList2=f2,subReadThreads=8L,shortreadRAM=3e8,genomeRefFile=data_args$genomic_fna_path, genomeAnnotFile=data_args$genomic_gtf_path, alignPairedOutput=gsub(basename(paired_files),pattern ="\\.fastq\\.gz",replacement = "\\.bam")) #RNAseq
  
  # Renombrar el archivo "rnaFeatureCount.rds" de las muestras paired
  file.rename(list.files(path = data_args$fastq_files_path,pattern='rnaFeatureCount.rds'), 'rnaFeatureCount_paired.rds')
  
  if (length(single_files)>0){
    # EJECUCION PARA LAS MUESTRAS SINGLE
    # Contar los genes. Genera archivos _trimed.fastq.gz, bam files, index y rnaFeatureCount.rds

	#Dependiendo de los recursos informaticos de su equipo, los parámetros subReadThreads y shortreadRAM deben ser modificados
	# Donde:
	# 	subReadThreads es el número de hilos que se utilizará para el procesamiento
	#	shortreadRAM es la cantidad de ram que el programa puede utilizar
	# Para ambos casos se recomienda utilizar la mitad de los recursos disponibles (mitad de hilos del procesador y mitad de RAM de la computadora) 
    seqW(fileList1=single_files,fileList2=NULL,subReadThreads=8L,shortreadRAM=3e8,genomeRefFile=data_args$genomic_fna_path, genomeAnnotFile=data_args$genomic_gtf_path) #RNAseq
    
    # Generar un archivo csv del conteo de genes de las muestras paired y single
    rnaFeatureCount_paired = data.frame(readRDS(paste(data_args$fastq_files_path,'rnaFeatureCount_paired.rds',sep='/'))$counts)
    rnaFeatureCount = data.frame(readRDS(paste(data_args$fastq_files_path,'rnaFeatureCount.rds',sep='/'))$counts)
    gene_counts = merge(rnaFeatureCount_paired,rnaFeatureCount,by='row.names')
    row.names(gene_counts) = gene_counts$Row.names
    gene_counts$Row.names = NULL
    write.csv(gene_counts,file = paste(data_args$fastq_files_path,'gene_counts.csv',sep='/'))
  }else{
    # Generar un archivo csv del conteo de genes de las muestras paired
    rnaFeatureCount_paired = readRDS(paste(data_args$fastq_files_path,'rnaFeatureCount_paired.rds',sep='/'))
    write.csv(rnaFeatureCount_paired$counts,file = paste(data_args$fastq_files_path,'gene_counts.csv',sep='/'))
  }
  
}else{
  # Muestras SINGLE
  f2 = NULL
  # Contar los genes. Genera archivos _trimed.fastq.gz, bam files, index y rnaFeatureCount.rds
	
	#Dependiendo de los recursos informaticos de su equipo, los parámetros subReadThreads y shortreadRAM deben ser modificados
	# Donde:
	# 	subReadThreads es el número de hilos que se utilizará para el procesamiento
	#	shortreadRAM es la cantidad de ram que el programa puede utilizar
	# Para ambos casos se recomienda utilizar la mitad de los recursos disponibles (mitad de hilos del procesador y mitad de RAM de la computadora) 
  seqW(fileList1=f1,fileList2=f2,subReadThreads=8L,shortreadRAM=3e8,genomeRefFile=data_args$genomic_fna_path, genomeAnnotFile=data_args$genomic_gtf_path) #RNAseq
  
  # Generar un archivo csv del conteo de genes
  rnaFeatureCount = readRDS(paste(data_args$fastq_files_path,'rnaFeatureCount.rds',sep='/'))
  write.csv(rnaFeatureCount$counts,file = paste(data_args$fastq_files_path,'gene_counts.csv',sep='/'))
}

