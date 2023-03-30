'''
-----------------------------------------------------
Autor Versión 1 : Lozada Sánchez Alan Omar

Actualización a Versión 2: López García Juan Ángel
                           Morales Mendoza Fernando
                           
Actualización a Versión Server: López García Juan Ángel
                                Morales Mendoza Fernando
                           
Date : 23/06/2022
-----------------------------------------------------
Requirements :
python3.9.x
R-4.1.3
sratoolkit3.0.0
'''


import re
import gzip
import shutil
import os
import io
import GEOparse
from pysradb.sraweb import SRAweb
from subprocess import call
import logging
import warnings


#----------- Configuraciones para logging -----------#
# Mensajes de INFO, DEBUG, WARNING, ERROR, CRITICAL
logger = logging.getLogger('APP')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(fmt='[%(name)s] %(asctime)s %(levelname)s - %(message)s',datefmt='%d-%m-%y %H:%M:%S')
_handler = logging.StreamHandler()
_handler.setLevel(logging.DEBUG)
_handler.setFormatter(formatter)
logger.handlers.clear() # Borra los handlers anteriores 
logger.addHandler(_handler)


#----------- ÁREA DE TRABAJO -----------#
platform = 'GPL17225'
# Directorio de trabajo principal.
work_dir = f'/home/lopezj/LabRedesBiologicas/Schizosaccharomyces_pombe/{platform}'
os.makedirs(work_dir) if not os.path.exists(work_dir) else None # Crea la carpeta si no existe.
os.chdir (work_dir)
# Directorio para archivos soft (metadata).
metadata_dir = f'{work_dir}/metadata'
# Directorio para los experimentos (series).
experiments_dir = f'{work_dir}/experiments'
# Directorio para documentos (pdfs, csv).
documents_dir = f'{work_dir}/documents'
# Directorio del script de R (genecount.R) para contar genes.
genecount_script = '/home/lopezj/LabRedesBiologicas/Code/gene_count.R'
# Directorio del archivo de referencia del genoma (.fna.gz).
genomaref_fnafile = '/home/lopezj/LabRedesBiologicas/Schizosaccharomyces_pombe/GCF_000002945.1_ASM294v2_genomic.fna.gz'
# Directorio del archivo de anotaciones del genoma (.gtf.gz).
genomaannot_gtffile = '/home/lopezj/LabRedesBiologicas/Schizosaccharomyces_pombe/GCF_000002945.1_ASM294v2_genomic.gtf.gz'
# Forma de ejecutar un script de R.
rscript_dir = 'Rscript'


#----------- RESUMEN ÁREA DE TRABAJO -----------#
logger.info('Ejecución esblecida para Linux')
logger.info(f'Plataforma {platform}')
logger.info(f'Directorio principal {work_dir}')
logger.info(f'Directorio para metadata {metadata_dir}')
logger.info(f'Directorio para las series {experiments_dir}')
logger.info(f'Directorio para documentos {documents_dir}')
logger.info(f'Script de R para el conteo de genes {genecount_script}')
logger.info(f'Archivo de referencia del genoma {genomaref_fnafile}')
logger.info(f'Archivo de anotaciones del genoma {genomaannot_gtffile}')
logger.info(f'Directorio de R {rscript_dir}') 


#---------- DESCARGA METADATA ---------#
logger.info(f'Descarga de metadata para la plataforma {platform}')
# Descarga de metadata de la plataforma.
gpl = GEOparse.get_GEO(geo = platform, destdir = metadata_dir)
# Nombres de las series de la plataforma.
series = dict.fromkeys(gpl.metadata['series_id'], None)
logger.info(f'Descarga de metadata para la plataforma {platform} finalizada')
logger.info(f'Descarga de metadata para las series de la plataforma {platform}')
# Descarga de metadata de las series
for serie in series.keys():
    series[ serie ] = GEOparse.get_GEO(geo = serie, destdir = metadata_dir)
logger.info(f'Descarga de metadata para las series de la plataforma {platform} finalizada')


#----------- FUNCIÓN DE BORRADO DE EJECUCIONES ANTERIORES -----------#
def borrarArchvivosDeEjecucionesAnteriores(path):
    '''Borra los archivos residuales de ejecuciones anteriores del conteo de genes
    Parameters
    ----------
    path : string
        Dirección de la carpeta de una serie con los archvios residuales de ejecuciones de conteo de genes anteriores.  
    Returns
    -------
    void
    '''
    if (os.path.isdir(f'{path}/bam')):
        shutil.rmtree(f'{path}/bam')
    if (os.path.isfile(f'{path}/rnaFeatureCount.rds')):
        os.remove(f'{path}/rnaFeatureCount.rds')
    if (os.path.isfile(f'{path}/rnaFeatureCount_paired.rds')):
        os.remove(f'{path}/rnaFeatureCount_paired.rds')
    for dire in os.listdir(path):
        if (dire.startswith('my_index.')):
            os.remove(f'{path}/{dire}')
        if (dire.endswith('_trimed.fastq.gz')):
            os.remove(f'{path}/{dire}')
        if (re.fullmatch('[0-9]+\.txt|[0-9]+-[0-9]+\.txt',dire)):
            os.remove(f'{path}/{dire}')
    logger.info(f'Carpeta {path} limpiada correctamente')


#----------- CONTEO DE GENES CON SCRIPT DE R -----------#
logger.info(f'Conteo de genes para las muestras de cada serie de la plataforma {platform}')
# Expresion regular para conseguir el id GEO de archivos fastq.gz.
re_fastqfile = re.compile('(GSM[0-9]*)\.fastq\.gz')
# Conteo de genes para todas las series descargadas.
countSerie = [0, len(series.keys())]
for serie in ['GSE97982']: #series.keys():
    countSerie[0] += 1
    print(f'{"*"*(len(serie)+14)}\n****** {serie} ****** [{countSerie[0]}/{countSerie[1]}] \n{"*"*(len(serie)+14)}\n')
    path = f'{experiments_dir}/{serie}'
    # Filtra las muestras para que sean solo de la plataforma elegida.
    gsmsKeys = [x for x in series[serie].gsms.keys() if series[serie].gsms[x].metadata['platform_id'][0]==platform] 
    try:
        os.chdir(path)
        response = 'si'
        if (os.path.isfile(f'{path}/gene_counts.csv')):
            print(f'Ya existe un conteo de genes para la serie {serie}')
            response = str(input('¿Aún asi esea realizar un nuveo conteo? Teclea "Si" o "No": '))
        inexistFastqFiles = [file for file in gsmsKeys if file not in re_fastqfile.findall(" ".join(os.listdir(path)))]
        if (len(inexistFastqFiles) > 0 and response.lower() == 'si'):
            print(f'\nEn la carpeta de la serie {serie} faltan {len(inexistFastqFiles)} archivos fastq.gz de las muestras: {", ".join(inexistFastqFiles)}')
            response = str(input('¿Aún asi desea continuar con el conteo de genes de las muestras existentes? Teclea "Si" o "No": '))
        if (response.lower() == 'si'):
            logger.info(f'Ejecutando script para la serie {serie}')
            borrarArchvivosDeEjecucionesAnteriores(path)
            # Crear un txt con las muestas que hay menos las que faltan.
            with open(f'{len(gsmsKeys)}-{len(inexistFastqFiles)}.txt' if len(inexistFastqFiles)>0 else f'{len(gsmsKeys)}.txt' ,'w') as handle: 
                for sample in inexistFastqFiles:
                    handle.write(f'{sample}\n')
                logger.info(f'Documento de texto con las muestras descartadas creado correctamente ({handle.name})')
            # Ejecucion del script de R.
            print(call([rscript_dir,genecount_script,path,genomaref_fnafile,genomaannot_gtffile]))
        else:
            logger.info(f'Se ha cancelado el conteo de genes para la serie {serie}')
    except FileNotFoundError:
        logger.error(f'No existe la carpeta {path}')
logger.info('Conteo de genes finalizado')