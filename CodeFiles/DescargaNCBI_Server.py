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

#----------- BIBLIOTECAS -----------#
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
# Plataforma sobre la que se va a trabajar.
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
# Directorio configurado en SRA Toolkit para almacenamiento de archivos Prefetch
tempSRA = '/home/lopezj/LabRedesBiologicas/Code/sratoolkit.3.0.0-ubuntu64/temp/sra'


#----------- RESUMEN ÁREA DE TRABAJO -----------#
logger.info('Ejecución esblecida para Linux')
logger.info(f'Plataforma {platform}')
logger.info(f'Directorio principal {work_dir}')
logger.info(f'Directorio para metadata {metadata_dir}')
logger.info(f'Directorio para las series {experiments_dir}')
logger.info(f'Directorio para documentos {documents_dir}')


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



#---------- DESCARGA SRATOOLKIT ---------#
db = SRAweb()
# Ejecutar hasta que no haya errores de descarga o se considere finalizado.
logger.info('Descarga de archvios sra y transformacion a fastq.gz')
tryAgain = 'si'
while (tryAgain == 'si'):
    errores = {}
    countSerie = [0,len(series.keys())]
    for serie in ['GSE97982']: #series.keys():
        countSerie[0] += 1
        print(f'{"*"*(len(serie)+14)}\n****** {serie} ****** [{countSerie[0]}/{countSerie[1]}] \n{"*"*(len(serie)+14)}\n')
        # Control de errores de descarga.
        errores[serie] = []
        # Filtra las muestras para que sean solo de la plataforma elegida.
        gsmsKeys = [x for x in series[serie].gsms.keys() if series[serie].gsms[x].metadata['platform_id'][0]==platform]
        countMuestra = [0,len(gsmsKeys)]
        for muestra in gsmsKeys:
            countMuestra[0] += 1
            print(f'****** {muestra} ****** [{countMuestra[0]}/{countMuestra[1]}]')
            path = f'{experiments_dir}/{serie}'
            if os.path.isfile(f'{path}/{muestra}.fastq.gz'):
                print(f'El archivo {muestra}.fastq.gz ya existe\n')
            else:
                has_error = True
                num_error = 0
                while (has_error and num_error < 10):
                    try:
                        # Obtiene la metadata con los id SRR.
                        sampleMetadata = db.sra_metadata(muestra)
                        pathTemp = f'{experiments_dir}/{serie}/{muestra}'
                        os.makedirs(pathTemp) if not os.path.exists(pathTemp) else None # Crea la carpeta si no existe.
                        for srr in sampleMetadata['run_accession']:
                    
                            # Descargar el archivo sra.
                            logger.info(f'Descarga sra: {srr}')
                            call(['prefetch','-p',srr])
                            logger.info(f'Descarga finalizada: {srr}')
                            
                            # Transformacion a fastq
                            logger.info('Transformación a fastq iniciada.')
                            call(['fasterq-dump',srr,'-p','-O',pathTemp])
                            logger.info('Transformacion a fastq finalizada.')
                            
                        # compresión gzip y renombre de archivos fastq
                        logger.info('Compresión gz iniciada.')
                        fastqFiles = [f'{pathTemp}/{x}' for x in sorted(os.listdir(pathTemp))]
                        index = 0
                        for file in fastqFiles:
                            logger.info(f'Archivo {file} comprimiendose...')
                            if index == 0:
                                with open(file, 'rb') as f_in:
                                    with gzip.open(f'{path}/{muestra}.fastq.gz', 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                index += 1
                            else:
                                if (file.endswith('_2.fastq')):
                                    with open(file, 'rb') as f_in:
                                        with gzip.open(f'{path}/{muestra}{f"_{index}" if index>1 else ""}_paired.fastq.gz', 'wb') as f_out:
                                            shutil.copyfileobj(f_in, f_out)
                                else:
                                    index += 1
                                    with open(file, 'rb') as f_in:
                                        with gzip.open(f'{path}/{muestra}_{index}.fastq.gz', 'wb') as f_out:
                                            shutil.copyfileobj(f_in, f_out)
                        logger.info('Compresión gz finalizada.')            
                        
                        #Eliminando carpetas de archivos innecesarios
                        shutil.rmtree(tempSRA)
                        shutil.rmtree(pathTemp)
                        
                        has_error = False
                    except Exception as e:
                        logger.error(f'{e}.  Intento {num_error+1}')
                        has_error = True
                        num_error += 1
                if has_error:
                    logger.error(f'Ocurrio un error de descarga para la muestra {muestra} :c')
                    errores[serie].append(muestra)

    # ERRORES.
    for serie in errores.keys():
        if len(errores[serie]) > 0:
            logger.warning(f'La serie {serie} tuvo {len(errores[serie])} errores: {errores[serie]}')
        else:
            logger.info(f'La serie {serie} tuvo {len(errores[serie])} errores: {errores[serie]}')
    # Si hay errores de descarga se pregunta si desea repetir las descargas.
    if sum([len(errores[x]) for x in errores.keys()]) > 0:
        tryAgain = str(input('¿Desea ejecutar nuevamente? Teclea "Si" o "No": ')).lower()
    else:
        tryAgain = 'no'
logger.info('Descarga de archvios sra y transformacion a fastq.gz finalizada')