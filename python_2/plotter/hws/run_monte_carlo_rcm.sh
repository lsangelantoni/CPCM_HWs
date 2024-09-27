# !/bin/bash 


#BSUB -J run_monte_carlo_nonCP
#BSUB -n 20
#BSUB -R "rusage[mem=4G]"
#BSUB -q p_long 
#BSUB -o stdout%J.out
#BSUB -e stdout%J.err
#BSUB -P 0564

MODELS_nonCP=('BCCR-AUTH' 'BTU' 'CMCC' 'CNRM' 'ETHZ' 'FZJ-IDL' 'HCLIM' 'ICTP' 'KIT' 'KNMI' 'UKMO')
MODELS_CP=('BCCR-AUTH' 'BTU' 'CMCC' 'CNRM' 'ETHZ' 'FZJ-IDL' 'HCLIM' 'ICTP' 'KIT' 'KNMI' 'UKMO' 'JLU')

python do_monte_carlo_RCM.py > log.monte_carlo_rcm

