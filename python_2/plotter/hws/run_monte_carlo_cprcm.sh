# !/bin/bash 


#BSUB -J run_monte_carlo_CP
#BSUB -n 144
#BSUB -R "rusage[mem=1G]"
#BSUB -q p_long 
#BSUB -o stdout%J.out
#BSUB -e stdout%J.err
#BSUB -P 0564

MODELS_nonCP=('BCCR-AUTH' 'BTU' 'CMCC' 'CNRM' 'ETHZ' 'FZJ-IDL' 'HCLIM' 'ICTP' 'KIT' 'KNMI' 'UKMO')
MODELS_CP=('BCCR-AUTH' 'BTU' 'CMCC' 'CNRM' 'ETHZ' 'FZJ-IDL' 'HCLIM' 'ICTP' 'KIT' 'KNMI' 'UKMO' 'JLU')

W=5
python do_monte_carlo_CPRCM.py $W > log.monte_carlo_cprcm

