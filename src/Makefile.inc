#
#   version 
#
###########################################
###########################################
#   SETUP WHERE IS THE linear algebra procedures
###########################################
###########################################
OPTIM=../OPTIMIZATION
###########################################

###########################################
###########################################
#   DONE     TENTATIVELY YOU SHOULD NOT MODIFY BELOW
###########################################
###########################################

RM = /bin/rm -f

CC =  gcc
FC =  gfortran
FL =  gfortran

F77 = gfortran
F90 = gfortran
ICC = gcc

AR = ar vr
RANLIB  = echo

#OPTIONS COMPILE
#OPTF =  -C -g -DEBUG 
#OPTC =  -C -g  -DEBUG 
#OPTL =  -C -g -DEBUG 
OPTF =  -O3 -openmp
OPTC =  -O3 -openmp
OPTL =  -O3 -openmp
#OPTFF = ${OPTF}

LOPTIM=-L $(OPTIM)/lib -lSEISCOPE_OPTIM
########################################################################
#INCPAR =  -I$(LIB)/COMMON/include
##############################################################################

INCPAR =-I $(OPTIM)/COMMON/include
LIBPAR=$(LOPTIM) 


INC = $(INCPAR)
LIB = $(LIBPAR)

