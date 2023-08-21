########################################################################################
# Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
# https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

#OBJS specifies which files to compile as part of the project
OBJS = simulate.c network.c
#CC specifies compiler 
CC = gcc 
COMPILER_FLAGS = -Wall -O2
#LINKER_FLAGS specifies libraries 
LINKER_FLAGS =  -I/usr/local/lib/ -L/usr/local/lib/ -lgsl -lgslcblas -lm
#OBJ_NAME specifies the name of exectuable
OBJ_NAME = ../simulate
#This is the target that compiles executable
all : $(OBJS)
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)
