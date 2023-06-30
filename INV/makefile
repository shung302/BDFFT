FLAGS = -ffixed-line-length-none -O3 -mcmodel=medium -g
LIB	= 
OBJ	= face_lift_sub.o
BIN	=
%:	%.f
	gfortran $(FLAGS) face_lift_sub.f -c
	gfortran $(FLAGS) $(OBJ) $< -o $*
clean:	
	\rm $<
