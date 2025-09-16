CC       = g++
CFLAGS   = 
INCLUDE  = /usr/local/include
MY_APP	 = san
LIB	 = /usr/local/lib

run:	mouse0D.c
	${CC} ${CFLAGS} -c mouse0D.c -o mouse0D.o
	${CC} ${CFLAGS} mouse0D.o -lm -lsundials_cvode -lsundials_nvecserial -lblas -llapack -o ${MY_APP} 
	
clean: 
	rm -f *.o san *.dat *.out *.vtk
	
1d:		mouse1d.c
	${CC} ${CFLAGS} -c mouse1d.c -o mouse1d.o -fopenmp
	${CC} ${CFLAGS} mouse1d.o -lm -fopenmp -lsundials_cvode -lsundials_nvecserial -lblas -llapack -o ${MY_APP}
	
2d:		mouse2d.c
	${CC} ${CFLAGS} -c mouse2d.c -o mouse2d.o -fopenmp
	${CC} ${CFLAGS} mouse2d.o -lm -fopenmp -lsundials_cvode -lsundials_nvecserial -lblas -llapack -o ${MY_APP}



2d_op:		mouse2d_operator_splitting.c
	${CC} ${CFLAGS} -c mouse2d_operator_splitting.c -o mouse2d_operator_splitting.o -fopenmp
	${CC} ${CFLAGS} mouse2d_operator_splitting.o -lm -fopenmp -lsundials_cvode -lsundials_nvecserial -lblas -llapack -o ${MY_APP}
