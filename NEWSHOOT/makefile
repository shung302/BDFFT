# makefile for shootray_evt and shootray_df_evt
# compile in gluon machine
#
#FC      = /usr/local/pgi/linux86/bin/pgf77
#FC      = pgf90
#FFLAGS  =  -Mextend -O3
FC      = gfortran
FFLAGS = -ffixed-line-length-none -fbounds-check -g -mcmodel=medium -fdollar-ok
#FFLAGS  =  -g -O3

#FC	= /home/userlibs/intel9.1/fce/bin/ifort
#FC	= ifort
#FFLAGS   = -132 -O3
#LIBS	= -L/home/shung/response_stn/evalresp -levresp $$SACLIB/sacio.a
# specify libraries
PROG1    = shootray_evt
PROG2    = shootray_df_evt

# list object files
OBJS	=	azdel.o bisec.o bpt.o coordl.o \
                decomq.o derivs.o derivs_pq.o dfcns.o dfint.o dgauleg.o \
	        dintw_kd.o dintw_kntbl.o driver_rk4.o driver_pq_rk4.o \
                euler_rot.o euler_subs.o eval_kntstr.o fftsubs.o findr.o \
                getfid.o getvmdl.o getv_cmdl.o getv_smdl.o \
                hessian.o intrprayp.o locates.o \
                pkv.o prjlag.o projectv.o projqn.o \
                ray2xyz.o raylag.o rayptbl.o rfft.o rotate.o rtp2xyz.o \
                shootout.o signature.o \
                tracer.o tracer_pq.o \
                velodisc.o velomodel.o xyz2rtp.o
OBJSD	= 	raydiff.o derivs_df.o projq_df.o productv.o tracer_df.o

all:	${PROG1} ${PROG2}

$(PROG1):	shootray_evt.f $(OBJS)
	$(FC) shootray_evt.f $(FFLAGS) -o $(PROG1) $(OBJS) $(LIBS)

$(PROG2):	shootray_df_evt.f $(OBJS) $(OBJSD)
	$(FC) shootray_df_evt.f $(FFLAGS) -o $(PROG2) $(OBJS) $(OBJSD) $(LIBS)

clean:
	rm -f *.o
