# makefile 
### g77 compiler offeres more info if problems!
FC = gfortran
FFLAGS = -ffixed-line-length-none -fbounds-check -g -mcmodel=medium
XDIR	 = .

# list object files
OBJS	=	azdel.o bisec.o bpt.o coordl.o \
                calc_stfgs.o calc_stfcn.o calc_stf1p.o calc_stfsac.o \
                decomq.o derivs.o derivs_pq.o dfcns.o dfint.o dgauleg.o \
	        dintw_kd.o dintw_kntbl.o driver_rk4.o driver_pq_rk4.o \
                euler_rot.o euler_subs.o eval_kntstr.o \
                fftsubs.o findr.o getfid.o \
                hessian.o intrprayp.o kernblk.o \
                locates.o pkv.o prjlag.o productv.o projectv.o \
                projq.o projqn.o \
                ray2xyz.o raylag.o           rfft.o rotate.o rtp2xyz.o \
                shootout.o signature.o \
                tracer.o tracer_pq.o \
                velodisc.o velomodel.o xyz2rtp.o \
                kernelr.o projq_r.o \
                spline.o splint.o kernblk27.o kernblk27n.o \
                kernblk27_vel.o kernblk27n_vel.o tran.o rayptbl.o

all:	kernbg27_both kernrtp kernbg27_both_vel
kernbg27_both:	kernbg27_both.f $(OBJS)
	$(FC) kernbg27_both.f $(FFLAGS) -o $(XDIR)/kernbg27_both $(OBJS) $(LIBS)

kernrtp:	kernrtp.f $(OBJS) rayptbl.o
	$(FC) kernrtp.f $(FFLAGS) -o $(XDIR)/kernrtp $(OBJS) $(LIBS)

kernbg27_both_vel:	kernbg27_both_vel.f $(OBJS)
	$(FC) kernbg27_both_vel.f $(FFLAGS) -o $(XDIR)/kernbg27_both_vel $(OBJS) $(LIBS)

*.f.o:
	$(FC) $(FFLAGS) -c $*.f

clean:
	rm -f *.o
