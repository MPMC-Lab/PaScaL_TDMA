include ../Makefile.inc

LIB = libpascal_tdma.a

OBJDIR = ./obj

RANGE = $(OBJDIR)/para_range.o
TDMAS = $(OBJDIR)/tdmas.o
PTDMA  = $(OBJDIR)/pascal_tdma.o

TDMAS_CUDA = $(OBJDIR)/tdmas_cuda.o
PTDMA_CUDA  = $(OBJDIR)/pascal_tdma_cuda.o

OBJS   = $(RANGE) $(TDMAS) $(PTDMA)

ifdef USE_CUDA
    OBJS += $(TDMAS_CUDA) $(PTDMA_CUDA)
endif

all: $(LIB) inc

$(LIB): $(OBJS)
	ar qc $@ $(OBJS)
	cp $@ ../lib

inc: $(LIB)
	cp $(OBJDIR)/*.mod ../include

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@  -module $(OBJDIR) -I$(OBJDIR)
$(OBJDIR)/%.o: %.cuf
	$(FC) $(FLAG) -c $^ -o $@  -module $(OBJDIR) -I$(OBJDIR)
$(OBJDIR)/%.o:$(OBJDIR)/%.mod

%(PTDMA)	 : %(RANGE) %(TDMAS)

clean:
	rm -rf ../lib ../include *__genmod.f90 $(LIB)
	rm -rf $(OBJDIR)

