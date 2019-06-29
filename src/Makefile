include ../Makefile.inc

LIB = libpascal_tdma.a

OBJDIR = ./obj

TDMAS = $(OBJDIR)/tdmas.o
RANGE = $(OBJDIR)/para_range.o
PTDMA  = $(OBJDIR)/pascal_tdma.o

OBJS   = $(TDMAS) $(RANGE) $(PTDMA)

all: $(LIB) inc

$(LIB): $(OBJS)
	ar qc $@ $(OBJS)
	cp $@ ../lib

inc: $(LIB)
	cp $(OBJDIR)/*.mod ../include

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@  -module $(OBJDIR) -I$(OBJDIR)
$(OBJDIR)/%.o:$(OBJDIR)/%.mod

%(PTDMA)	 : %(RANGE) %(TDMAS)

clean:
	rm -rf ../lib ../include *__genmod.f90 $(LIB)
	rm -rf $(OBJDIR)

