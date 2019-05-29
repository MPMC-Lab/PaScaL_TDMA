include ../Makefile.inc

LIB = libp3ta.a

OBJDIR = ./obj

TDMAS = $(OBJDIR)/tdmas.o
RANGE = $(OBJDIR)/para_range.o
P3TA  = $(OBJDIR)/module_stdma.o

OBJS   = $(TDMAS) $(RANGE) $(P3TA)

all: $(LIB) inc

$(LIB): $(OBJS)
	ar qc $@ $(OBJS)
	cp $@ ../lib

inc: $(LIB)
	cp $(OBJDIR)/*.mod ../include

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@  -module $(OBJDIR) -I$(OBJDIR)
$(OBJDIR)/%.o:$(OBJDIR)/%.mod

%(P3TA)	 : %(RANGE) %(TDMAS)

clean:
	rm -rf ../lib ../include *__genmod.f90 $(LIB)
	rm -rf $(OBJDIR)

