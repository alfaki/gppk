SYSTEM    = x86-64_rhel4.0_3.4
LIBFORMAT = static_pic

CCC = g++

#-------------------------------------------------------------------------------
# compiler options
#-------------------------------------------------------------------------------

CCOPT = -O -fPIC -fexceptions -DIL_STD -g

CPLEXDIR   = $(ILOG_HOME)/cplex102
CONCERTDIR = $(ILOG_HOME)/concert24

GMSLIBDIR     = $(GAMS_HOME)/apifiles
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

GMSINCDIR     = $(GMSLIBDIR)/gdx
CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert \
            -lm -lpthread
CCFLAGS   = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
GMSFLAGS  = -I$(GMSINCDIR) -I$(GMSLIBDIR)/common -ldl

SRCDIR = src
INCDIR = include
DATDIR = dpp

GPP_PRJ = gppbnd.o gppgraph.o gpplib.o gppmain.o gppmodels.o

PRJ_NAME = gppk

#-------------------------------------------------------------------------------
# compile, execute, clean files.
#-------------------------------------------------------------------------------

all: $(PRJ_NAME) printgppk

execute: $(PRJ_NAME)
	./$(PRJ_NAME) dsc $(DATDIR)/RT2.dpp
	
clean :
	/bin/rm -rf *.o $(PRJ_NAME) dscgppk printgppk gms/lst/*.lst
	/bin/rm -rf 225* lps/*.lp lps/*.log *.lst *.log *.opt out/*.out
	/bin/rm -rf *.tab *.gdx

# ------------------------------------------------------------------------------
# compilations
# ------------------------------------------------------------------------------
$(PRJ_NAME): $(GPP_PRJ)
	$(CCC) $(CCFLAGS) $(GMSFLAGS) $(GPP_PRJ) $(GMSINCDIR)/gdxco.cpp $(GMSINCDIR)/gdxcc.c -o $(PRJ_NAME) $(CCLNFLAGS)

gppbnd.o: $(SRCDIR)/gppbnd.cpp $(SRCDIR)/$(PRJ_NAME).h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/gppbnd.cpp -o gppbnd.o
	
gppgraph.o: $(SRCDIR)/gppgraph.cpp $(SRCDIR)/$(PRJ_NAME).h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/gppgraph.cpp -o gppgraph.o
	
gpplib.o: $(SRCDIR)/gpplib.cpp $(SRCDIR)/$(PRJ_NAME).h
	 $(CCC) -c $(CCFLAGS) $(SRCDIR)/gpplib.cpp -o gpplib.o

gppmain.o: $(SRCDIR)/gppmain.cpp $(SRCDIR)/$(PRJ_NAME).h
	 $(CCC) -c $(CCFLAGS) $(SRCDIR)/gppmain.cpp -o gppmain.o
	 
gppmodels.o: $(SRCDIR)/gppmodels.cpp $(SRCDIR)/$(PRJ_NAME).h
	$(CCC) -c $(CCFLAGS) $(GMSFLAGS) $(SRCDIR)/gppmodels.cpp -o gppmodels.o
	
printgppk: printgppk.o
	$(CCC) printgppk.o -o printgppk
printgppk.o: $(SRCDIR)/printgppk.cpp
	$(CCC) -c $(SRCDIR)/printgppk.cpp -o printgppk.o
# ------------------------------------------------------------------------------
# eof
# ------------------------------------------------------------------------------
