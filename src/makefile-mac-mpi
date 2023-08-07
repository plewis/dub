SRCDIR		= .
O		    = o
CC		    = /usr/local/bin/mpic++
CLINKER		= /usr/local/bin/mpic++
CFLAGS		= -std=c++11 
INCLUDES    = -I/Users/plewis/Documents/software/libraries/boost_1_73_0 -I/Users/plewis/Documents/software/libraries/ncl-2.1.18
LIBDIR      = /Users/plewis/lib/static
LIBS        = $(LIBDIR)/libboost_filesystem.a $(LIBDIR)/libboost_system.a $(LIBDIR)/libboost_program_options.a $(LIBDIR)/libncl.a
OPTF		= -O3
MATHLIB		= -lm
EXECS		= mpipversion
OBJS		= main.o

default: mpipversion

all: $(EXECS)

mpipversion: FORCE
	$(MAKE) dubmpi

dubmpi: main.o proj.hpp
	$(CLINKER) $(CFLAGS) -L$(LIBDIR) -o dubmpi $(OBJS) -lm $(LIBS)

main.o: $(SRCDIR)/main.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(OPTF) -c $(SRCDIR)/main.cpp

FORCE:
