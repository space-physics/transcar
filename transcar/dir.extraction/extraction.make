INPUT=.
BUT=.
#
FLAGS=  -c -g
#
$(BUT)/extraction.out:   $(BUT)/extraction.o 
#
	f77 -O -o $(BUT)/extraction.out $(BUT)/extraction.o ${LIBGR}
#
$(BUT)/extraction.o:	$(INPUT)/extraction.f $(INCCINE) 
	cd $(BUT); f77 $(FLAGS) $(INPUT)/extraction.f
