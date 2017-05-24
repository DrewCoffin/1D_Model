COMP=mpiifort
FLAGS=-g -Ofast
OBJECTS=functions.o debug.o reactions.o varTypes.o global.o cfit.o rrfit.o dielectronic.o timeStep.o readEmis.o inputs.o ftmix.o para.o output.o
MODS=functions.mod debug.mod reactions.mod vartypes.mod global.mod fitc.mod fitr.mod dielectronic.mod timeStep.mod readEmis.mod inputs.mod ftmix.mod para.mod output.mod
EXE= torus

all: onebox

onebox: $(OBJECTS)
	$(COMP) $(FLAGS) -o $(EXE) onebox.f90 $(OBJECTS)

inputs.o:
	$(COMP) $(FLAGS) -c inputs.f90

timeStep.o: functions.o
	$(COMP) $(FLAGS) -c timeStep.f90

ftmix.o:
	$(COMP) $(FLAGS) -c ftmix.f90

readEmis.o: varTypes.o
	$(COMP) $(FLAGS) -c readEmis.f90

functions.o: reactions.o varTypes.o global.o cfit.o rrfit.o dielectronic.o debug.o para.o
	$(COMP) $(FLAGS) -c functions.f90

para.o:
	$(COMP) $(FLAGS) -o para.o -c ParallelVars.f90

debug.o: varTypes.o
	$(COMP) $(FLAGS) -c debug.f90

varTypes.o: 
	$(COMP) $(FLAGS) -c varTypes.f90

reactions.o:
	$(COMP) $(FLAGS) -c reactions.f90

global.o:
	$(COMP) $(FLAGS) -c global.f90

cfit.o:
	$(COMP) $(FLAGS) -c cfit.f90

rrfit.o:
	$(COMP) $(FLAGS) -c rrfit.f90

dielectronic.o: para.o
	$(COMP) $(FLAGS) -c dielectronic.f90

output.o:
	$(COMP) $(FLAGS) -c output.f90

clean:
	\rm *.o *.mod $(EXE)

