ct: 
	clear
	nvcc --generate-line-info SOURCES/crecimientoTumoral.cu -o BIN/ct
	
ct2: 
	nvcc --generate-line-info SOURCES/crecimientoTumoralTiempos.cu -o BIN/ct

ctDLB: 
	clear
	nvcc --generate-line-info SOURCES/crecimientoTumoralDLB.cu -o BIN/ct
	
ctDLBbits:
	nvcc --generate-line-info SOURCES/crecimientoTumoralDLBbits.cu -o BIN/ct
	
ctbits:
	nvcc --generate-line-info SOURCES/crecimientoTumoralbits.cu -o BIN/ct
	
Gom:
	nvcc --generate-line-info SOURCES/Gompertz.cu -o BIN/ct	
	
rm:
	rm -f TABLES/*	
	rm -f IMAGES/*	
	rm -f OUTPUT/*
	rm -f TIMES/*

exe:
	./BIN/ct
	python3 SOURCES/grafiquear.py	
	
run: ct rm exe

runTiempos: ct2 rm exe
	
runDLB: ctDLB rm exe

runDLBbits: ctDLBbits rm exe

runbits: ctbits rm exe
	
runGom: Gom rm exe
	
clean: rm
	rm -f BIN/*
