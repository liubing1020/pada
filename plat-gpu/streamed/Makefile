bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu: ${METHOD}/src/strbal/*.java ${METHOD}/declarations.cu ${METHOD}/orchestration.cu config/${MODEL}_${ARCH}_${REAL}.${METHOD}
	javac ${METHOD}/src/strbal/*.java
	rm -f bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu
	java -ea -cp ${METHOD}/src strbal.Algorithm models/${MODEL}.model `grep nrWarps config/${MODEL}_${ARCH}_${REAL}.${METHOD} | awk '{print $$3}'` bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation_auto.cu  ${REAL} `grep nrAccess config/${MODEL}_${ARCH}_${REAL}.${METHOD} | awk '{print $$3}'` `grep fullThreshold config/${MODEL}_${ARCH}_${REAL}.${METHOD} | awk '{print $$3}'`
	cp -f ${METHOD}/declarations.cu bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu
	cat bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation_auto.cu >> bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu
	cat ${METHOD}/orchestration.cu >> bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu
