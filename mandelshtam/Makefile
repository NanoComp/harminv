FFLAGS = -O2 -ff90 # g77 flags

OUTPUT = parameters ReSp

all: $(OUTPUT)

generate: generate.f
	$(F77) $(FFLAGS) generate.f -o $@

signal param: generate
	./generate

fd_qz: fd_qz99.f
	$(F77) $(FFLAGS) fd_qz99.f -o $@

$(OUTPUT): signal param fd_qz in_qz
	./fd_qz < in_qz

signal.dat: signal
	sed 's/  */ /g' signal | tail +2 | cut -d" " -f3-4 | sed 's/ /+/;s/$$/i/' > $@

signal.mat: signal
	(echo "# Created from signal file"; \
	 echo "# name: s"; \
	 echo "# type: complex matrix"; \
	 echo "# rows: "`tail +2 signal | wc -l signal`; \
	 echo "# columns: 1"; \
	 sed 's/  */ /g' signal | tail +2 | cut -d" " -f3-4 | sed 's/^/(/;s/ /,/;s/$$/)/') > $@


clean:
	rm -f signal param fd_qz generate $(OUTPUT) signal.dat signal.mat
