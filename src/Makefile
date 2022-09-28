all: custom_hawkes.so hawkes.so lgcp.so marked_lgcp.so neg_alpha_custom_hawkes.so neg_alpha_hawkes.so spatial_hawkes.so spde_hawkes.so 

custom_hawkes.so: custom_hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('custom_hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

hawkes.so: hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

lgcp.so: lgcp.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('lgcp.cpp', safeunload=FALSE, safebounds=FALSE)"

marked_lgcp.so: marked_lgcp.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('marked_lgcp.cpp', safeunload=FALSE, safebounds=FALSE)"

neg_alpha_custom_hawkes.so: neg_alpha_custom_hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('neg_alpha_custom_hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

neg_alpha_hawkes.so: neg_alpha_hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('neg_alpha_hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

spatial_hawkes.so: spatial_hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('spatial_hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

spde_hawkes.so: spde_hawkes.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('spde_hawkes.cpp', safeunload=FALSE, safebounds=FALSE)"

clean:
	rm -rf *o
