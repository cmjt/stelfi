all: custom_hawkes.so hawkes.so lgcp.so marked_lgcp.so neg_alpha_custom_hawkes.so neg_alpha_hawkes.so spatial_hawkes.so spde_hawkes.so 

custom_hawkes.so: custom_hawkes.cpp
    Rscript --vanilla -e "TMB::compile('custom_hawkes.cpp')"

hawkes.so: hawkes.cpp
    Rscript --vanilla -e "TMB::compile('hawkes.cpp')"

lgcp.so: lgcp.cpp
    Rscript --vanilla -e "TMB::compile('lgcp.cpp')"

marked_lgcp.so: marked_lgcp.cpp
    Rscript --vanilla -e "TMB::compile('marked_lgcp.cpp')"

neg_alpha_custom_hawkes.so: neg_alpha_custom_hawkes.cpp
    Rscript --vanilla -e "TMB::compile('neg_alpha_custom_hawkes.cpp')"

neg_alpha_hawkes.so: neg_alpha_hawkes.cpp
    Rscript --vanilla -e "TMB::compile('neg_alpha_hawkes.cpp')"

spatial_hawkes.so: spatial_hawkes.cpp
    Rscript --vanilla -e "TMB::compile('spatial_hawkes.cpp')"

spde_hawkes.so: spde_hawkes.cpp
    Rscript --vanilla -e "TMB::compile('spde_hawkes.cpp')"

clean:
    rm -rf *o
