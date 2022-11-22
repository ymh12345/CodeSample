CC=gcc

obj=msg_init_and_fina.o msg_io.o msg_aux.o msg_main.o msg_cal_f_Verlet.o msg_time_marching_Verlet.o msg_pre_contact.o fft.o msg_RandMedia.o

target=mass_spring_global_Verlet
FLAG=-fopenmp -std=c99 -lm -O3

gsl_include=/public1/home/scb5004/gsl/include
gsl_lib=/public1/home/scb5004/gsl/lib
GSL_FLAG=-lgsl -lgslcblas -L$(gsl_lib) -I$(gsl_include)

fftw_include=/public1/home/scb5004/fftw/include
fftw_lib=/public1/home/scb5004/fftw/lib
FFTW_FLAG=-lfftw3 -I$(fftw_include) -L$(fftw_lib)

$(target): $(obj)
	$(CC) $(obj) $(FLAG) $(GSL_FLAG) $(FFTW_FLAG) -o $(target)

%.o: %.c
	$(CC) -c $<  $(FLAG) $(GSL_FLAG) $(FFTW_FLAG) -o $@

clean:
	rm -rf *.gch
	rm -rf $(target) $(obj)
install:
	cp $(target) ../bin/
	cp msg_plot_wavefield.sh ../bin/
