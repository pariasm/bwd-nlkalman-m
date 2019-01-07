MEX     = mex                # for matlab
# MEX     = mkoctfile --mex    # for octave

LDLIBS  = -lpng -ltiff -ljpeg
CFLAGS  = -fPIC
MXOPTS  = -O CFLAGS="$(CFLAGS) -std=c99"

BIN = iio_read.mexa64 iio_write.mexa64 # for matlab
# BIN = iio_read.mex iio_write.mex # for octave

default : $(BIN)

%.mexa64: %.c iio.o ; $(MEX) $(MXOPTS) $< iio.o $(LDLIBS)

clean   :           ; $(RM) *.o *.mex *.mexa64

test    : $(BIN)    ; octave example_iio.m
