all: CImg.h p3.cpp
	g++ counter.cpp -o counter -lX11 -lpthread

clean:
	rm counter