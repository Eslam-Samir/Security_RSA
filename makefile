CPPFLAGS = -std=c++11 -o
CXX = g++
	
RSA: RSA.cpp
	$(CXX) -O2 $? $(CPPFLAGS) RSA.o
	
clean: 
	rm -f *~ 
	rm -f *.o
