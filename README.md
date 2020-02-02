# acse-5-assignment-nan

# Group Member
* Sokratis Anagnostopoulos
* Lingaona Zhu
* Hao Lu

# BLAS Pre-requisites
* Have OpenBLAS installed to recogonise <cblas.h>
* To make sure it compiles, the path of installed OpenBLAS need to be exported first:
	* Below is the example used in macOS
	brew install OpenBLAS
	export LDFLAGS="-L/usr/local/opt/openblas/lib"
  	export CPPFLAGS="-I/usr/local/opt/openblas/include"
* Then in the command line compile using: 
	gcc-9 -lstdc++ -g -I/usr/local/opt/openblas/include -L/usr/local/opt/openblas/lib -lopenblas try_blas.cpp 
* If using Windows please have OpenBLAS installed first, or other open libraries that include cblas.h
