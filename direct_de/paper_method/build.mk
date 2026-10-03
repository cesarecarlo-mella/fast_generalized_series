# make -j2 all
CXX=g++
FLAGS=-std=c++17 -O3 -march=native -ffp-contract=off -DNDEBUG -DBOOST_UBLAS_NDEBUG
GENFLAGS=-std=c++17 -O2 -march=native -ffp-contract=off -DNDEBUG
QD=-Iqdpatch -DDIFFEQS_QD
GENSRC=$(wildcard gen/fill_*.cpp gen/mv_*.cpp)
OBJ_D=$(patsubst gen/%.cpp,obj_double/%.o,$(GENSRC))
OBJ_DD=$(patsubst gen/%.cpp,obj_dd/%.o,$(GENSRC))

all: h_family_c_double h_family_c_dd h_family_gen_double h_family_gen_dd h_family_double h_family_dd frobenius_double frobenius_dd

obj_double/%.o: gen/%.cpp gen/hgen.hpp
	@mkdir -p obj_double; $(CXX) $(GENFLAGS) -DREAL=double -c -o $@ $<
obj_dd/%.o: gen/%.cpp gen/hgen.hpp
	@mkdir -p obj_dd; $(CXX) $(GENFLAGS) $(QD) -DREAL=dd_real -c -o $@ $<

h_family_gen_double: h_family_gen.cpp diffeqs.hpp gen/hgen.hpp gen/hgen_letters.hpp $(OBJ_D)
	$(CXX) $(FLAGS) -o $@ h_family_gen.cpp $(OBJ_D)
h_family_gen_dd: h_family_gen.cpp diffeqs.hpp gen/hgen.hpp gen/hgen_letters.hpp $(OBJ_DD)
	$(CXX) $(FLAGS) $(QD) -DREAL=dd_real -o $@ h_family_gen.cpp $(OBJ_DD) -lqd

h_family_double: h_family.cpp diffeqs.hpp h_data.hpp
	$(CXX) $(FLAGS) -o $@ h_family.cpp
h_family_dd: h_family.cpp diffeqs.hpp h_data.hpp
	$(CXX) $(FLAGS) $(QD) -DREAL=dd_real -o $@ h_family.cpp -lqd
frobenius_double: frobenius.cpp diffeqs.hpp h_data.hpp
	$(CXX) $(FLAGS) -o $@ frobenius.cpp
frobenius_dd: frobenius.cpp diffeqs.hpp h_data.hpp
	$(CXX) $(FLAGS) $(QD) -DREAL=dd_real -o $@ frobenius.cpp -lqd

obj_double/hgen_data.o: gen/hgen_data.cpp gen/hgen_data.hpp
	@mkdir -p obj_double; $(CXX) -std=c++17 -O1 -c -o $@ $<
h_family_c_double: h_family_c.cpp diffeqs.hpp gen/hgen_data.hpp gen/hgen_letters.hpp obj_double/hgen_data.o
	$(CXX) $(FLAGS) -o $@ h_family_c.cpp obj_double/hgen_data.o
h_family_c_dd: h_family_c.cpp diffeqs.hpp gen/hgen_data.hpp gen/hgen_letters.hpp obj_double/hgen_data.o
	$(CXX) $(FLAGS) $(QD) -DREAL=dd_real -o $@ h_family_c.cpp obj_double/hgen_data.o -lqd

h_family_c_uv_double: h_family_c.cpp diffeqs.hpp gen/hgen_data.hpp gen/hgen_letters.hpp obj_double/hgen_data.o
	$(CXX) $(FLAGS) -DUV_CHART -o $@ h_family_c.cpp obj_double/hgen_data.o
h_family_c_uv_dd: h_family_c.cpp diffeqs.hpp gen/hgen_data.hpp gen/hgen_letters.hpp obj_double/hgen_data.o
	$(CXX) $(FLAGS) $(QD) -DUV_CHART -DREAL=dd_real -o $@ h_family_c.cpp obj_double/hgen_data.o -lqd
