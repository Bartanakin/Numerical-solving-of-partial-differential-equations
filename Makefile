clean:
	rm -rf build
	mkdir build

configure:
	cmake -B build -S . -G "Unix Makefiles"

build_Elliptic:
	cmake --build build --target Elliptic --config Debug
run_Elliptic:
	./build/src/EllipticEquations/main/Elliptic.exe

build_EllipticErrors:
	cmake --build build --target EllipticErrors --config Debug
run_EllipticErrors:
	./build/src/EllipticEquations/EllipticErrors/EllipticErrors.exe

build_Parabolic:
	cmake --build build --target Parabolic --config Debug
run_Parabolic:
	./build/src/ParabolicEquations/main/Parabolic.exe

build_2Dbase:
	cmake --build build --target 2Dbase --config Debug
run_2Dbase:
	./build/src/Charts/2Dbase/2Dbase.exe $(f)
