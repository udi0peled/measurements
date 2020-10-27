Measure_Name := main

App_C_Flags := -g -O0 -Wall -Wextra -Wvla -Wno-unknown-pragmas -I.
App_Cpp_Flags := $(App_C_Flags) -std=c++14
App_Link_Flags := -lcrypto

all: $(Measure_Name)

main.o: main.c
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(Measure_Name): main.o
	@$(CXX) $^ -o $@ $(App_Link_Flags)
	@echo "LINK =>  $@"

clean:
	@rm -rf $(Measure_Name) *.o