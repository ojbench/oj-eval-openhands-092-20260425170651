#include <iostream>
#include "src.hpp"

int main(){
    const int interface_size = 3;
    const int connection_size = 3;
    int from[connection_size] = {1, 1, 2};
    int to[connection_size] = {2, 3, 3};
    fraction resistance[connection_size] = {fraction(1, 2), fraction(1, 4), fraction(2)};
    fraction current[interface_size] = {fraction(2), fraction(1), fraction(-3)};
    fraction voltage[interface_size] = {fraction(1), fraction(2), fraction(1, 2)};

    resistive_network instance(interface_size, connection_size, from, to, resistance);
    std::cout << instance.get_equivalent_resistance(1, 2) << "\n";
    std::cout << instance.get_voltage(2, current) << "\n";
    std::cout << instance.get_power(voltage) << "\n";
    return 0;
}
