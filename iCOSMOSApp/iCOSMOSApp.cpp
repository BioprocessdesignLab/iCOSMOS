#include <iostream>
#include "./bin/include/iCOSMOSdllcall.h"
#include "./InputFile/equation.h"

int main()
{
    std::cout << "\nHello iCOSMOS!\n";
    iCOSMOSmain(formula);
    std::cout << "\nGoodbye iCOSMOS!\n";
    return 0;
}