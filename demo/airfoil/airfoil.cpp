//cpp file to make arrays
#include <airfoil.h>



int main() {    

    Airfoil airfoil(0.1, 10, 0.3, 0.3);
    airfoil.printToVtk();
    
    return 0;
}