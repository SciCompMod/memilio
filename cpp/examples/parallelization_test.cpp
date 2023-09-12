#include <iostream>
#include <vector>
#include <omp.h>

int main() {
    const int vectorSize = 1000000; // Größe der Vektoren anpassen

    std::vector<double> vectorA(vectorSize);
    std::vector<double> vectorB(vectorSize);
    std::vector<double> vectorC(vectorSize);
    std::vector<double> vectorD(vectorSize);

    // Zufällige Werte für die Vektoren generieren
    #pragma omp parallel for
    for (int i = 0; i < vectorSize; ++i) {
        vectorA[i] = static_cast<double>(rand()) / RAND_MAX;
        vectorB[i] = static_cast<double>(rand()) / RAND_MAX;
        vectorC[i] = static_cast<double>(rand()) / RAND_MAX;
        vectorD[i] = 0;
    }

    // Vektortriade parallel berechnen
    #pragma omp parallel for
    for (int i = 0; i < vectorSize; ++i) {
        vectorD[i] = vectorA[i] + vectorB[i] - vectorC[i];
    }

    std::cout << "Vektortriade abgeschlossen." << std::endl;

    return 0;
}
