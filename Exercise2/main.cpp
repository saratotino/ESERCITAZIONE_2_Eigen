#include "Eigen/Eigen"
#include <iostream>


using namespace std;
using namespace Eigen;

Vector2d decomposizionePALU(const Matrix2d& A, const Vector2d& b);
Vector2d decomposizioneQR(const Matrix2d& A, const Vector2d& b);
void checkSoluzione(const Vector2d& solPALU, const Vector2d& solQR, const Vector2d& soluzione, double& erroreRelPALU, double& erroreRelQR);

int main()
{
    Vector2d soluzione(-1.0000e+0, -1.0000e+00);

    //Primo sistema
    Matrix2d A1 = Matrix2d::Zero();
    Vector2d b1 = Vector2d::Zero();
    Vector2d x1PALU;
    Vector2d x1QR;
    double erroreRelPALU1;
    double erroreRelQR1;

    A1 << 5.547001962252291e-01, -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    b1 << -5.169911863249772e-01, 1.672384680188350e-01;

    x1PALU = decomposizionePALU(A1, b1);
    x1QR = decomposizioneQR(A1, b1);

    checkSoluzione(x1PALU, x1QR, soluzione, erroreRelPALU1, erroreRelQR1);

    cout << "Soluzione usando decomposizione PALU del primo sistema: [" << x1PALU << "] con errore realativo: " << erroreRelPALU1 << endl;
    cout << "Soluzione usando decomposizione QR del primo sistema: [" << x1QR << "] con errore realativo: " << erroreRelQR1 << endl;


    //Secondo sistema
    Matrix2d A2 = Matrix2d::Zero();
    Vector2d b2 = Vector2d::Zero();
    Vector2d x2PALU;
    Vector2d x2QR;
    double erroreRelPALU2;
    double erroreRelQR2;

    A2 << 5.547001962252291e-01, -5.540607316466765e-01,
        8.320502943378437e-01, -8.324762492991313e-01;
    b2 << -6.394645785530173e-04, 4.259549612877223e-04;

    x2PALU = decomposizionePALU(A2, b2);
    x2QR = decomposizioneQR(A2, b2);

    checkSoluzione(x2PALU, x2QR, soluzione, erroreRelPALU2, erroreRelQR2);

    cout << "Soluzione usando decomposizione PALU del secondo sistema: [" << x2PALU << "] con errore realativo: " << erroreRelPALU2 << endl;
    cout << "Soluzione usando decomposizione QR del secondo sistema: [" << x2QR << "] con errore realativo: " << erroreRelQR2 << endl;


    //Terzo sistema
    Matrix2d A3 = Matrix2d::Zero();
    Vector2d b3 = Vector2d::Zero();
    Vector2d x3PALU;
    Vector2d x3QR;
    double erroreRelPALU3;
    double erroreRelQR3;

    A3 << 5.547001962252291e-01, -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    b3 << -5.169911863249772e-01, 1.672384680188350e-01;

    x3PALU = decomposizionePALU(A3, b3);
    x3QR = decomposizioneQR(A3, b3);

    checkSoluzione(x3PALU, x3QR, soluzione, erroreRelPALU3, erroreRelQR3);

    cout << "Soluzione usando decomposizione PALU del terzo sistema: [" << x3PALU << "] con errore realativo: " << erroreRelPALU3 << endl;
    cout << "Soluzione usando decomposizione QR del terzo sistema: [" << x3QR << "] con errore realativo: " << erroreRelQR3 << endl;



    return 0;
}


//Funzione per calcolare la decomposizione PALU
Vector2d decomposizionePALU(const Matrix2d& A, const Vector2d& b)
{
    Vector2d x;
    x = A.lu().solve(b);

    return x;
}


//Funzione per calcolare la decomposizione QR
Vector2d decomposizioneQR(const Matrix2d& A, const Vector2d& b)
{
    Vector2d x;
    x = A.householderQr().solve(b);

    return x;
}


//Funzione per calcolare l'errore realtivo
void checkSoluzione(const Vector2d& solPALU, const Vector2d& solQR, const Vector2d& soluzione, double& erroreRelPALU, double& erroreRelQR)
{
    erroreRelPALU = (solPALU - soluzione).norm() / soluzione.norm();
    erroreRelQR = (solQR - soluzione).norm() / soluzione.norm();
}
