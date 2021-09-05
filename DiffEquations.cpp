#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

class EquationBase
{
public:
    EquationBase(double x0, double u0, double v0, double h, int A, int B) :
        x0(x0), u0(u0), v0(v0), h(h), A(A), B(B)
    {}

    void Runge_Kutta()
    {
        double k11, k21, k12, k22, k13, k23, k14, k24;
        double u, v;

        for (int i = A; i < B + 1; i++)
        {
            k11 = h * f1(x0, u0, v0);
            k21 = h * f2(x0, u0, v0);

            k12 = h * f1(x0 + h * 0.5, u0 + k11 * 0.5, v0 + k21 * 0.5);
            k22 = h * f2(x0 + h * 0.5, u0 + k11 * 0.5, v0 + k21 * 0.5);

            k13 = h * f1(x0 + h * 0.5, u0 + k12 * 0.5, v0 + k22 * 0.5);
            k23 = h * f2(x0 + h * 0.5, u0 + k12 * 0.5, v0 + k22 * 0.5);

            k14 = h * f1(x0 + h, u0 + k13, v0 + k23);
            k24 = h * f2(x0 + h, u0 + k13, v0 + k23);

            u0 += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
            v0 += (k21 + 2 * k22 + 2 * k23 + k24) / 6;

            x0 += h;

            Print();
        }
    }

    virtual void PrintCoeffs() = 0;

protected:
    double x0, u0, v0, h;

private:
    virtual double f1(double x, double y, double z) = 0;
    virtual double f2(double x, double y, double z) = 0;
    virtual void Print() = 0;

    int A, B;
};

class Equation_1 : public EquationBase
{
public:
    Equation_1(double x0, double u0, double v0, double h, int A, int B) :
        EquationBase(x0, u0, v0, h, A, B)
    {}

    virtual void PrintCoeffs() override
    {
        cout << "h = " << h;
        cout << " x0 = " << x0;
        cout << " y0 = " << u0;
        cout << " z0 = " << v0 << endl;

        cout << setw(2) << "x";
        cout << setw(8) << "y" << endl;
    }

private:

    virtual double f1(double x, double y, double z) override
    {
        return z;
    }

    virtual double f2(double x, double y, double z) override
    {
        double a = -z / x - (1 - (5 / (x * x))) * y;
        return a;
    }

    virtual void Print() override
    {
        cout << setw(2) << round(x0 * 1000) / 1000;
        cout << setw(8) << round(u0 * 1000) / 1000 << endl;
    }
};

class Equation_2 : public EquationBase
{
public:
    Equation_2(double x0, double y0, double z0, double h, int A, int B) :
        EquationBase(x0, y0, z0, h, A, B)
    {
        l1 = (sqrt(17) - 1) * 0.5;
        l2 = (-sqrt(17) - 1) * 0.5;
        c2 = (1 - 2 * l1) / (l2 - l1);
        c1 = 2 - c2;
    }

    virtual void PrintCoeffs() override
    {
        cout << "h = " << h;
        cout << " x0 = " << x0;
        cout << " y0 = " << u0;
        cout << " z0 = " << v0 << endl;

        cout << setw(3) << "t";
        cout << setw(10) << "u";
        cout << setw(10) << "u_exact";
        cout << setw(10) << "delta u";
        cout << setw(10) << "v";
        cout << setw(10) << "v_exact";
        cout << setw(10) << "delta v";
        cout << setw(15) << "expression #1";
        cout << setw(15) << "expression #2" << endl;
    }

private:
    virtual double f1(double x, double y, double z) override
    {
        return - 2 * y - 2 * z;
    }

    virtual double f2(double x, double y, double z) override
    {
        return - y + z;
    }

    //точное решение v
    double v_exact(double t)
    {
        return c1 * exp(l1 * t) + c2 * exp(l2 * t);
    }
    //точное v'
    double dv_exact(double t)
    {
        return c1 * l1 * exp(l1 * t) + c2 * l2 * exp(l2 * t);
    }
    //точное v'' 
    double d2v_exact(double t)
    {
        return c1 * l1 * l1 * exp(l1 * t) + c2 * l2 * l2 * exp(l2 * t);
    }
    //точное u
    double u_exact(double t)
    {
        return (v_exact(t) - dv_exact(t));
    }
    //точное u'
    double du_exact(double t)
    {
        return (dv_exact(t) - d2v_exact(t));
    }

    //первое выражение системы (для проверки аналитического решения, оно должно будет = 0)
    double eq1(double t)
    {
        return (du_exact(t) + dv_exact(t) + 3 * u_exact(t) + v_exact(t));
    }
    //второе выражение системы (для проверки аналитического решения, оно должно будет = 0)
    double eq2(double t)
    {
        return (du_exact(t) - dv_exact(t) + 3 * v_exact(t) + u_exact(t));
    }

    virtual void Print() override
    {
        cout << setw(3) << round(x0 * 1000) / 1000;
        cout << setw(10) << round(u0 * 1000) / 1000;
        cout << setw(10) << round(u_exact(x0) * 1000) / 1000;
        cout << setw(10) << round(abs(u0 - u_exact(x0)) * 1000) / 1000;
        cout << setw(10) << round(v0 * 1000) / 1000;
        cout << setw(10) << round(v_exact(x0) * 1000) / 1000;
        cout << setw(10) << round(abs(v0 - v_exact(x0)) * 1000) / 1000;
        cout << setw(15) << round(eq1(x0) * 1000) / 1000;
        cout << setw(15) << round(eq2(x0) * 1000) / 1000 << endl;
    }

    double l1, l2, c1, c2;
};

int main()
{
    // Первое уравнение и Второе уравнение (рассмотрела отрезок от 0 до 5)
    const int N = 2;
    EquationBase* equations[N] = {new Equation_1(1, 0, 1, 1, 2, 20), new Equation_2(0, 1, 2, 0.1, 0, 49)};

    for (int i = 0; i < N; i++)
    {
        equations[i]->PrintCoeffs();
        equations[i]->Runge_Kutta();
    }

    return 0;
}