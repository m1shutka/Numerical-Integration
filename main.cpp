#include <iostream>
#include <iomanip>
#include <functional>
#include "Point.h"
#include "Integration_Scheme_Interval.h"

int main()
{
	//подынтегральная функция f(x) = -1/x^2
	std::function<double(const Com_Methods::Point& P)> f =
		[](const Com_Methods::Point& P) { return -1 / pow(P.x(), 2); };
	//первообразная F(x) = 1/x
	std::function<double(const Com_Methods::Point& P)> F =
		[](const Com_Methods::Point& P) { return 1 / P.x(); };

	//квадратурная формула Гаусс-1
	Com_Methods::Integration_Scheme_Interval Quadrature_Formula1(Com_Methods::Integration_Scheme::Gauss1);
	Com_Methods::Integration_Scheme_Interval Quadrature_Formula2(Com_Methods::Integration_Scheme::Gauss4);

	//начало и конец отрезка интегрирования
	auto Begin = Com_Methods::Point(1.5, 0, 0);
	auto End = Com_Methods::Point(2.1, 0, 0);

	//число сегментов
	const int Num_Segments = 1;

	//точное значение интеграла (ф. Ньютона-Лейбница)
	double I_True = F(End) - F(Begin);

	//численное значение интеграла
	double I[3];
	double Ih[3];

	for (int i = 0; i < 3; i++) {
		I[i] = Quadrature_Formula2.Calculate_Integral(Begin, End, Num_Segments * pow(2, i), f);
		Ih[i] = Quadrature_Formula2.Calculate_Integral(Begin, End, (Num_Segments * pow(2, i)) * 2, f);

		//для таблицы 1
		std::cout << i + 1<< ")********************" << std::endl;
		std::cout << std::scientific;
		std::cout << "h = " << (End.x() - Begin.x()) / (Num_Segments * pow(2, i)) << std::setw(15) << std::endl;
		std::cout << "I = " << I[i] << std::endl;
		std::cout << "|I - I_True| = " << fabs(I[i] - I_True) << std::endl;
		int k = round(log2(fabs((1 + (Ih[i] - I[i]) / (I_True - Ih[i])))));
		std::cout << "k = " << k << std::endl;
		std::cout << std::endl;
		
		//для таблицы 2
		std::cout << "I_True - I = " << I_True - I[i] << std::endl;
		std::cout << "I_True - I/I_True - Ih = " << (I_True - I[i])/(I_True - Ih[i]) << std::endl;
		std::cout << "Ih - I/2^k - 1 = " << (Ih[i] - I[i]) / (pow(2, k) - 1) << std::endl;
		double IR = Ih[i] + ((Ih[i] - I[i]) / (pow(2, k) - 1));
		std::cout << "IR = " << IR << std::endl;
		std::cout << "|I_True - IR| = " << I_True - IR << std::endl;
		std::cout << std::endl;
	}
}