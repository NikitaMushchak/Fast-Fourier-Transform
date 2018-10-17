#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <vector>

#if !defined M_PI
    #define M_PI 3.14159265358979323846264338327950288
#endif

//глобальные переменные, необходимые в разных фрагментах расчётной программы

//Type
#define INITIAL -666	//Инициализация значения типа элемента

#define OUTSIDE 0
#define RIBBON  1
#define CHANNEL 2

//XDerivative:
#define LEFT   -1
#define ALONE   0
#define RIGHT   1
#define REGULAR 2

//YDerivative:
#define BOTTOM -1
#define ALONE   0
#define TOP     1
#define REGULAR 2


extern int regime;

extern double xStar0;
extern double dx;
extern double dy;
extern double bn;
extern double Q0;
extern double E;
extern double C;
extern double Kic;
extern double alpha;
extern double Amu;
extern double epsilon;
extern double dt_step;		///<Переменная для сохранения шага по времени для передачи в расчет проппанта

/// Координата источника по абсциссе
extern size_t i00;
/// Координата источника по ординате
extern size_t j00;

/// Версия сборки
extern std::string buildVersion;

/// \brief Класс граничного элемента
/// \details Класс, задающий граничный элемент (координаты – индексы ячейки) с 
/// оператором сравнения двух граничных элементов
class Ribbon{
    public:
        size_t i;
        size_t j;
        
        Ribbon(const size_t i, const size_t j): i(i), j(j){};
        
        bool operator == (const Ribbon ribbon) const{
            return this->i == ribbon.i && this->j == ribbon.j;
        }
};

/// \brief Класс ячейки сетки
/// \details Класс, задающий ячейку расчётной сетки (координаты центра, тип, 
/// тип производной вдоль обоих орт, время активации)
class Cell{
    public:
        double x = 0;
        double y = 0;
        
        int type = INITIAL;
		int xDerivative = INITIAL;
        int yDerivative = INITIAL;
        
        double activationTime = -1.;
        
        void setCoordinates(const double x, const double y){
            this->x = x;
            this->y = y;
        }
};

/// \brief Функция возвращает название режима
std::string regimeName();

/// \brief Функция выводит лого программы
void printLogo();

/// \brief Основная расчётная функция
int planar3D(
    const double T1,
    const double mu,
    const double theta,
    const double ts,
    const double cellSize,
    const double meshSize,
    const std::string pathToBarriersFile,
    const std::string pathToInjectionFile,
    const bool runningFromGUI,
    const bool saveSteps
);
