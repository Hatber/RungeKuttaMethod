#include <iostream> //Консольный вывод
#include <math.h>   //Мат. функции: возведение в степень, корни и т.д.

#include <fstream> //Запись в файл

#include "RungeKutta.hpp"

using namespace std;
using RK::MersonMethod;


const double pi = 3.14159265359;
const double g = 9.80665;

const double TreeCrownsHeight = 10;

const double PressureEnviroment = 101325.0;   //В паскалях
const double AirUniversalGasConstant = 8.314/(6.022*pow((double)(10.0),(double)(-4.0))*(46.5+3.3+53.2)/3.0);

const double alpha = 0.6; //Коэффициент вовлечения воздуха из приземного слоя атмосферы в термик
const double betta = 0.6; //Коэффициент вовлечения продуктов горения в термик (0..1);

double BurningRate = 1; //скорость горения лесных горючих материалов (ЛГМ) [кг/с];

double TemperatureEnviroment = 288;

//*** Next is Magic ***//
double Wg = 1;
double Tg = 293;
double cd = 1;
double cp = 1;

double M; //Масса термика
double W; //Скорость термика
double T = 323; //Температура термика
double Z; //Высота подъема термика


double ThermalMass(double CurrentStep, double CurrentValue) {
    return alpha*CurrentValue*W + betta*BurningRate;
}

double ThermalSpeed(double CurrentStep, double CurrentValue) {
    double Ro = PressureEnviroment/(AirUniversalGasConstant*T);
    return ((T - TemperatureEnviroment)*g)/TemperatureEnviroment +
           (betta*BurningRate*(Wg-CurrentValue))/M - alpha*CurrentValue*CurrentValue -
           ((pi*cd)/(2*M))*pow((3*M)/(4*pi), 2./3) * pow(Ro, 1./3) * CurrentValue*CurrentValue;
}

double ThermalTemperature(double CurrentStep, double CurrentValue) {
    return ((BurningRate*betta)/M)*(Tg + TemperatureEnviroment - 2*CurrentValue) +
           (TemperatureEnviroment - CurrentValue)*alpha*W - (g*W*CurrentValue)/(cp*TemperatureEnviroment);
}

double ThermalLift(double CurrentStep, double CurrentValue) {
    return CurrentValue + W*CurrentStep;
}

typedef MersonMethod RungeType;
//typedef EulerMethod RungeType;

int main() {
    const double DeltaX = 0.0001;  //Шаг по координате
    const double StartX = 0;
    const double DestinationX = 20;
    const double Epsilon = 0.0001;

    //** Задание начальных условий **//
    M = (pi*PressureEnviroment*(TreeCrownsHeight/2))/(6*AirUniversalGasConstant*323)*pow(TreeCrownsHeight, 3);
    cout << "Start Thermal Mass: " << M << endl;

    W = sqrt((g*pow(TreeCrownsHeight, 2)*(323-288))/288);
    cout << "Start Thermal Speed: " << W << endl;

    T = 323;
    cout << "Start Thermal Temperature: " << T << endl;

    Z = 0;
    cout << "Start Thermal Lift: " << Z << endl;

    //** Сюда записываются результаты вычислений **//
    ofstream MASS("MASS.txt");
    ofstream SPEED("SPEED.txt");
    ofstream TEMPERATURE("TEMPERATURE.txt");
    ofstream MOVEMENT("MOVEMENT.txt");
    ofstream ALL("ALL.txt");

    //** Всего у нас 3 дифура **//
    RungeType MM(StartX, DestinationX, DeltaX, Epsilon, M, ThermalMass);
    RungeType MS(StartX, DestinationX, DeltaX, Epsilon, W, ThermalSpeed);
    RungeType MT(StartX, DestinationX, DeltaX, Epsilon, T, ThermalTemperature);
    RungeType ML(StartX, DestinationX, DeltaX, Epsilon, Z, ThermalLift);

    int endStatus;
    bool Recalc[3] = {false};
    double RecommendDeltaStep[3];
    double NewDeltaStep;

    size_t StepCounter = 0;
    size_t AllStepCounter = 0;

    //*** Непосредственно счет **//
    while(!endStatus) {
        endStatus = 0;
        AllStepCounter++;

        endStatus |= MM.MakeStep();
        endStatus |= MS.MakeStep();
        endStatus |= MT.MakeStep();
        ML.MakeStep();

        RecommendDeltaStep[0] = MM.CalculatingDeltaStep(Recalc[0]);
        RecommendDeltaStep[1] = MS.CalculatingDeltaStep(Recalc[1]);
        RecommendDeltaStep[2] = MT.CalculatingDeltaStep(Recalc[2]);

        NewDeltaStep = min(RecommendDeltaStep[0], min(RecommendDeltaStep[1], RecommendDeltaStep[2]));

        MM.SetDeltaStep(NewDeltaStep);
        MS.SetDeltaStep(NewDeltaStep);
        MT.SetDeltaStep(NewDeltaStep);

        //*** Необходимо, чтобы выполнялась заданная точность ***//
        if(Recalc[0] || Recalc[1] || Recalc[2]) {
                cout << "FUS!" << endl;
            continue;
        }

        MM.ConfirmStep();
        MS.ConfirmStep();
        MT.ConfirmStep();
        ML.ConfirmStep();

        M = MM.GetCurrentValue();
        W = MS.GetCurrentValue();
        T = MT.GetCurrentValue();
        Z = ML.GetCurrentValue();
        if(Z<=0) {
            break;
        }

        ALL << M << "\t" << W << "\t" << T << "\t" << Z << endl;

        MASS         << M << " " << MM.GetCurrentStep() << endl;
        SPEED        << W << " " << MS.GetCurrentStep() << endl;
        TEMPERATURE  << T << " " << MT.GetCurrentStep() << endl;
        MOVEMENT     << Z << " " << ML.GetCurrentStep() << endl;
        StepCounter++;
    }

    cout << endl;
    cout << "End Thermal Mass: " << M << endl;
    cout << "End Thermal Speed: " << W << endl;
    cout << "End Thermal Temperature: " << T << endl;
    cout << "End Thermal Lift: " << Z << endl;

    cout << endl << "StepCounter: " << StepCounter << "\t" << "AllStepCounter: " << AllStepCounter << endl;

    MASS.close();
    SPEED.close();
    TEMPERATURE.close();
    MOVEMENT.close();

    return 0;
}
