#include <iostream> //Консольный вывод
#include <math.h>   //Мат. функции: возведение в степень, корни и т.д.

#include <fstream> //Запись в файл

#include "RungeKutta.hpp"
#include "Out/GNUPlot/GNUPlotOut.hpp"

using namespace std;
using RK::MersonMethod;


const double pi = 3.14159265359;
const double g = 9.80665;

const double TreeCrownsHeight = 10;

const double PressureEnviroment = 101325.0;   //В паскалях
const double AirUniversalGasConstant = 8.314/(6.022*pow((double)(10.0),(double)(-4.0))*(46.5+3.3+53.2)/3.0);

const double alpha = 0.02; //Коэффициент вовлечения воздуха из приземного слоя атмосферы в термик
const double betta = 0.5; //Коэффициент вовлечения продуктов горения в термик (0..1);

double BurningRate = 60; //скорость горения лесных горючих материалов (ЛГМ) [кг/с];

double TemperatureEnviroment = 288;
double ThermalStartTemperature = 400;

//*** Next is Magic ***//
double Wg = 7;
double Tg = 1200;
double cd = 1; //коэффициент сопротивления?
double cp = 1; //еплоемкость термика?

double M; //Масса термика
double W; //Скорость термика
double T = ThermalStartTemperature; //Температура термика
double Z; //Высота подъема термика

double Stratification() {
    double Te;

//    if(Z<1) { Te=1200; }
//    else if(Z>=2&&Z<3) { Te=800; }
//    else if(Z>=3&&Z<6) { Te=600; }
//    else if(Z>=6&&Z<9) { Te=400; }
//    else { Te=273; }

    Te = TemperatureEnviroment; - Z * 0.06;
    //if(Te < 283) { Te = 283; }

    return Te;//-TemperatureEnviroment/44308.0*Z;;
}


double ThermalMass(double CurrentStep, double CurrentValue) {
    return alpha*CurrentValue*fabs(W) + betta*BurningRate;
}

double ThermalSpeed(double CurrentStep, double CurrentValue) {
    double Ro = PressureEnviroment/(AirUniversalGasConstant*T);
    return
            //1/(alpha*CurrentValue + 1) *
            //(
                ((T - Stratification())*g)/Stratification() //?
                +(betta*BurningRate*(Wg-CurrentValue))/M    //Влет продуктов горения
                - alpha*CurrentValue*CurrentValue
                - ((pi*cd)/(2*M))*pow((3*M)/(4*pi), 2./3) * pow(Ro, 1./3) * CurrentValue*CurrentValue;
           // );
}

double ThermalTemperature(double CurrentStep, double CurrentValue) {
    return ((BurningRate*betta)/M)*(Tg + Stratification() - 2*CurrentValue) //Вдув от пожара
            + (Stratification() - CurrentValue)*alpha*W //Приток воздуха
            - (g*W*CurrentValue)/(cp*Stratification());
}

double ThermalLift(double CurrentStep, double CurrentValue) {
    return W;
}

typedef MersonMethod RungeType;
//typedef EulerMethod RungeType;

int main() {
    const double DeltaX = 0.0001;  //Шаг по координате
    const double StartX = 0;
    const double DestinationX = 60*5;
    const double Epsilon = 0.0001;

    const double writeStep = 0.1;
    double lastWrite = -writeStep;

    //** Задание начальных условий **//
    M = (pi*PressureEnviroment*(TreeCrownsHeight/2))/(6*AirUniversalGasConstant*ThermalStartTemperature)*pow(TreeCrownsHeight, 3);
    cout << "Start Thermal Mass: " << M << endl;

    W = sqrt((g*pow(TreeCrownsHeight, 2)*(ThermalStartTemperature-TemperatureEnviroment))/TemperatureEnviroment);
    cout << "Start Thermal Speed: " << W << endl;

    T = ThermalStartTemperature;
    cout << "Start Thermal Temperature: " << T << endl;

    Z = 0;
    cout << "Start Thermal Lift: " << Z << endl;

    //** Сюда записываются результаты вычислений **//

    GNUPlotOut MASS;
    GNUPlotOut SPEED;
    GNUPlotOut TEMPERATURE;
    GNUPlotOut MOVEMENT;

    //** Всего у нас 3 дифура **//
    RungeType MM(StartX, DestinationX, DeltaX, Epsilon, M, ThermalMass);
    RungeType MS(StartX, DestinationX, DeltaX, Epsilon, W, ThermalSpeed);
    RungeType MT(StartX, DestinationX, DeltaX, Epsilon, T, ThermalTemperature);
    RungeType ML(StartX, DestinationX, DeltaX, Epsilon, Z, ThermalLift);

    int endStatus;
    bool Recalc[4] = {false};
    double RecommendDeltaStep[4];
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
        endStatus |= ML.MakeStep();

        RecommendDeltaStep[0] = MM.CalculatingDeltaStep(Recalc[0]);
        RecommendDeltaStep[1] = MS.CalculatingDeltaStep(Recalc[1]);
        RecommendDeltaStep[2] = MT.CalculatingDeltaStep(Recalc[2]);
        RecommendDeltaStep[3] = ML.CalculatingDeltaStep(Recalc[3]);

        NewDeltaStep = min(RecommendDeltaStep[0], min(min(RecommendDeltaStep[1], RecommendDeltaStep[2]), RecommendDeltaStep[3]));

        MM.SetDeltaStep(NewDeltaStep);
        MS.SetDeltaStep(NewDeltaStep);
        MT.SetDeltaStep(NewDeltaStep);
        ML.SetDeltaStep(NewDeltaStep);

        //*** Необходимо, чтобы выполнялась заданная точность ***//
        if(Recalc[0] || Recalc[1] || Recalc[2] || Recalc[3]) {
                cout << "FUS!" << endl;
            continue;
        }

        MM.ConfirmStep();
        MS.ConfirmStep();
        MT.ConfirmStep();
        ML.ConfirmStep();

        M = MM.GetCurrentValue();
        W = MS.GetCurrentValue();
        //if(W < 0) break;
        T = MT.GetCurrentValue();
        Z = ML.GetCurrentValue();
        if(Z < 0) break;

        //BurningRate = BurningRate/(120 - MM.GetCurrentStep());
        //if(BurningRate < 0) BurningRate = 0;


        if(MM.GetCurrentStep() - lastWrite > writeStep) {
            lastWrite = MM.GetCurrentStep();

            MASS.addCoordinate(M, MM.GetCurrentStep());
            SPEED.addCoordinate(W, MS.GetCurrentStep());
            TEMPERATURE.addCoordinate(T, MT.GetCurrentStep());
            MOVEMENT.addCoordinate(Z, ML.GetCurrentStep());

            cout << MM.GetCurrentStep() << endl;
        }

        StepCounter++;
    }

    cout << endl;
    cout << "End Thermal Mass: " << M << endl;
    cout << "End Thermal Speed: " << W << endl;
    cout << "End Thermal Temperature: " << T << endl;
    cout << "End Thermal Lift: " << Z << endl;

    cout << endl << "StepCounter: " << StepCounter << "\t" << "AllStepCounter: " << AllStepCounter << endl;

    MASS.out("MASS.txt");
    SPEED.out("SPEED.txt");
    TEMPERATURE.out("TEMPERATURE.txt");
    MOVEMENT.out("MOVEMENT.txt");

    cin.get();

    return 0;
}
