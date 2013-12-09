#include <iostream>
#include <iomanip>
#include <math.h>

#include <ostream>
#include <fstream>

using namespace std;


const int CalculationNotRequired = 1;

typedef double(*DifferentialFunction)(double CurrentStep, double CurrentValue, double DeltaStep);

class RungeKuttaMethod {
public:
    virtual int MakeStep() = 0;
    virtual void GetSolution(ostream &stream) = 0;

    RungeKuttaMethod(const double StartStep, const double EndStep, double DeltaStep,
                     const double Epsilon,
                     DifferentialFunction TargetFunction);

protected:
    virtual double CalculatingDeltaStep(bool &Recalc) = 0;

    const double StartStep;
    const double EndStep;
    double DeltaStep;
    double CurrentStep;
    const double Epsilon; //Допустимая погрешность

    DifferentialFunction TargetFunction;
};

RungeKuttaMethod::RungeKuttaMethod(const double StartStep, const double EndStep, double DeltaStep,
                                   const double Epsilon,
                                   DifferentialFunction TargetFunction)
        : StartStep(StartStep), EndStep(EndStep), DeltaStep(DeltaStep),
          CurrentStep(StartStep + DeltaStep), Epsilon(Epsilon),
          TargetFunction(TargetFunction) { }



//Метод Рунге - Кутты первого порядка точности
class EulerMethod : public RungeKuttaMethod
{
public:
    int MakeStep();
    void GetSolution(ostream &stream);

    EulerMethod(const double StartStep, const double EndStep, double DeltaStep,
                     const double Epsilon,
                     DifferentialFunction TargetFunction);

    void PrintStatistic();

private:
    double CalculatingDeltaStep(bool &Recalc);

    double CurrentValue, NextValue;

    int StepCount;
    int Low, High;
};

EulerMethod::EulerMethod(const double StartStep, const double EndStep, double DeltaStep,
                        const double Epsilon,
                        DifferentialFunction TargetFunction)
    : RungeKuttaMethod(StartStep, EndStep, DeltaStep, Epsilon, TargetFunction)
{
    CurrentValue = 1;
    StepCount = 0;
    Low = 0;
    High = 0;
}

int EulerMethod::MakeStep() {
    int StepStatus = 0;

    if((CurrentStep + DeltaStep) >= EndStep) {
            DeltaStep = EndStep - CurrentStep;
            StepStatus = CalculationNotRequired;
    }

    NextValue = CurrentValue + DeltaStep * TargetFunction(CurrentStep, CurrentValue, DeltaStep);
    StepCount++;

    return StepStatus;
}

void EulerMethod::GetSolution(ostream &stream) {
    double NewDeltaStep;
    bool Recalc = false;

    stream << CurrentValue << " " << CurrentStep - DeltaStep << endl;

    while(true) {
        do {
            if(MakeStep() == CalculationNotRequired) {
                stream << NextValue << " " << CurrentStep + DeltaStep << endl;
                return;
            }

            NewDeltaStep = CalculatingDeltaStep(Recalc);
            DeltaStep = NewDeltaStep;
        } while(Recalc);

        CurrentStep+=DeltaStep;
        stream << NextValue << " " << CurrentStep << " " << DeltaStep << endl;

        CurrentValue = NextValue;
    }
}

double EulerMethod::CalculatingDeltaStep(bool &Recalc) {
    double Sigma = 0.5*DeltaStep*fabs(NextValue - CurrentValue); //Погрешность метода
    double CorrectionStep = sqrt(Epsilon/Sigma);

    if(CorrectionStep < 1) {
        Low++;
        Recalc = true;
    } else {
        High++;
        Recalc = false;
    }

    return CorrectionStep*DeltaStep/1.1;
}


void EulerMethod::PrintStatistic() {
    cout << "StepCount: " << StepCount << endl;
    cout << "\t" << "LOW: " << Low << "" << endl;
    cout << "\t" << "HIGH: " << High << "" << endl;
    cout << "\t" << "Result: " << CurrentValue << endl;
    cout << "\t" << "LastDelta: " << DeltaStep << endl << endl;
}





class MersonMethod : public RungeKuttaMethod
{
public:
    int MakeStep();
    void GetSolution(ostream &stream);

    double GetCurrentValue() {
        return CurrentValue;
    }

    double GetCurrentStep() {
        return CurrentStep;
    }

    void ConfirmStep() {
        CurrentStep += DeltaStep;
        CurrentValue = NextValue;
    }

    void SetDeltaStep(double DeltaStep) {
        this->DeltaStep = DeltaStep;
    }

    MersonMethod(const double StartStep, const double EndStep, double DeltaStep,
           const double Epsilon,
           const double StartCondition,
           DifferentialFunction TargetFunction);

    double CalculatingDeltaStep(bool &Recalc);

    void PrintStatistic();

private:
    double k1, k2, k3, k4, k5; //Этапы вычисления по Мерсону

    double CurrentValue, NextValue;

    int StepCount;
    int Low, High;
};

MersonMethod::MersonMethod(const double StartStep, const double EndStep, double DeltaStep,
               const double Epsilon,
               const double StartCondition,
               DifferentialFunction TargetFunction)
    : RungeKuttaMethod(StartStep, EndStep, DeltaStep, Epsilon, TargetFunction),
      CurrentValue(StartCondition)
{
    StepCount = 0;
    Low = 0;
    High = 0;
}

int MersonMethod::MakeStep() {
    int StepStatus = 0;

    if((CurrentStep + DeltaStep) >= EndStep) {
            DeltaStep = EndStep - CurrentStep;
            StepStatus = CalculationNotRequired;
    }

    k1 = 1./3 * DeltaStep * TargetFunction(CurrentStep,                    CurrentValue,                            DeltaStep);
    k2 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./3 * DeltaStep, CurrentValue + k1,                       DeltaStep);
    k3 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./3 * DeltaStep, CurrentValue + 1./2*k1 + 1./2*k2,        DeltaStep);
    k4 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./2 * DeltaStep, CurrentValue + 3./8*k1 + 9./8*k3,        DeltaStep);
    k5 = 1./3 * DeltaStep * TargetFunction(CurrentStep +        DeltaStep, CurrentValue + 3./2*k1 - 9./2*k3 + 6*k4, DeltaStep);

    NextValue = CurrentValue + 1./2*(k1 + 4*k4 + k5);

    StepCount++;

    return StepStatus;
}

void MersonMethod::GetSolution(ostream &stream) {
    double NewDeltaStep;
    bool Recalc = false;

    stream << CurrentValue << " " << CurrentStep - DeltaStep << endl;

    while(true) {
        do {
            if(MakeStep() == CalculationNotRequired) {
                stream << NextValue << " " << CurrentStep + DeltaStep << endl;
                return;
            }

            NewDeltaStep = CalculatingDeltaStep(Recalc);
            DeltaStep = NewDeltaStep;
        } while(Recalc);

        CurrentStep+=DeltaStep;
        stream << NextValue << " " << CurrentStep << " " << DeltaStep << endl;

        CurrentValue = NextValue;
    }
}

double MersonMethod::CalculatingDeltaStep(bool &Recalc) {
    //double Sigma = k1 - (9./2)*k3 + 4*k4 - (1./2)*k5; //Погрешность метода
    //double Sigma = 1.0/30.0*(2*k1-9*k3+8.0*k4-k5);
    double Sigma = 0.2*DeltaStep*(NextValue - CurrentValue);
    Sigma = fabs(Sigma);
    //cout << Sigma << endl;

    if(Sigma > Epsilon*5) {
        Low++;
        Recalc = true;
        return DeltaStep/2;
    } else if(Sigma < 5./32*Epsilon){
        Recalc = false;
        High++;
        return DeltaStep*2;
    } else {
        Recalc = false;
        return DeltaStep;
    }

    //return DeltaStep*=pow(Epsilon/Sigma, 1.0/6.0);
}


void MersonMethod::PrintStatistic() {
    cout << "StepCount: " << StepCount << endl;
    cout << "\t" << "LOW: " << Low << "" << endl;
    cout << "\t" << "HIGH: " << High << "" << endl;
    cout << "\t" << "Result: " << CurrentValue << endl;
    cout << "\t" << "LastDelta: " << DeltaStep << endl << endl;
}

const double pi = 3.14159265359;
const double g = 9.80665;

const double TreeCrownsHeight = 10;

const double PressureEnviroment = 101325.0;   //В паскалях
const double AirUniversalGasConstant = 8.314/(6.022*pow((double)(10.0),(double)(-4.0))*(46.5+3.3+53.2)/3.0);

const double alpha = 0.6; //Коэффициент вовлечения воздуха из приземного слоя атмосферы в термик
const double betta = 0.6; //Коэффициент вовлечения продуктов горения в термик (0..1);

double BurningRate = 1; //скорость горения лесных горючих материалов (ЛГМ) [кг/с];

double TemperatureEnviroment = 288;
double Wg = 1;
double Tg = 293;
double cd = 1;
double cp = 1;

double M; //Масса термика
double W; //Скорость термика
double T = 323; //Температура термика

double TestFunction(double CurrentStep, double CurrentValue, double DeltaStep) {
    return pow(CurrentStep, 3) - CurrentValue;
}

double ThermalMass(double CurrentStep, double CurrentValue, double DeltaStep) {
    return alpha*CurrentValue*W + betta*BurningRate;
}

double ThermalSpeed(double CurrentStep, double CurrentValue, double DeltaStep) {
    double Ro = PressureEnviroment/(AirUniversalGasConstant*T);
    return ((T - TemperatureEnviroment)*g)/TemperatureEnviroment +
           (betta*BurningRate*(Wg-CurrentValue))/M - alpha*CurrentValue*CurrentValue -
           ((pi*cd)/(2*M))*pow((3*M)/(4*pi), 2./3) * pow(Ro, 1./3) * CurrentValue*CurrentValue;
}

double ThermalTemperature(double CurrentStep, double CurrentValue, double DeltaStep) {
    return ((BurningRate*betta)/M)*(Tg + TemperatureEnviroment - 2*CurrentValue) +
           (TemperatureEnviroment - CurrentValue)*alpha*W - (g*W*CurrentValue)/(cp*TemperatureEnviroment);
}

int main() {
    const double DeltaX = 0.0001;  //Шаг по координате
    const double StartX = 0;
    const double DestinationX = 6;
    const double Epsilon = 0.0001;

    //Задание начальных условий
    M = (pi*PressureEnviroment*(TreeCrownsHeight/2))/(6*AirUniversalGasConstant*323)*pow(TreeCrownsHeight, 3);
    cout << "Start Thermal Mass: " << M << endl;

    W = sqrt((g*pow(TreeCrownsHeight, 2)*(323-288))/288);
    cout << "Start Thermal Speed: " << W << endl;

    T = 323;
    cout << "Start Thermal Temperature: " << T << endl;

    //ofstream Euler("Euler.txt");
    ofstream MASS("MASS.txt");
    ofstream SPEED("SPEED.txt");
    ofstream TEMPERATURE("TEMPERATURE.txt");
    ofstream ALL("ALL.txt");

    //EulerMethod EM(StartX, DestinationX, DeltaX, Epsilon, ThermalMass);
    //EM.GetSolution(Euler);
    //EM.PrintStatistic();

    MersonMethod MM(StartX, DestinationX, DeltaX, Epsilon, M, ThermalMass);
    MersonMethod MS(StartX, DestinationX, DeltaX, Epsilon, W, ThermalSpeed);
    MersonMethod MT(StartX, DestinationX, DeltaX, Epsilon, T, ThermalTemperature);

    int endStatus;
    bool Recalc[3] = {false};
    double RecommendDeltaStep[3];
    double NewDeltaStep;

    size_t StepCounter = 0;
    size_t AllStepCounter = 0;

    while(!endStatus) {
        endStatus = 0;
        AllStepCounter++;

        endStatus |= MM.MakeStep();
        endStatus |= MS.MakeStep();
        endStatus |= MT.MakeStep();

        RecommendDeltaStep[0] = MM.CalculatingDeltaStep(Recalc[0]);
        RecommendDeltaStep[1] = MS.CalculatingDeltaStep(Recalc[1]);
        RecommendDeltaStep[2] = MT.CalculatingDeltaStep(Recalc[2]);

        NewDeltaStep = min(RecommendDeltaStep[0], min(RecommendDeltaStep[1], RecommendDeltaStep[2]));

        MM.SetDeltaStep(NewDeltaStep);
        MS.SetDeltaStep(NewDeltaStep);
        MT.SetDeltaStep(NewDeltaStep);

        if(Recalc[0] || Recalc[1] || Recalc[2]) {
            continue;
        }



        MM.ConfirmStep();
        MS.ConfirmStep();
        MT.ConfirmStep();

        M = MM.GetCurrentValue();
        W = MS.GetCurrentValue();
        T = MT.GetCurrentValue();

        ALL << M << "\t" << W << "\t" << T << endl;

        MASS         << M << " " << MM.GetCurrentStep()*1000 << endl;
        SPEED        << W << " " << MS.GetCurrentStep()*1000 << endl;
        TEMPERATURE  << T << " " << MT.GetCurrentStep()*1000 << endl;
        StepCounter++;
    }

    cout << endl;
    cout << "End Thermal Mass: " << M << endl;
    cout << "End Thermal Speed: " << W << endl;
    cout << "End Thermal Temperature: " << T << endl;

    cout << endl << "StepCounter: " << StepCounter << "\t" << "AllStepCounter: " << AllStepCounter << endl;

    return 0;
}
