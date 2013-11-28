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
        stream << NextValue << " " << CurrentStep << endl;

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

const double DeltaX = 0.0001;  //Шаг по координате
const double StartX = 0;
const double DestinationX = 2;
const double Epsilon = 0.0001;

double TestFunction(double CurrentStep, double CurrentValue, double DeltaStep) {
    return pow(CurrentStep, 3) - CurrentValue;
}

int main() {
    ofstream Test("Test.txt");

    EulerMethod EM(StartX, DestinationX, DeltaX, Epsilon, TestFunction);
    EM.GetSolution(Test);
    EM.PrintStatistic();

    return 0;
}
