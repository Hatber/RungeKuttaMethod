#include <ostream> //Потоковая запись

namespace RK {
    using std::ostream;
    using std::endl;

    //*** Возвращается методом "int MakeStep()" когда вычисления законченны
    const int CalculationNotRequired = 1;

    //*** Прототип для дифференциальных уравнений
    typedef double(*DifferentialFunction)(double CurrentStep, double CurrentValue);

    /********************************************
     *  RungeKuttaMethod  -  виртуальный класс. *
     *  Он описывает интерфейс к  Явным методам *
     *  Рунге  -  Кутты  произвольного  порядка *
     *  точности.                               *
     *                                          *
     *  Для этого необходимы реализация методов *
     *  MakeStep() и  CalculatingDeltaStep()  в *
     *  производных классах.                    *
     *                                          *
     *  Также  он  реализует  шаблонный   метод *
     *  GetSolution(),   позволяющий   получить *
     *  решение  одного  ДифУр'a  полностью,  с *
     *  записью результатов в поток.            *
     ********************************************/
    class RungeKuttaMethod {
    public:
        virtual int MakeStep() = 0;
        virtual double CalculatingDeltaStep(bool &Recalc) = 0;

        RungeKuttaMethod(const double StartStep, const double EndStep, double DeltaStep,
                         const double Epsilon, const double StartCondition,
                         DifferentialFunction TargetFunction):
            StartStep(StartStep), EndStep(EndStep), DeltaStep(DeltaStep),
            CurrentStep(StartStep + DeltaStep), Epsilon(Epsilon), CurrentValue(StartCondition),
            TargetFunction(TargetFunction) { };

        void GetSolution(ostream &stream) {
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

        //** Текущее значение
        double GetCurrentValue() {
            return CurrentValue;
        }

        //** Пройденное растояние по временной координате
        double GetCurrentStep() {
            return CurrentStep;
        }

        //** Устанавливает новое значение шага по времени
        void SetDeltaStep(double DeltaStep) {
            this->DeltaStep = DeltaStep;
        }

        //** Подтвердить выполнение шага
        void ConfirmStep() {
            CurrentStep += DeltaStep;
            CurrentValue = NextValue;
        }

    protected:
        const double StartStep;
        const double EndStep;
        double DeltaStep;
        double CurrentStep;
        const double Epsilon; //Допустимая погрешность

        double CurrentValue, NextValue;

        DifferentialFunction TargetFunction;
    };

    //Метод Рунге - Кутты первого порядка точности
    class EulerMethod : public RungeKuttaMethod
    {
    public:
        EulerMethod(const double StartStep, const double EndStep, double DeltaStep,
                         const double Epsilon, const double StartCondition,
                         DifferentialFunction TargetFunction) :
            RungeKuttaMethod(StartStep, EndStep, DeltaStep, Epsilon, StartCondition, TargetFunction) { }

        int MakeStep() {
            int StepStatus = 0;

            if((CurrentStep + DeltaStep) >= EndStep) {
                    DeltaStep = EndStep - CurrentStep;
                    StepStatus = CalculationNotRequired;
            }

            NextValue = CurrentValue + DeltaStep * TargetFunction(CurrentStep, CurrentValue);

            return StepStatus;
        }

        double CalculatingDeltaStep(bool &Recalc) {
            double Sigma = 0.5*DeltaStep*fabs(NextValue - CurrentValue); //Погрешность метода
            double CorrectionStep = sqrt(Epsilon/Sigma);

            if(CorrectionStep < 1) {
                Recalc = true;
            } else {
                Recalc = false;
            }

            return CorrectionStep*DeltaStep/1.1;
        }
    };

    //Метод Рунге - Кутты пятого порядка точности
    class MersonMethod : public RungeKuttaMethod
    {
    public:
        MersonMethod(const double StartStep, const double EndStep, double DeltaStep,
               const double Epsilon,
               const double StartCondition,
               DifferentialFunction TargetFunction) :
            RungeKuttaMethod(StartStep, EndStep, DeltaStep, Epsilon, StartCondition, TargetFunction) { }

        int MakeStep() {
            int StepStatus = 0;

            if((CurrentStep + DeltaStep) >= EndStep) {
                    DeltaStep = EndStep - CurrentStep;
                    StepStatus = CalculationNotRequired;
            }

            k1 = 1./3 * DeltaStep * TargetFunction(CurrentStep,                    CurrentValue);
            k2 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./3 * DeltaStep, CurrentValue + k1);
            k3 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./3 * DeltaStep, CurrentValue + 1./2*k1 + 1./2*k2);
            k4 = 1./3 * DeltaStep * TargetFunction(CurrentStep + 1./2 * DeltaStep, CurrentValue + 3./8*k1 + 9./8*k3);
            k5 = 1./3 * DeltaStep * TargetFunction(CurrentStep +        DeltaStep, CurrentValue + 3./2*k1 - 9./2*k3 + 6*k4);

            NextValue = CurrentValue + 1./2*(k1 + 4*k4 + k5);

            return StepStatus;
        }

        double CalculatingDeltaStep(bool &Recalc) {
            //double Sigma = k1 - (9./2)*k3 + 4*k4 - (1./2)*k5; //Погрешность метода
            //double Sigma = 1.0/30.0*(2*k1-9*k3+8.0*k4-k5);
            double Sigma = 0.2*DeltaStep*(NextValue - CurrentValue);
            Sigma = fabs(Sigma);

            if(Sigma > Epsilon*5) {
                Recalc = true;
                return DeltaStep/2;
            } else if(Sigma < 5./32*Epsilon){
                Recalc = false;
                return DeltaStep*2;
            } else {
                Recalc = false;
                return DeltaStep;
            }

            //return DeltaStep*=pow(Epsilon/Sigma, 1.0/6.0);
        }

    private:
        double k1, k2, k3, k4, k5; //Этапы вычисления по Мерсону
    };
} //namespace RK
