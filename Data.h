#ifndef DATA_H
#define DATA_H

#include <QFile>
#include <QVector>
#include <cmath>
#include <QGenericMatrix>
#include <QTextStream>
#include <iostream>
#include <iomanip>


using namespace std;



struct Solutions
{
    double _Hz;
    double _HzA;
    double _H2;
    double _H2A;

    bool operator <(const Solutions& sol) const;
    Solutions operator +(const Solutions &sol);
    Solutions operator -(const Solutions &sol);
    Solutions operator *(const Solutions &sol);
};




class Data
{
public:
    Data();
    ~Data();
    void parse(QFile &file, QFile &Gfile, QFile &Lfile, QFile &Sfile);

    void SetHz(Solutions& sol ,double Hz) { sol._Hz=Hz; }
    void SetHzA(Solutions& sol ,double HzA) { sol._HzA=HzA; }
    void SetH2(Solutions& sol ,double H2) { sol._H2=H2; }
    void SetH2A(Solutions& sol ,double H2A) { sol._H2A=H2A; }

    double Hz(Solutions sol) { return sol._Hz; }
    double HzA(Solutions sol) { return sol._HzA; }
    double H2(Solutions sol) { return sol._H2; }
    double H2A(Solutions sol) { return sol._H2A; }


    QVector<double> SolveQubicEquationCO(QGenericMatrix<3,3,double> M);
    QVector<double> SolveQubicEquationSF(QGenericMatrix<3,3,double> M);
    QVector<double> SolveQubicEquationSS(QGenericMatrix<3,3,double> M);
    void PrintSolutions(double n, double T, QFile &Lfile, QFile &Sfile);
    void PrintGlobalSolutions(double n, double T, QFile &Gfile);
    double Norm(double &x1, double &y1, double &z1, double &w1, double &x2, double &y2, double &z2, double &w2);
    void MiniSign(double &n, double &T, Solutions point1, Solutions point2);
    void Sign(double &n, double &T, Solutions Min, Solutions Max);
    Solutions Solve(double &n, double &T, double &D, double &HzMin, double &HzMax, double &HzA, double &H2, double &H2A);
    QGenericMatrix<4,4,double> FindInverse(QGenericMatrix<4,4,double>& M);
    bool IsEqual(double n, double T, Solutions sol);
    QGenericMatrix<3, 3, double> IsMinimumNO(double T, Solutions sol);
    QGenericMatrix<3, 3, double> IsMinimumCO(double T, Solutions sol);
    QGenericMatrix<3, 3, double> IsMinimumSF(double T, Solutions sol);
    QGenericMatrix<3, 3, double> IsMinimumSS(double T, Solutions sol);
    void FindSolutions(double n, double T);
    Solutions RealSolve(double &n, double &T, Solutions sol);


private:
    double _nMin;
    double _nMax;
    double _TMax;
    double _TMin;
    double _D;
    double _V;
    double _tb;
    QVector<double> *_Data;
    QVector<double> *_FEnergy;
    QVector<Solutions> *_Sign;
    QVector<Solutions> *_Solutions;

    double _ExtrType;
    double _Accuracy;
    double _GlobalFEnergy;

    //SupplyMentAry

    double SqP(double &T, Solutions sol)
    {
        return sqrt((Hz(sol)+HzA(sol))*(Hz(sol)+HzA(sol))+(H2(sol)+H2A(sol))*(H2(sol)+H2A(sol)))/T;
    }

    double SqM(double &T, Solutions sol)
    {
        return sqrt((Hz(sol)-HzA(sol))*(Hz(sol)-HzA(sol))+(H2(sol)-H2A(sol))*(H2(sol)-H2A(sol)))/T;
    }

    double ZcP(double &T, Solutions sol)
    {
        return 2*(1+exp(-_D/T)*cosh(SqP(T,sol)));
    }

    double ZcM(double &T, Solutions sol)
    {
        return 2*(1+exp(-_D/T)*cosh(SqM(T,sol)));
    }

    double PhiP(double &T, Solutions sol)
    {   
        double ok;

        double ret1=sinh(SqP(T,sol))/(SqP(T,sol)*(1+exp(-_D/T)*cosh(SqP(T,sol))));;

        double ret2=exp(_D/T)*tanh(SqP(T,sol))/SqP(T,sol);

        if (qIsNaN(ret1) || qIsInf(ret1)) ok=ret2;
        else ok=ret1;

        return ok;
    }

    double PhiM(double &T, Solutions sol)
    {
        double ok;

        double ret1=sinh(SqM(T,sol))/(SqM(T,sol)*(1+exp(-_D/T)*cosh(SqM(T,sol))));;

        double ret2=exp(_D/T)*tanh(SqM(T,sol))/SqM(T,sol);

                if (qIsNaN(ret1) || qIsInf(ret1)) ok=ret2;
                else ok=ret1;

        return ok;
    }

    double ThetaP(double &T, Solutions sol)
    {
        double ok;

        double ret1=(cosh(SqP(T,sol))+exp(-_D/T)-sinh(SqP(T,sol))*(1+exp(-_D/T)*cosh(SqP(T,sol)))/SqP(T,sol))
                /
                (SqP(T,sol)*SqP(T,sol)*(1+exp(-_D/T)*cosh(SqP(T,sol)))*(1+exp(-_D/T)*cosh(SqP(T,sol))));

        double ret2=exp(2*_D/T)/(SqP(T,sol)*SqP(T,sol)*cosh(SqP(T,sol)))-exp(_D/T)*tanh(SqP(T,sol))/(SqP(T,sol)*SqP(T,sol)*SqP(T,sol));

                if (qIsNaN(ret1) || qIsInf(ret1)) ok=ret2;
                else ok=ret1;

        return ok;
    }

    double ThetaM(double &T, Solutions sol)
    {
        double ok;

        double ret1=(cosh(SqM(T,sol))+exp(-_D/T)-sinh(SqM(T,sol))*(1+exp(-_D/T)*cosh(SqM(T,sol)))/SqM(T,sol))
                /
                (SqM(T,sol)*SqM(T,sol)*(1+exp(-_D/T)*cosh(SqM(T,sol)))*(1+exp(-_D/T)*cosh(SqM(T,sol))));

        double ret2=exp(2*_D/T)/(SqM(T,sol)*SqM(T,sol)*cosh(SqM(T,sol)))-exp(_D/T)*tanh(SqM(T,sol))/(SqM(T,sol)*SqM(T,sol)*SqM(T,sol));

                if (qIsNaN(ret1) || qIsInf(ret1)) ok=ret2;
                else ok=ret1;

        return ok;
    }

    double PiP(double &T, Solutions sol)
    {
        return (Hz(sol)+HzA(sol))*(Hz(sol)+HzA(sol))*ThetaP(T,sol)/(T*T)+PhiP(T,sol);
    }

    double PiM(double &T, Solutions sol)
    {
        return (Hz(sol)-HzA(sol))*(Hz(sol)-HzA(sol))*ThetaM(T,sol)/(T*T)+PhiM(T,sol);
    }

    double PsiP(double &T, Solutions sol)
    {
        return (Hz(sol)+HzA(sol))*(H2(sol)+H2A(sol))*ThetaP(T,sol)/(T*T);
    }

    double PsiM(double &T, Solutions sol)
    {
        return (Hz(sol)-HzA(sol))*(H2(sol)-H2A(sol))*ThetaM(T,sol)/(T*T);
    }

    double DeltaP(double &T, Solutions sol)
    {
        return (H2(sol)+H2A(sol))*(H2(sol)+H2A(sol))*ThetaP(T,sol)/(T*T)+PhiP(T,sol);
    }

    double DeltaM(double &T, Solutions sol)
    {
        return (H2(sol)-H2A(sol))*(H2(sol)-H2A(sol))*ThetaM(T,sol)/(T*T)+PhiM(T,sol);
    }

    //Second Derivatives

    double aHz(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(PiP(T,sol)-PiM(T,sol))/T;
    }

    double aHzA(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(PiP(T,sol)+PiM(T,sol))/T;
    }

    double aH2(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(PsiP(T,sol)-PsiM(T,sol))/T;
    }

    double aH2A(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(PsiP(T,sol)+PsiM(T,sol))/T;
    }

    double BH2(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(DeltaP(T,sol)+DeltaM(T,sol))/T;
    }

    double BH2A(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(DeltaP(T,sol)-DeltaM(T,sol))/T;
    }

    double f2Hz(double &T, Solutions sol)
    {
        return -4*_V*aHz(T,sol)*aHz(T,sol)-2*_tb*(aH2A(T,sol)*aH2A(T,sol)-aH2(T,sol)*aH2(T,sol))-aHzA(T,sol);
    }

    double f2HzHzA(double &T, Solutions sol)
    {
        return -4*_V*aHz(T,sol)*aHzA(T,sol);
    }

    double f2HzH2(double &T, Solutions sol)
    {
        return -4*_V*aHz(T,sol)*aH2(T,sol)-2*_tb*(aH2A(T,sol)*BH2(T,sol)-aH2(T,sol)*BH2A(T,sol));
    }

    double f2HzH2A(double &T, Solutions sol)
    {
        return -4*_V*aHz(T,sol)*aH2A(T,sol)-2*_tb*(aH2A(T,sol)*BH2A(T,sol)-aH2(T,sol)*BH2(T,sol));
    }

    double f2HzA(double &T, Solutions sol)
    {
        return -4*_V*aHzA(T,sol)*aHzA(T,sol)-2*_tb*(aH2(T,sol)*aH2(T,sol)-aH2A(T,sol)*aH2A(T,sol))+aHzA(T,sol);
    }

    double f2HzAH2(double &T, Solutions sol)
    {
        return -4*_V*aHzA(T,sol)*aH2(T,sol)-2*_tb*(aH2(T,sol)*BH2(T,sol)-aH2A(T,sol)*BH2A(T,sol))+aH2(T,sol);
    }

    double f2HzAH2A(double &T, Solutions sol)
    {
        return -4*_V*aHzA(T,sol)*aH2A(T,sol)-2*_tb*(aH2(T,sol)*BH2A(T,sol)-aH2A(T,sol)*BH2(T,sol))+aH2A(T,sol);
    }

    double f2H2(double &T, Solutions sol)
    {
        return -4*_V*aH2(T,sol)*aH2(T,sol)-2*_tb*(BH2(T,sol)*BH2(T,sol)-BH2A(T,sol)*BH2A(T,sol))+BH2(T,sol);
    }

    double f2H2A(double &T, Solutions sol)
    {
        return -4*_V*aH2A(T,sol)*aH2A(T,sol)-2*_tb*(BH2A(T,sol)*BH2A(T,sol)-BH2(T,sol)*BH2(T,sol))+BH2(T,sol);
    }

    double f2H2H2A(double &T, Solutions sol)
    {
        return -4*_V*aH2(T,sol)*aH2A(T,sol)+BH2A(T,sol);
    }

    QGenericMatrix<3,3,double> FEnergySecondDerivative(double &T, Solutions sol)
    {
        QGenericMatrix<3,3,double> FEn2;

        FEn2.fill(0);

        FEn2(0,0)=f2HzA(T,sol)+aHz(T,sol)*aHz(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-2*aHz(T,sol)*f2HzHzA(T,sol)/aHzA(T,sol);
        FEn2(0,1)=f2HzAH2(T,sol)+aHz(T,sol)*aH2A(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aHz(T,sol)*f2HzH2(T,sol)/aHzA(T,sol)-aH2A(T,sol)*f2HzHzA(T,sol)/aHzA(T,sol);
        FEn2(0,2)=f2HzAH2A(T,sol)+aHz(T,sol)*aH2(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aHz(T,sol)*f2HzH2A(T,sol)/aHzA(T,sol)-aH2(T,sol)*f2HzHzA(T,sol)/aHzA(T,sol);

        FEn2(1,0)=f2HzAH2(T,sol)+aHz(T,sol)*aH2A(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aHz(T,sol)*f2HzH2(T,sol)/aHzA(T,sol)-aH2A(T,sol)*f2HzHzA(T,sol)/aHzA(T,sol);
        FEn2(1,1)=f2H2(T,sol)+aH2A(T,sol)*aH2A(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-2*aH2A(T,sol)*f2HzH2(T,sol)/aHzA(T,sol);
        FEn2(1,2)=f2H2H2A(T,sol)+aH2(T,sol)*aH2A(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aH2A(T,sol)*f2HzH2A(T,sol)/aHzA(T,sol)-aH2(T,sol)*f2HzH2(T,sol)/aHzA(T,sol);

        FEn2(2,0)=f2HzAH2A(T,sol)+aHz(T,sol)*aH2(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aHz(T,sol)*f2HzH2A(T,sol)/aHzA(T,sol)-aH2(T,sol)*f2HzHzA(T,sol)/aHzA(T,sol);
        FEn2(2,1)=f2H2H2A(T,sol)+aH2(T,sol)*aH2A(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-aH2A(T,sol)*f2HzH2A(T,sol)/aHzA(T,sol)-aH2(T,sol)*f2HzH2(T,sol)/aHzA(T,sol);
        FEn2(2,2)=f2H2A(T,sol)+aH2(T,sol)*aH2(T,sol)*f2Hz(T,sol)/(aHzA(T,sol)*aHzA(T,sol))-2*aH2(T,sol)*f2HzH2A(T,sol)/aHzA(T,sol);

        return FEn2;
    }


    //Equations

    double a(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*((Hz(sol)+HzA(sol))*PhiP(T,sol)-(Hz(sol)-HzA(sol))*PhiM(T,sol))/T;
    }

    double B(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*((H2(sol)+H2A(sol))*PhiP(T,sol)+(H2(sol)-H2A(sol))*PhiM(T,sol))/T;
    }

    double Ba(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*((H2(sol)+H2A(sol))*PhiP(T,sol)-(H2(sol)-H2A(sol))*PhiM(T,sol))/T;
    }

    double EqN(double &n, double &T, Solutions sol)
    {
        return n-0.5*exp(-_D/T)*((Hz(sol)+HzA(sol))*PhiP(T,sol)+(Hz(sol)-HzA(sol))*PhiM(T,sol))/T;
    }

    double EqCO(double &n,double &T, Solutions sol)
    {
        return  HzA(sol)-4*_V*a(T,sol);
    }

    double EqSF(double &n,double &T, Solutions sol)
    {
        return  H2(sol)-2*_tb*B(T,sol);
    }

    double EqASF(double &n,double &T, Solutions sol)
    {
        return  H2A(sol)+2*_tb*Ba(T,sol);
    }

    //FEnergy

    double FEnergy(double n,double &T, Solutions sol)
    {
        double okP;
        double okM;

        double ret1P=-0.5*T*log(ZcP(T,sol));
        double ret1M=-0.5*T*log(ZcM(T,sol));

        double ret2P=-0.5*T*(SqP(T,sol)-_D/T);
        double ret2M=-0.5*T*(SqM(T,sol)-_D/T);

                if (qIsNaN(ret1P) || qIsInf(ret1P)) okP=ret2P;
                else okP=ret1P;

                if (qIsNaN(ret1M) || qIsInf(ret1M)) okM=ret2M;
                else okM=ret1M;

        return  okP
                +
                okM
                +
                2*_V*(n*n-a(T,sol)*a(T,sol))
                -
                _tb*(B(T,sol)*B(T,sol)-Ba(T,sol)*Ba(T,sol))
                +
                H2(sol)*B(T,sol)
                +
                H2A(sol)*Ba(T,sol)
                +
                HzA(sol)*a(T,sol)
                +
                Hz(sol)*n;
    }

};

#endif // DATA_H
