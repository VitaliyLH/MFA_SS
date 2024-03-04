#include "Data.h"

Data::Data()
{
    _Data = new QVector<double>;
    _FEnergy = new QVector<double>;
    _Sign = new QVector<Solutions>;
    _Solutions = new QVector<Solutions>;
}

Data::~Data()
{
    delete _Data;
    delete _Sign;
    delete _FEnergy;
    delete _Solutions;
}

bool Solutions::operator <(const Solutions& sol) const
{
    bool lessthan=true;

    if (fabs(_HzA-sol._HzA)>0.0001 && fabs(_H2-sol._H2)<0.0001 && fabs(_H2A-sol._H2A)<0.0001)
        lessthan=(_HzA <= sol._HzA);

    if (fabs(_HzA-sol._HzA)<0.0001 && fabs(_H2-sol._H2)>0.0001 && fabs(_H2A-sol._H2A)<0.0001)
        lessthan=(_H2 <= sol._H2);

    if (fabs(_HzA-sol._HzA)>0.0001 && fabs(_H2-sol._H2)>0.0001 && fabs(_H2A-sol._H2A)>0.0001)
        lessthan=(_HzA <= sol._HzA && _H2 <= sol._H2 && _H2A <= sol._H2A);

    return lessthan;
}


Solutions Solutions::operator +(const Solutions& sol)
{
    Solutions solp;
    solp._Hz=_Hz+sol._Hz;
    solp._HzA=_HzA+sol._HzA;
    solp._H2=_H2+sol._H2;
    solp._H2A=_H2A+sol._H2A;

    return solp;
}

Solutions Solutions::operator -(const Solutions &sol)
{
    Solutions solm;
    solm._Hz=_Hz-sol._Hz;
    solm._HzA=_HzA-sol._HzA;
    solm._H2=_H2-sol._H2;
    solm._H2A=_H2A-sol._H2A;

    return solm;
}

Solutions Solutions::operator *(const Solutions& sol)
{
    Solutions solp;
    solp._Hz=_Hz*sol._Hz;
    solp._HzA=_HzA*sol._HzA;
    solp._H2=_H2*sol._H2;
    solp._H2A=_H2A*sol._H2A;

    return solp;
}



void Data::parse(QFile &file, QFile &Gfile, QFile &Lfile, QFile &Sfile)
{
    QTextStream stream(&file);
    while (!stream.atEnd()) {

        double data;
        stream>>data;
        _Data->append(data);
    }



    _nMin=_Data->at(0);
    _nMax=_Data->at(1);
    _TMin=_Data->at(2);
    _TMax=_Data->at(3);
    _D=_Data->at(4);
    _V=_Data->at(5);
    _tb=_Data->at(6);

    _Accuracy=0.0000000001; //0.0000000000001;

    double Stepn=0.001;
    double StepT=0.001;

    Solutions Min;
    Solutions Max;

    for (double n=_nMin; n <= _nMax+0.01*Stepn; n+=Stepn)
    {

        for (double T=_TMin; T <= _TMax+0.01*StepT; T+=StepT)
        {
            SetHz(Min,0);
            SetHzA(Min,-0.001);
            SetH2(Min,-0.001);
            SetH2A(Min,-0.001);

            SetHz(Max,10.0*T+2*n);
            SetHzA(Max,0.001);
            SetH2(Max,0.001);
            SetH2A(Max,0.001);

            Sign(n,T, Min, Max);

            SetHz(Min,0);
            SetHzA(Min,0);
            SetH2(Min,0);
            SetH2A(Min,0);

            SetHz(Max,5);
            SetHzA(Max,4*_V+0.1);
            SetH2(Max,2*_tb+0.1);
            SetH2A(Max,2*_tb+0.1);

            Sign(n,T, Min, Max);

            FindSolutions(n,T);

            PrintGlobalSolutions(n, T, Gfile);

            //PrintSolutions(n, T, Lfile, Sfile);

            _Sign->clear();
            _Solutions->clear();
        }

    }

}







QVector<double> Data::SolveQubicEquationSS(QGenericMatrix<3,3,double> M)
{

    QVector<double> Lambda;

    double A11=M(0,0);
    double A12=M(0,1);
    double A13=M(0,2);

    double A22=M(1,1);
    double A23=M(1,2);
    double A33=M(2,2);

    double d11=A11;
    double d21=A22;
    double d31=A33;

    double d12=A22*A33-A23*A23;
    double d22=A11*A33-A13*A13;
    double d32=A11*A22-A12*A12;

    double d=A11*(A22*A33-A23*A23)-A12*(A12*A33-A13*A23)+A13*(A12*A23-A13*A22);

    if (fabs(d) < 0.001)
    {
        double b=(d11+d21+d31);
        double c=(d12+d22+d32);

        double D=b*b-4*c;

        double L1=0.5*(b+sqrt(D));
        double L2=0.5*(b-sqrt(D));
        double L3=0;

        double Eq1=L1*L1-b*L1+c;
        double Eq2=L2*L2-b*L2+c;

        if (D>=0 && fabs(Eq1) < 0.0000001 && fabs(Eq2) < 0.0000001)
        {
            Lambda.append(L1);
            Lambda.append(L2);
            Lambda.append(L3);
        }
    }

    return Lambda;
}

QVector<double> Data::SolveQubicEquationSF(QGenericMatrix<3,3,double> M)
{

    QVector<double> Lambda;

    double A11=M(0,0);
    double A13=M(0,2);
    double A22=M(1,1);
    double A33=M(2,2);

    double D=(A11+A33)*(A11+A33)-4*(A11*A33-A13*A13);

    double L1=0.5*((A11+A33)+sqrt(D));
    double L2=0.5*((A11+A33)-sqrt(D));
    double L3=A22;

    double Eq1=L1*L1-(A11+A33)*L1+(A11*A33-A13*A13);
    double Eq2=L2*L2-(A11+A33)*L2+(A11*A33-A13*A13);

    if (D>=0 && fabs(Eq1) < 0.0000001 && fabs(Eq2) < 0.0000001)
    {
        Lambda.append(L1);
        Lambda.append(L2);
        Lambda.append(L3);
    }

    return Lambda;
}


QVector<double> Data::SolveQubicEquationCO(QGenericMatrix<3,3,double> M)
{

    QVector<double> Lambda;

    double A11=M(0,0);

    double A22=M(1,1);
    double A23=M(1,2);
    double A33=M(2,2);

    double D=(A22+A33)*(A22+A33)-4*(A22*A33-A23*A23);

    double L1=0.5*((A22+A33)+sqrt(D));
    double L2=0.5*((A22+A33)-sqrt(D));
    double L3=A11;

    double Eq1=L1*L1-(A22+A33)*L1+(A22*A33-A23*A23);
    double Eq2=L2*L2-(A22+A33)*L2+(A22*A33-A23*A23);

    if (D>=0 && fabs(Eq1) < 0.0000001 && fabs(Eq2) < 0.0000001)
    {
        Lambda.append(L1);
        Lambda.append(L2);
        Lambda.append(L3);
    }

    return Lambda;
}



QGenericMatrix<3,3,double> Data::IsMinimumNO(double T, Solutions sol)
{
    _ExtrType=0;
    QGenericMatrix<3,3,double> DM=FEnergySecondDerivative(T,sol);

    if (DM(0,0)>=-_Accuracy && DM(1,1)>=-_Accuracy && DM(2,2)>=-_Accuracy) _ExtrType=0;
    else _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=3;

    //    if (DM(0,0)>=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=1;
    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;

    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=2;

    return DM;
}



QGenericMatrix<3,3,double> Data::IsMinimumCO(double T, Solutions sol)
{
    QVector<double> Lambda;
    QGenericMatrix<3,3,double> U;
    QGenericMatrix<3,3,double> M=FEnergySecondDerivative(T,sol);
    QGenericMatrix<3,3,double> DM;

    U.fill(0);

    double A11=M(0,0);
    double A22=M(1,1);
    double A23=M(1,2);
    double A33=M(2,2);

    if (fabs(A23) > 0.001)
    {

        Lambda=SolveQubicEquationCO(M);

        DM(0,0)=Lambda.at(0);
        DM(1,1)=Lambda.at(1);
        DM(2,2)=Lambda.at(2);

    }
    else
    {
        DM=M;
    }


    if (DM(0,0)>=-_Accuracy && DM(1,1)>=-_Accuracy && DM(2,2)>=-_Accuracy) _ExtrType=0;
    else _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=3;

    //    if (DM(0,0)>=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=1;
    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;

    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=2;



    for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++)
        {
            if (i!=j && fabs(DM(i,j))>0.001) _ExtrType=-1;
        }

    return DM;

}


QGenericMatrix<3,3,double> Data::IsMinimumSF(double T, Solutions sol)
{
    QVector<double> Lambda;
    QGenericMatrix<3,3,double> U;
    QGenericMatrix<3,3,double> M=FEnergySecondDerivative(T,sol);
    QGenericMatrix<3,3,double> DM;

    U.fill(0);

    double A11=M(0,0);
    double A13=M(0,2);
    double A22=M(1,1);
    double A33=M(2,2);

    if (fabs(A13) > 0.001)
    {

        Lambda=SolveQubicEquationSF(M);

        DM(0,0)=Lambda.at(0);
        DM(1,1)=Lambda.at(1);
        DM(2,2)=Lambda.at(2);

    }
    else
    {
        DM=M;
    }

    if (DM(0,0)>=-_Accuracy && DM(1,1)>=-_Accuracy && DM(2,2)>=-_Accuracy) _ExtrType=0;
    else _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=3;

    //    if (DM(0,0)>=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=1;
    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;

    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=2;

    for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++)
        {
            if (i!=j && fabs(DM(i,j))>0.001) _ExtrType=-1;
        }

    return DM;
}



QGenericMatrix<3,3,double> Data::IsMinimumSS(double T, Solutions sol)
{
    QVector<double> Lambda;
    QGenericMatrix<3,3,double> U;
    QGenericMatrix<3,3,double> U2;
    QGenericMatrix<3,3,double> U3;
    QGenericMatrix<3,3,double> M=FEnergySecondDerivative(T,sol);
    QGenericMatrix<3,3,double> DM;

    U.fill(0);
    U2.fill(0);

    double A11=M(0,0);
    double A12=M(0,1);
    double A13=M(0,2);

    double A22=M(1,1);
    double A23=M(1,2);
    double A33=M(2,2);

    Lambda=SolveQubicEquationSS(M);

    DM(0,0)=Lambda.at(0);
    DM(1,1)=Lambda.at(1);
    DM(2,2)=Lambda.at(2);


    if (DM(0,0)>=-_Accuracy && DM(1,1)>=-_Accuracy && DM(2,2)>=-_Accuracy) _ExtrType=0;
    else _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=3;

    //    if (DM(0,0)>=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=1;
    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)>=-0.001) _ExtrType=1;

    //    if (DM(0,0)>=-0.001 && DM(1,1)<=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)>=-0.001 && DM(2,2)<=-0.001) _ExtrType=2;
    //    if (DM(0,0)<=-0.001 && DM(1,1)<=-0.001 && DM(2,2)>=-0.001) _ExtrType=2;

    for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++)
        {
            if (i!=j && fabs(DM(i,j))>0.001) _ExtrType=-1;
        }

    return DM;

}









bool Data::IsEqual(double n, double T, Solutions sol)
{
    bool is_equal=false;

    if (_Solutions->size()!=0)
    {
        for (int i=0; i<_Solutions->size();i++)
        {
            double FEn1=FEnergy(n,T,sol);
            double FEn2=FEnergy(n,T,_Solutions->at(i));

            if (
                    fabs(FEn2-FEn1) < _Accuracy
                    &&
                    fabs(fabs(HzA(_Solutions->at(i)))-fabs(HzA(sol))) < _Accuracy
                    &&
                    fabs(fabs(H2(_Solutions->at(i)))-fabs(H2(sol))) < _Accuracy
                    &&
                    fabs(fabs(H2A(_Solutions->at(i)))-fabs(H2A(sol))) < _Accuracy
                    ) is_equal=true;
        }
    }

    if (
            fabs(EqN(n,T,sol)) > _Accuracy
            ||
            fabs(EqCO(n,T,sol)) > _Accuracy
            ||
            fabs(EqSF(n,T,sol)) > _Accuracy
            ||
            fabs(EqASF(n,T,sol)) > _Accuracy
            ) is_equal=true;

    if (qIsNaN(FEnergy(n,T,sol)) || qIsNaN(Hz(sol)) || fabs(Hz(sol)) > 40.0 || fabs(HzA(sol)) > 4.0 || fabs(H2(sol)) > 3.0 || fabs(H2A(sol)) > 3.0 ) is_equal=true;

    return is_equal;
}



void Data::FindSolutions(double n, double T)
{

    for (int i=0; i < _Sign->size(); i++)
    {
        Solutions sol;

        sol=RealSolve(n, T, _Sign->at(i));

        if (!IsEqual(n,T,sol)) _Solutions->append(sol);
    }

}


void Data::PrintSolutions(double n, double T, QFile &Lfile, QFile &Sfile)
{
    int FACU=12;

    QTextStream PrintLocal( &Lfile );
    QTextStream PrintSecond( &Sfile );

    QGenericMatrix<3,3,double> DM;

    Solutions NO;
    QVector<Solutions>* CO;
    QVector<Solutions>* SF;
    QVector<Solutions>* SS;

    CO = new QVector<Solutions>;
    SF = new QVector<Solutions>;
    SS = new QVector<Solutions>;

    for (int i=0; i<_Solutions->size();i++)
    {
        int Type1=0;

        if (fabs(HzA(_Solutions->at(i))) > 0.0001 && fabs(H2(_Solutions->at(i))) < 0.0001 && fabs(H2A(_Solutions->at(i))) < 0.0001)
        {
            Type1=1;
            CO->append(_Solutions->at(i));
        }

        if (fabs(HzA(_Solutions->at(i))) < 0.0001 && fabs(H2(_Solutions->at(i))) > 0.0001 && fabs(H2A(_Solutions->at(i))) < 0.0001)
        {
            Type1=2;
            SF->append(_Solutions->at(i));
        }

        if (fabs(HzA(_Solutions->at(i))) > 0.0001 && fabs(H2(_Solutions->at(i))) > 0.0001 && fabs(H2A(_Solutions->at(i))) > 0.0001)
        {
            Type1=3;
            SS->append(_Solutions->at(i));
        }

        if (Type1==0) NO=_Solutions->at(i);
    }

    DM=IsMinimumNO(T,NO);

    int Type=(_ExtrType==0 ? 1 : 0);

    int TypeNO=_ExtrType;

    PrintSecond << 0
                << "\t"
                << _ExtrType
                << "\t"
                << n
                << "\t"
                << T
                << "\t"
                << qSetRealNumberPrecision(FACU) << DM(0,0)
                << "\t"
                << qSetRealNumberPrecision(FACU) << DM(1,1)
                << "\t"
                << qSetRealNumberPrecision(FACU) << DM(2,2)
                << "\n";

    if ( !(fabs(_GlobalFEnergy-FEnergy(n,T,NO)) < _Accuracy)) Type*=-1;

    PrintLocal << 0
               << "\t"
               << Type
               << "\t"
               << TypeNO
               << "\t"
               << n
               << "\t"
               << T
               << "\t"
               << qSetRealNumberPrecision(FACU) << Hz(NO)
               << "\t"
               << qSetRealNumberPrecision(FACU) << HzA(NO)
               << "\t"
               << qSetRealNumberPrecision(FACU) << H2(NO)
               << "\t"
               << qSetRealNumberPrecision(FACU) << H2A(NO)
               << "\t"
               << qSetRealNumberPrecision(FACU) << FEnergy(n,T,NO)
               << "\n";


    qSort(*CO);
    qSort(*SF);
    qSort(*SS);

    Type=0;

    QVector<long double> FEn;
    QVector<int> FEn2;


    for (int i=0; i<CO->size();i++)
    {
        FEn.append(FEnergy(n,T,CO->at(i)));
        DM=IsMinimumCO(T,CO->at(i));
        FEn2.append(_ExtrType);

        PrintSecond << 1
                    << "\t"
                    << _ExtrType
                    << "\t"
                    << n
                    << "\t"
                    << T
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(0,0)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(1,1)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(2,2)
                    << "\n";
    }

    if (CO->size()==1)
    {
        bool is_ok = (TypeNO != 0 && FEn2.at(0) == 0);

        if (is_ok) if (FEn.at(0) <= FEnergy(n,T,NO)) Type=2;
    }

    if (CO->size()==2)
    {
        bool is_ok = (TypeNO == 0 && FEn2.at(0) != 0 && FEn2.at(1) == 0);

        if (is_ok) Type=( FEnergy(n,T,NO) <= FEn.at(1) ? 3 : 4 );
    }


    if (CO->size()==3)
    {
        bool is_ok = (TypeNO != 0 && FEn2.at(0) == 0 && FEn2.at(1) != 0 && FEn2.at(2) == 0);

        if (is_ok) Type=( FEn.at(0) <= FEn.at(2) ? 5 : 6 );
    }

    if (Type==0) for (int i=0; i<CO->size();i++) if (FEn2.at(i)==0) Type=1;

    if (FEn2.at(CO->size()-1)==0) if ( !(fabs(_GlobalFEnergy-FEn.at(CO->size()-1)) < _Accuracy)) Type*=-1;

    for (int i=0; i<CO->size();i++)
    {
        PrintLocal << 1
                   << "\t"
                   << Type
                   << "\t"
                   << FEn2.at(i)
                   << "\t"
                   << n
                   << "\t"
                   << T
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << Hz(CO->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << HzA(CO->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2(CO->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2A(CO->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << FEnergy(n,T,CO->at(i))
                   << "\n";
    }


    Type=0;

    FEn.clear();
    FEn2.clear();


    for (int i=0; i<SF->size();i++)
    {
        FEn.append(FEnergy(n,T,SF->at(i)));
        DM=IsMinimumSF(T,SF->at(i));
        FEn2.append(_ExtrType);

        PrintSecond << 2
                    << "\t"
                    << _ExtrType
                    << "\t"
                    << n
                    << "\t"
                    << T
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(0,0)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(1,1)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(2,2)
                    << "\n";
    }

    if (SF->size()==1)
    {
        bool is_ok = (TypeNO != 0 && FEn2.at(0) == 0);

        if (is_ok) if (FEn.at(0) <= FEnergy(n,T,NO)) Type=2;
    }

    if (SF->size()==2)
    {
        bool is_ok = (FEn2.at(0) == 0 && FEn2.at(1) == 0);

        if (is_ok) Type=3;
    }

    if (Type==0) for (int i=0; i<SF->size();i++) if (FEn2.at(i)==0) Type=1;

    if (Type==3)
    {
        if (!(fabs(_GlobalFEnergy-FEn.at(0)) < _Accuracy) && !(fabs(_GlobalFEnergy-FEn.at(1)) < _Accuracy)) Type*=-1;
    }
    else
    {
        for (int i=0; i<SF->size();i++)
        {
            if (FEn2.at(i)==0) if (!(fabs(_GlobalFEnergy-FEn.at(i)) < _Accuracy)) Type*=-1;
        }
    }

    for (int i=0; i<SF->size();i++)
    {
        PrintLocal << 2
                   << "\t"
                   << Type
                   << "\t"
                   << FEn2.at(i)
                   << "\t"
                   << n
                   << "\t"
                   << T
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << Hz(SF->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << HzA(SF->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2(SF->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2A(SF->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << FEnergy(n,T,SF->at(i))
                   << "\n";
    }



    Type=0;

    FEn.clear();
    FEn2.clear();



    for (int i=0; i<SS->size();i++)
    {
        FEn.append(FEnergy(n,T,SS->at(i)));
        DM=IsMinimumSS(T,SS->at(i));
        FEn2.append(_ExtrType);

        PrintSecond << 3
                    << "\t"
                    << _ExtrType
                    << "\t"
                    << n
                    << "\t"
                    << T
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(0,0)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(1,1)
                    << "\t"
                    << qSetRealNumberPrecision(FACU) << DM(2,2)
                    << "\n";
    }

    if (SS->size()==1)
    {
        bool is_ok = (TypeNO != 0 && FEn2.at(0) == 0);

        if (is_ok) if (FEn.at(0) <= FEnergy(n,T,NO)) Type=2;
    }

    if (SS->size()==2)
    {
        bool is_ok = (FEn2.at(0) == 0 && FEn2.at(1) == 0);

        if (is_ok) Type=3;
    }

    if (Type==0) for (int i=0; i<SS->size();i++) if (FEn2.at(i)==0) Type=1;

    if (Type==3)
    {
        if (!(fabs(_GlobalFEnergy-FEn.at(0)) < _Accuracy) && !(fabs(_GlobalFEnergy-FEn.at(1)) < _Accuracy)) Type*=-1;
    }
    else
    {
        for (int i=0; i<SS->size();i++)
        {
            if (FEn2.at(i)==0) if (!(fabs(_GlobalFEnergy-FEn.at(i)) < _Accuracy)) Type*=-1;
        }
    }

    for (int i=0; i<SS->size();i++)
    {

        PrintLocal << 3
                   << "\t"
                   << Type
                   << "\t"
                   << FEn2.at(i)
                   << "\t"
                   << n
                   << "\t"
                   << T
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << Hz(SS->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << HzA(SS->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2(SS->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << H2A(SS->at(i))
                   << "\t"
                   << qSetRealNumberPrecision(FACU) << FEnergy(n,T,SS->at(i))
                   << "\n";
    }


    delete CO;
    delete SF;
    delete SS;

}



void Data::PrintGlobalSolutions(double n, double T, QFile &Gfile)
{
    int FACU=12;

    QTextStream PrintGlobal( &Gfile );

    QVector<double> FEn;

    QGenericMatrix<3,3,double> DM;

    for (int i=0; i<_Solutions->size();i++)
        FEn.append(FEnergy(n,T,_Solutions->at(i)));

    qSort(FEn);

    _GlobalFEnergy=FEn.at(0);

    for (int i=0; i<_Solutions->size();i++)
    {
        if (fabs(FEnergy(n,T,_Solutions->at(i))-FEn.at(0)) < _Accuracy)
        {
            int Type=0;

            if (fabs(HzA(_Solutions->at(i))) > 0.0001 && fabs(H2(_Solutions->at(i))) < 0.0001 && fabs(H2A(_Solutions->at(i))) < 0.0001) Type=1;
            if (fabs(HzA(_Solutions->at(i))) < 0.0001 && fabs(H2(_Solutions->at(i))) > 0.0001 && fabs(H2A(_Solutions->at(i))) < 0.0001) Type=2;
            if (fabs(HzA(_Solutions->at(i))) > 0.0001 && fabs(H2(_Solutions->at(i))) > 0.0001 && fabs(H2A(_Solutions->at(i))) > 0.0001) Type=3;

            if (Type==0) DM=IsMinimumNO(T,_Solutions->at(i));
            if (Type==1) DM=IsMinimumCO(T,_Solutions->at(i));
            if (Type==2) DM=IsMinimumSF(T,_Solutions->at(i));
            if (Type==3) DM=IsMinimumSS(T,_Solutions->at(i));

            if (_ExtrType==0)
            {
                PrintGlobal << Type
                            << "\t"
                            << T
                            << "\t"
                            << n
                            << "\t"
                            << qSetRealNumberPrecision(FACU) << FEnergy(n,T,_Solutions->at(i))
                            << "\n";
            }
        }
    }

}



double Data::Norm(double& x1, double& y1, double& z1, double& w1, double& x2, double& y2, double& z2, double& w2)
{
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)+(w2-w1)*(w2-w1));
}

QGenericMatrix<4,4,double> Data::FindInverse(QGenericMatrix<4,4,double>& M)
{
    double A11=M(0,0);
    double A12=M(0,1);
    double A13=M(0,2);
    double A14=M(0,3);

    double A21=M(1,0);
    double A22=M(1,1);
    double A23=M(1,2);
    double A24=M(1,3);

    double A31=M(2,0);
    double A32=M(2,1);
    double A33=M(2,2);
    double A34=M(2,3);

    double A41=M(3,0);
    double A42=M(3,1);
    double A43=M(3,2);
    double A44=M(3,3);

    double d1=A33*A44-A43*A34;
    double d2=A32*A44-A42*A34;
    double d3=A32*A43-A42*A33;
    double d4=A31*A43-A41*A33;
    double d5=A31*A42-A41*A32;
    double d6=A31*A44-A41*A34;

    double dd1=A23*A44-A43*A24;
    double dd2=A22*A44-A42*A24;
    double dd3=A22*A43-A42*A23;
    double dd4=A21*A44-A41*A24;
    double dd5=A21*A43-A41*A23;
    double dd6=A21*A42-A41*A22;

    double ddd1=A23*A34-A33*A24;
    double ddd2=A22*A34-A32*A24;
    double ddd3=A22*A33-A32*A23;
    double ddd4=A21*A34-A31*A24;
    double ddd5=A21*A33-A31*A23;
    double ddd6=A21*A32-A31*A22;



    double d11=A22*d1-A23*d2+A24*d3;
    double d12=A21*d1-A23*d6+A24*d4;
    double d13=A21*d2-A22*d6+A24*d5;
    double d14=A21*d3-A22*d4+A23*d5;

    double d21=A12*d1-A13*d2+A14*d3;
    double d22=A11*d1-A13*d6+A14*d4;
    double d23=A11*d2-A12*d6+A14*d5;
    double d24=A11*d3-A12*d4+A13*d5;

    double d31=A12*dd1-A13*dd2+A14*dd3;
    double d32=A11*dd1-A13*dd4+A14*dd5;
    double d33=A11*dd2-A12*dd4+A14*dd6;
    double d34=A11*dd3-A12*dd5+A13*dd6;

    double d41=A12*ddd1-A13*ddd2+A14*ddd3;
    double d42=A11*ddd1-A13*ddd4+A14*ddd5;
    double d43=A11*ddd2-A12*ddd4+A14*ddd6;
    double d44=A11*ddd3-A12*ddd5+A13*ddd6;

    double d=A11*d11-A12*d12+A13*d13-A14*d14;

    QGenericMatrix<4,4,double> InvM;

    InvM.fill(0);

    InvM(0,0)=d11/d;
    InvM(0,1)=-d21/d;
    InvM(0,2)=d31/d;
    InvM(0,3)=-d41/d;

    InvM(1,0)=-d12/d;
    InvM(1,1)=d22/d;
    InvM(1,2)=-d32/d;
    InvM(1,3)=d42/d;

    InvM(2,0)=d13/d;
    InvM(2,1)=-d23/d;
    InvM(2,2)=d33/d;
    InvM(2,3)=-d43/d;

    InvM(3,0)=-d14/d;
    InvM(3,1)=d24/d;
    InvM(3,2)=-d34/d;
    InvM(3,3)=d44/d;

    return InvM;

}

Solutions Data::RealSolve(double &n, double &T, Solutions sol)
{
    Solutions point=sol;

    QGenericMatrix<1,4,double> X0;
    QGenericMatrix<1,4,double> Xi;
    QGenericMatrix<1,4,double> F;
    QGenericMatrix<4,4,double> J;
    QGenericMatrix<4,4,double> InvJ;
    QGenericMatrix<4,4,double> Test;


    X0.fill(0);
    Xi.fill(1);
    F.fill(0);
    J.fill(0);
    InvJ.fill(0);
    Test.fill(0);

    int count=0;

    while( Norm(X0(0,0),X0(1,0),X0(2,0),X0(3,0),Xi(0,0),Xi(1,0),Xi(2,0),Xi(3,0)) > _Accuracy )
    {
        X0(0,0)=Hz(point);
        X0(1,0)=HzA(point);
        X0(2,0)=H2(point);
        X0(3,0)=H2A(point);

        F(0,0)=EqN(n, T, point);
        F(1,0)=EqCO(n, T, point);
        F(2,0)=EqSF(n, T, point);
        F(3,0)=EqASF(n, T, point);


        J(0,0)=-0.5*exp(-_D/T)*(PiP(T,point)+PiM(T,point))/T;
        J(0,1)=-0.5*exp(-_D/T)*(PiP(T,point)-PiM(T,point))/T;
        J(0,2)=-0.5*exp(-_D/T)*(PsiP(T,point)+PsiM(T,point))/T;
        J(0,3)=-0.5*exp(-_D/T)*(PsiP(T,point)-PsiM(T,point))/T;


        J(1,0)=-2*_V*exp(-_D/T)*(PiP(T,point)-PiM(T,point))/T;
        J(1,1)=1-2*_V*exp(-_D/T)*(PiP(T,point)+PiM(T,point))/T;
        J(1,2)=-2*_V*exp(-_D/T)*(PsiP(T,point)-PsiM(T,point))/T;
        J(1,3)=-2*_V*exp(-_D/T)*(PsiP(T,point)+PsiM(T,point))/T;

        J(2,0)=-_tb*exp(-_D/T)*(PsiP(T,point)+PsiM(T,point))/T;
        J(2,1)=-_tb*exp(-_D/T)*(PsiP(T,point)-PsiM(T,point))/T;
        J(2,2)=1-_tb*exp(-_D/T)*(DeltaP(T,point)+DeltaM(T,point))/T;
        J(2,3)=-_tb*exp(-_D/T)*(DeltaP(T,point)-DeltaM(T,point))/T;

        J(3,0)=_tb*exp(-_D/T)*(PsiP(T,point)-PsiM(T,point))/T;
        J(3,1)=_tb*exp(-_D/T)*(PsiP(T,point)+PsiM(T,point))/T;
        J(3,2)=_tb*exp(-_D/T)*(DeltaP(T,point)-DeltaM(T,point))/T;
        J(3,3)=1+_tb*exp(-_D/T)*(DeltaP(T,point)+DeltaM(T,point))/T;


        InvJ=FindInverse(J);

        Test=J*InvJ;

        Xi=X0-InvJ*F;

        SetHz(point,Xi(0,0));
        SetHzA(point,Xi(1,0));
        SetH2(point,Xi(2,0));
        SetH2A(point,Xi(3,0));

        count++;

        if ( count==10000 || fabs(Xi(0,0)) > 40 || fabs(Xi(1,0)) > 4 || fabs(Xi(2,0)) > 3 || fabs(Xi(3,0)) > 3 )
        {
            X0(0,0)=Xi(0,0);
            X0(1,0)=Xi(1,0);
            X0(2,0)=Xi(2,0);
            X0(3,0)=Xi(3,0);
        }

    }

    SetHz(point,Xi(0,0));
    SetHzA(point,fabs(Xi(1,0)));
    SetH2(point,fabs(Xi(2,0)));
    SetH2A(point,fabs(Xi(3,0)));

    return point;

}



void Data::MiniSign(double &n,double &T,Solutions point1,Solutions point2)
{
    bool is_equal=false;
    double signconst=0.01;

    if (
            (
                EqN(n,T,point1)*EqN(n,T,point2)<signconst
                &&
                EqCO(n,T,point1)*EqCO(n,T,point2)<signconst
                &&
                EqSF(n,T,point1)*EqSF(n,T,point2)<=signconst
                &&
                EqASF(n,T,point1)*EqASF(n,T,point2)<=signconst
                )
            ||
            (
                EqN(n,T,point1)*EqN(n,T,point2)<signconst
                &&
                EqCO(n,T,point1)*EqCO(n,T,point2)<=signconst
                &&
                EqSF(n,T,point1)*EqSF(n,T,point2)<signconst
                &&
                EqASF(n,T,point1)*EqASF(n,T,point2)<=signconst
                )
            ||
            (
                EqN(n,T,point1)*EqN(n,T,point2)<signconst
                &&
                EqCO(n,T,point1)*EqCO(n,T,point2)<signconst
                &&
                EqSF(n,T,point1)*EqSF(n,T,point2)<signconst
                &&
                EqASF(n,T,point1)*EqASF(n,T,point2)<signconst
                )
            )
    {
        Solutions point=point1+point2;

        SetHz(point,0.5*Hz(point));
        SetHzA(point,0.5*HzA(point));
        SetH2(point,0.5*H2(point));
        SetH2A(point,0.5*H2A(point));

        if (_Sign->size()!=0)
        {
            for (int i=0; i<_Sign->size();i++)
            {
                double FEn1=FEnergy(n,T,point);
                double FEn2=FEnergy(n,T,_Sign->at(i));

                if (
                        fabs(FEn2-FEn1) < 0.01 // 0.1
                        &&
                        fabs(fabs(HzA(_Sign->at(i)))-fabs(HzA(point))) < 0.01 //0.1 //0.05
                        &&
                        fabs(fabs(H2(_Sign->at(i)))-fabs(H2(point))) < 0.01 //0.05
                        &&
                        fabs(fabs(H2A(_Sign->at(i)))-fabs(H2A(point))) < 0.01 //0.05

                        ) is_equal=true;
            }
        }

        if (!is_equal) _Sign->append(point);
    }
}


void Data::Sign(double &n, double &T, Solutions Min, Solutions Max)
{
    Solutions point1;
    Solutions point2;

    double SStep=10.0;

    double StepHz=(Hz(Max)-Hz(Min))/SStep;
    double StepHzA=(HzA(Max)-HzA(Min))/SStep;
    double StepH2=(H2(Max)-H2(Min))/SStep;
    double StepH2A=(H2A(Max)-H2A(Min))/SStep;


    for (double X1=Hz(Min);X1 < Hz(Max); X1+=StepHz)
        for (double X2=HzA(Min);X2 < HzA(Max); X2+=StepHzA)
            for (double X3=H2(Min);X3 < H2(Max); X3+=StepH2)
                for (double X4=H2A(Min);X4 < H2A(Max); X4+=StepH2A)
                {

                    SetHz(point1,X1);
                    SetHzA(point1,X2);
                    SetH2(point1,X3);
                    SetH2A(point1,X4);

                    if(SqP(T,point1)!=0.0 && SqM(T,point1)!=0.0)
                    {
                        SetHz(point2,X1+StepHz);
                        SetHzA(point2,X2);
                        SetH2(point2,X3);
                        SetH2A(point2,X4);

                        MiniSign(n,T,point1,point2);

                        SetHz(point2,X1);
                        SetHzA(point2,X2+StepHzA);
                        SetH2(point2,X3);
                        SetH2A(point2,X4);

                        MiniSign(n,T,point1,point2);

                        SetHz(point2,X1+StepHz);
                        SetHzA(point2,X2+StepHzA);
                        SetH2(point2,X3);
                        SetH2A(point2,X4);

                        MiniSign(n,T,point1,point2);


                        SetHz(point2,X1);
                        SetHzA(point2,X2);
                        SetH2(point2,X3+StepH2);
                        SetH2A(point2,X4);

                        MiniSign(n,T,point1,point2);

                        SetHz(point2,X1+StepHz);
                        SetHzA(point2,X2);
                        SetH2(point2,X3+StepH2);
                        SetH2A(point2,X4);

                        MiniSign(n,T,point1,point2);

                        SetHz(point2,X1);
                        SetHzA(point2,X2+StepHzA);
                        SetH2(point2,X3+StepH2);
                        SetH2A(point2,X4+StepH2A);

                        MiniSign(n,T,point1,point2);

                        SetHz(point2,X1+StepHz);
                        SetHzA(point2,X2+StepHzA);
                        SetH2(point2,X3+StepH2);
                        SetH2A(point2,X4+StepH2A);

                        MiniSign(n,T,point1,point2);
                    }

                }
}
