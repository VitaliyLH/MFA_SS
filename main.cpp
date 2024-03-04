#include"MFA_SS.h"
#include <QStringBuilder>
#include <iomanip>


int main(int argc, char *argv[])
{
    if (argc != 2) {
        puts("WTF!");
        return 1;
    }

    Tc* T;

    QFile file(argv[1]);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return EXIT_FAILURE;

    QString G("Global");
    QString L("LocalAndSaddle");
    QString S("SecondDerivatives");

    QString other(argv[1]);

    int index = other.indexOf(".dat");

    QString fileNameG(other);
    QString fileNameL(other);
    QString fileNameS(other);

    fileNameG.insert(index,G);
    fileNameL.insert(index,L);
    fileNameS.insert(index,S);

    QFile Gfile(fileNameG);
    QFile Lfile(fileNameL);
    QFile Sfile(fileNameS);

    Gfile.open(QIODevice::WriteOnly | QIODevice::Text);
    Lfile.open(QIODevice::WriteOnly | QIODevice::Text);
    Sfile.open(QIODevice::WriteOnly | QIODevice::Text);

    T = new Tc(file , Gfile , Lfile, Sfile);

    file.close();

    Gfile.close();
    Lfile.close();
    Sfile.close();

    return 0;

}
