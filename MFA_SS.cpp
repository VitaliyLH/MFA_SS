#include "MFA_SS.h"

Tc::Tc(QFile &file, QFile &Gfile, QFile &Lfile, QFile &Sfile)
{
   data.parse(file, Gfile, Lfile, Sfile);
}


Tc::~Tc()
{
  data.~Data();
}
