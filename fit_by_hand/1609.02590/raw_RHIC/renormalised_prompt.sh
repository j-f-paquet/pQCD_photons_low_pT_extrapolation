for i in pp pAu dAu He3Au AuAu; do echo "#pT (GeV) E dsigma/d^3p (GeV^-2 mb^2)" > ../${i}_RHIC_200GeV.dat; done

tail +2 Qs40/pp.dat    |  awk '{print $1,$2*1.81*1e-9}'  >> ../pp_RHIC_200GeV.dat
tail +2 Qs40/pAu.dat   |  awk '{print $1,$2*1.87*1e-9}'  >> ../pAu_RHIC_200GeV.dat
tail +2 Qs40/dAu.dat   |  awk '{print $1,$2*1.88*1e-9}'  >> ../dAu_RHIC_200GeV.dat
tail +2 Qs40/He3Au.dat |  awk '{print $1,$2*1.88*1e-9}'  >> ../He3Au_RHIC_200GeV.dat
tail +2 Qs40/AuAu.dat  |  awk '{print $1,$2*1.90*1e-9}'  >> ../AuAu_RHIC_200GeV.dat

