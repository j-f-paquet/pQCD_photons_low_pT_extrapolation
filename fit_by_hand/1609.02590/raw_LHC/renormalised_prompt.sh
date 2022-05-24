for i in pp pPb; do echo "#pT (GeV) E dsigma/d^3p (GeV^-2 mb^2)" > ../${i}_LHC_5020GeV.dat; done

tail +2 Qs40/pp.dat | awk '{print $1,$2*1.22*1e-9}' >> ../pp_LHC_5020GeV.dat
tail +2 Qs40/pPb.dat | awk '{print $1,$2*1.19*1e-9}' >> ../pPb_LHC_5020GeV.dat
