#include "Rivet/Particle.hh"
#include "Rivet/Projections/Beam.hh"

int main() {
  using namespace Rivet;
  using namespace std;

  const FourMomentum pbeam1a(sqrt(1E6 + 1)*GeV, 0, 0, 1*TeV), pbeam1b(sqrt(1E6 + 1)*GeV, 0, 0, -1*TeV);
  const Particle beam1a(PID::PROTON, pbeam1a), beam1b(PID::PROTON, pbeam1b);
  //
  const Vector3 beta1a = cmsBetaVec(pbeam1a, pbeam1b);
  const Vector3 beta1b = cmsBetaVec(beam1a, beam1b);
  cout << "Beta_symm = " << beta1a << " or " << beta1b << endl;
  //
  const Vector3 gamma1a = cmsGammaVec(pbeam1a, pbeam1b);
  const Vector3 gamma1b = cmsGammaVec(beam1a, beam1b);
  cout << "Gamma_symm = " << gamma1a << " or " << gamma1b << endl;
  //
  const LorentzTransform trfb1 = LorentzTransform::mkFrameTransformFromBeta(cmsBetaVec(pbeam1a, pbeam1b));
  cout << "Beta trf matrix = " << trfb1 << endl;
  const FourMomentum pbeam1c = trfb1.transform(pbeam1a);
  const FourMomentum pbeam1d = trfb1.transform(pbeam1b);
  cout << "Beta-boosted_symm = " << pbeam1c << " + " << pbeam1d << " = " << (pbeam1c + pbeam1d) << endl;
  //
  //const LorentzTransform trfy1 = LorentzTransform::mkFrameTransformFromGamma(cmsGammaVec(pbeam1a, pbeam1b));
  const LorentzTransform trfy1 = cmsTransform(pbeam1a, pbeam1b);
  cout << "Gamma trf matrix = " << trfy1 << endl;
  const FourMomentum pbeam1e = trfy1.transform(pbeam1a);
  const FourMomentum pbeam1f = trfy1.transform(pbeam1b);
  cout << "Gamma-boosted_symm = " << pbeam1e << " + " << pbeam1f << " = " << (pbeam1e + pbeam1f) << endl;


  cout << endl;


  const FourMomentum pbeam2a(sqrt(101)*GeV, 0, 0, 10*GeV), pbeam2b(2*GeV, 0, 0, 0);
  // const Particle beam2a(, pbeam1), beam2b(, pbeam1);
  cout << "Original_asymm = " << pbeam2a << " + " << pbeam2b << " = " << (pbeam2a + pbeam2b) << endl;
  //
  const Vector3 beta2a = cmsBetaVec(pbeam2a, pbeam2b);
  //const Vector3 beta2b = cmsBetaVec(beam2a, beam2b);
  cout << "Beta_asymm = " << beta2a << endl; // << " or " << beta2b << endl;
  //
  const Vector3 gamma2a = cmsGammaVec(pbeam2a, pbeam2b);
  //const Vector3 gamma2b = cmsGammaVec(beam2a, beam2b);
  cout << "Gamma_asymm = " << gamma2a << endl; // << " or " << gamma2b << endl;
  cout << "Gamma_asymm2 = " << (pbeam2a+pbeam2b).gammaVec() << endl; // << " or " << gamma2b << endl;
  //
  const LorentzTransform trfb2 = LorentzTransform::mkFrameTransformFromBeta(cmsBetaVec(pbeam2a, pbeam2b));
  cout << "Beta trf matrix = " << trfb2 << endl;
  const FourMomentum pbeam2c = trfb2.transform(pbeam2a);
  const FourMomentum pbeam2d = trfb2.transform(pbeam2b);
  cout << "Beta-boosted_asymm = " << pbeam2c << " + " << pbeam2d << " = " << (pbeam2c + pbeam2d) << endl;
  //
  //const LorentzTransform trfy2 = LorentzTransform::mkFrameTransformFromGamma(cmsGammaVec(pbeam2a, pbeam2b));
  const LorentzTransform trfy2 = cmsTransform(pbeam2a, pbeam2b);
  cout << "Gamma trf matrix = " << trfy2 << endl;
  const FourMomentum pbeam2e = trfy2.transform(pbeam2a);
  const FourMomentum pbeam2f = trfy2.transform(pbeam2b);
  cout << "Gamma-boosted_asymm = " << pbeam2e << " + " << pbeam2f << " = " << (pbeam2e + pbeam2f) << endl;

  return 0;
}
