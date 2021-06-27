#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using Rcpp::as;
using Eigen::EigenSolver;
using Eigen::sqrt;
using Eigen::LLT;
using Eigen::Map;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::DiagonalMatrix;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::BDCSVD;
using Eigen::ColPivHouseholderQR;

// [[Rcpp::export]]
Eigen::MatrixXd multipleCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y) {
  Eigen::MatrixXd Z = X * Y;
  return Z;
}

// [[Rcpp::export]]
Eigen::MatrixXd crossprodCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y) {
  Eigen::MatrixXd Z = X.adjoint() * Y;
  return Z;
}

// [[Rcpp::export]]
Eigen::MatrixXd tcrossprodCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y) {
  Eigen::MatrixXd Z = X * Y.adjoint();
  return Z;
}

// [[Rcpp::export]]
Eigen::MatrixXd selfcrossprodCpp(Eigen::Map<Eigen::MatrixXd> & A) {
  const int n(A.cols());
  Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
  return AtA;
}

// [[Rcpp::export]]
Eigen::MatrixXd selftcrossprodCpp(Eigen::Map<Eigen::MatrixXd> & A) {
  const int m(A.rows());
  Eigen::MatrixXd AAt(Eigen::MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(A));
  return AAt;
}

// [[Rcpp::export]]
List svdCpp(Eigen::Map<Eigen::MatrixXd> & X){
  Eigen::BDCSVD<MatrixXd> se(X, Eigen::ComputeThinU|Eigen::ComputeThinV);
  return List::create(Named("u") = se.matrixU(), Named("v") = se.matrixV(), Named("d") = se.singularValues());
}

// [[Rcpp::export]]
Eigen::MatrixXd QRsolveCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y){
  Eigen::MatrixXd se(X.colPivHouseholderQr().solve(Y));
  return se;
}

// [[Rcpp::export]]
Eigen::VectorXd lmCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::VectorXd> y){
  const int n(X.cols());
  LLT<MatrixXd> Z(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
  Eigen::VectorXd coef = Z.solve(X.adjoint()*y);
  return coef;
}

// [[Rcpp::export]]
List simplsCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y, int ncomp){
  std::list<Eigen::MatrixXd> coef;
  const int n = X.rows(), p = X.cols(), q = Y.cols();
  Eigen::MatrixXd V(Eigen::MatrixXd(p, ncomp).setZero()), R(Eigen::MatrixXd(p, ncomp)), tQ(Eigen::MatrixXd(ncomp, q));
  Eigen::MatrixXd S = X.adjoint() * Y, StS, SSt;
  Eigen::VectorXd qal(q), qan(q), pan(p), ra(p), ta(n), taj(n), van(p);
  Eigen::MatrixXd B(p, q);
  for(int i = 0; i < ncomp; i++){
    StS = Eigen::MatrixXd(q, q).setZero().selfadjointView<Eigen::Lower>().rankUpdate(S.adjoint());
    SSt = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(S);
    if(p > q) {
      Eigen::EigenSolver<Eigen::MatrixXd> qa(StS);
      qan = qa.pseudoEigenvectors().col(0);
    } else {
      Eigen::EigenSolver<Eigen::MatrixXd> qa(SSt);
      qal = S.adjoint() * qa.pseudoEigenvectors().col(0);
      qan = qal / sqrt(qal.adjoint() * qal);
    }
    ra = S * qan;
    ta = X * ra;
    for(int j = 0; j < n; j++){
      taj[j] = ta[j] - ta.mean();
    }
    ta = taj / sqrt(taj.adjoint() * taj);
    ra = ra / sqrt(taj.adjoint() * taj);
    pan = X.adjoint() * ta;
    qan = Y.adjoint() * ta;
    van = pan - V * (V.adjoint() * pan);
    van = van / sqrt(van.adjoint() * van);
    S = S - van * (van.adjoint() * S);
    R.col(i) = ra;
    tQ.row(i) = qan;
    V.col(i) = van;
    B = R * tQ;
    coef.push_back(B);
  }
  return List::create(Named("B") = B, Named("Coef") = coef);
}

// [[Rcpp::export]]
List svdpcrCpp(Eigen::Map<Eigen::MatrixXd> & X, Eigen::Map<Eigen::MatrixXd> & Y, int ncomp){
  std::list<Eigen::MatrixXd> coef;
  int n = X.rows(), p = X.cols(), q = Y.cols();
  Eigen::BDCSVD<MatrixXd> se(X, Eigen::ComputeFullU|Eigen::ComputeFullV);
  Eigen::VectorXd D = se.singularValues().head(ncomp);
  Eigen::MatrixXd D0 = D.asDiagonal();
  Eigen::MatrixXd TT = se.matrixU().block(0, 0, n, ncomp) * D0;
  Eigen::MatrixXd PP = se.matrixV().block(0, 0, p, ncomp);
  Eigen::MatrixXd QQ = TT.adjoint() * Y;
  Eigen::ArrayXXd tQ = QQ.array().colwise() / D.array();
  Eigen::MatrixXd B(p, q);
  B = PP * tQ.matrix();
  return List::create(Named("D") = D0, Named("TT") = TT, Named("PP") = PP, Named("tQ") = tQ, Named("QQ") = QQ);
}




