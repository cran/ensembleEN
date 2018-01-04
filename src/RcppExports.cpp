// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Ensemble_EN_Objective
double Ensemble_EN_Objective(const arma::mat& current_res, const arma::mat& beta, const double& lambda_sparsity, const double& lambda_diversity, const double& alpha);
RcppExport SEXP _ensembleEN_Ensemble_EN_Objective(SEXP current_resSEXP, SEXP betaSEXP, SEXP lambda_sparsitySEXP, SEXP lambda_diversitySEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type current_res(current_resSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_sparsity(lambda_sparsitySEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_diversity(lambda_diversitySEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(Ensemble_EN_Objective(current_res, beta, lambda_sparsity, lambda_diversity, alpha));
    return rcpp_result_gen;
END_RCPP
}
// Ensemble_EN_Grid
arma::cube Ensemble_EN_Grid(const arma::mat& x, const arma::vec& y, const int& which_lambda, const arma::vec& lambdas_grid, const double& lambda_fixed, const double& alpha, const arma::uword& num_groups, const double& tolerance, const arma::uword& max_iter);
RcppExport SEXP _ensembleEN_Ensemble_EN_Grid(SEXP xSEXP, SEXP ySEXP, SEXP which_lambdaSEXP, SEXP lambdas_gridSEXP, SEXP lambda_fixedSEXP, SEXP alphaSEXP, SEXP num_groupsSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type which_lambda(which_lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambdas_grid(lambdas_gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_fixed(lambda_fixedSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_groups(num_groupsSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Ensemble_EN_Grid(x, y, which_lambda, lambdas_grid, lambda_fixed, alpha, num_groups, tolerance, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// Main_Ensemble_EN
List Main_Ensemble_EN(const arma::mat& x_perm, const arma::vec& y_perm, const arma::uword num_lambdas_sparsity, const arma::uword num_lambdas_diversity, const double& alpha, const arma::uword& num_groups, const double& tolerance, const arma::uword& max_iter, const arma::uword& num_folds, const arma::uword& num_threads);
RcppExport SEXP _ensembleEN_Main_Ensemble_EN(SEXP x_permSEXP, SEXP y_permSEXP, SEXP num_lambdas_sparsitySEXP, SEXP num_lambdas_diversitySEXP, SEXP alphaSEXP, SEXP num_groupsSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP, SEXP num_foldsSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x_perm(x_permSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_perm(y_permSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type num_lambdas_sparsity(num_lambdas_sparsitySEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type num_lambdas_diversity(num_lambdas_diversitySEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_groups(num_groupsSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_folds(num_foldsSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Main_Ensemble_EN(x_perm, y_perm, num_lambdas_sparsity, num_lambdas_diversity, alpha, num_groups, tolerance, max_iter, num_folds, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ensembleEN_Ensemble_EN_Objective", (DL_FUNC) &_ensembleEN_Ensemble_EN_Objective, 5},
    {"_ensembleEN_Ensemble_EN_Grid", (DL_FUNC) &_ensembleEN_Ensemble_EN_Grid, 9},
    {"_ensembleEN_Main_Ensemble_EN", (DL_FUNC) &_ensembleEN_Main_Ensemble_EN, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_ensembleEN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
