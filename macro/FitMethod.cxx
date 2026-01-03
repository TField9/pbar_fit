/**
 * @brief Exponential-Gaussian-Exponential function.
 * 
 * This function calculates the value of an Exponential-Gaussian-Exponential function
 * given the input parameters and the variable x.
 * 
 * @param x Pointer to the variable x.
 * @param p Pointer to the array of parameters.
 * @return double The value of the function.
 */
double funcExpGausExp(double *x, double *p);

/**
 * @class expgausexp
 * @brief RooFit PDF for Exponential-Gaussian-Exponential distribution.
 * 
 * This class implements a custom PDF for the Exponential-Gaussian-Exponential distribution
 * in the RooFit framework.
 */
class expgausexp : public RooAbsPdf {
public:
    /**
     * @brief Default constructor.
     */
    expgausexp() {}

    /**
     * @brief Constructor with parameters.
     * 
     * @param name Name of the PDF.
     * @param title Title of the PDF.
     * @param _x Variable x.
     * @param _x0 Parameter x0.
     * @param _sigmaL Parameter sigmaL.
     * @param _alphaL Parameter alphaL.
     * @param _sigmaR Parameter sigmaR.
     * @param _alphaR Parameter alphaR.
     */
    expgausexp(const char *name, const char *title,
               RooAbsReal& _x,
               RooAbsReal& _x0,
               RooAbsReal& _sigmaL,
               RooAbsReal& _alphaL,
               RooAbsReal& _sigmaR,
               RooAbsReal& _alphaR);

    /**
     * @brief Copy constructor.
     * 
     * @param other The object to copy.
     * @param name Name of the new object.
     */
    expgausexp(const expgausexp& other, const char* name = 0);

    /**
     * @brief Clone method.
     * 
     * @param newname Name of the cloned object.
     * @return TObject* Pointer to the cloned object.
     */
    virtual TObject* clone(const char* newname) const;

    /**
     * @brief Destructor.
     */
    inline virtual ~expgausexp() {}

    /**
     * @brief Get analytical integral.
     * 
     * @param allVars All variables.
     * @param analVars Analytical variables.
     * @param rangeName Range name.
     * @return Int_t Integral code.
     */
    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;

    /**
     * @brief Calculate analytical integral.
     * 
     * @param code Integral code.
     * @param rangeName Range name.
     * @return Double_t Integral value.
     */
    Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

protected:
    RooRealProxy x;       ///< Variable x.
    RooRealProxy x0;      ///< Parameter x0.
    RooRealProxy sigmaL;  ///< Parameter sigmaL.
    RooRealProxy alphaL;  ///< Parameter alphaL.
    RooRealProxy sigmaR;  ///< Parameter sigmaR.
    RooRealProxy alphaR;  ///< Parameter alphaR.

    /**
     * @brief Evaluate the PDF.
     * 
     * @return Double_t PDF value.
     */
    Double_t evaluate() const;
};

/**
 * @brief Minimizing function for combined NLL.
 * 
 * This function minimizes the combined negative log-likelihood (NLL) for signal, background,
 * and negative data histograms using the RooFit framework.
 * 
 * @param hsig Pointer to the signal histogram.
 * @param hbkg Pointer to the background histogram.
 * @param hneg Pointer to the negative data histogram.
 * @return std::vector<double> Vector of optimized parameters and their errors.
 */
std::vector<double> RooCombFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg);

/**
 * @brief Iterative Weighted Least Squares (IWLS) function fitting.
 * 
 * This function performs IWLS fitting on signal, background, and negative data histograms.
 * 
 * @param hsig Pointer to the signal histogram.
 * @param hbkg Pointer to the background histogram.
 * @param hneg Pointer to the negative data histogram.
 * @param nsig Reference to the number of signal events.
 * @param nsig_error Reference to the error in the number of signal events.
 * @param iteration Number of iterations for the fitting process.
 */
void IWLS_FunctionFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig, double &nsig_error, Int_t iteration);

/**
 * @class QMinimizer
 * @brief Class for Q optimization.
 * 
 * This class implements a Q optimization algorithm for signal, background, and observed data histograms.
 */
class QMinimizer : public ROOT::Math::IBaseFunctionMultiDim {
public:
    /**
     * @brief Constructor.
     * 
     * @param signal Pointer to the signal histogram.
     * @param background Pointer to the background histogram.
     * @param observed Pointer to the observed data histogram.
     */
    QMinimizer(TH1* signal, TH1* background, TH1* observed);

    /**
     * @brief Destructor.
     */
    ~QMinimizer() override;

    /**
     * @brief Clone method.
     * 
     * @return QMinimizer* Pointer to the cloned object.
     */
    QMinimizer* Clone() const override;

    /**
     * @brief Get the number of dimensions.
     * 
     * @return unsigned int Number of dimensions.
     */
    unsigned int NDim() const override;

    /**
     * @brief Evaluate the function.
     * 
     * @param p Pointer to the parameters.
     * @return double Function value.
     */
    double DoEval(const double* p) const override;

    /**
     * @brief Run the optimization.
     */
    void Minimize();

    /**
     * @brief Get the optimized y1 parameter.
     * 
     * @return double Optimized y1 parameter.
     */
    double GetY1() const;

    /**
     * @brief Get the optimized y2 parameter.
     * 
     * @return double Optimized y2 parameter.
     */
    double GetY2() const;

    /**
     * @brief Get the error in the optimized y1 parameter.
     * 
     * @return double Error in the optimized y1 parameter.
     */
    double GetY1Error() const;

    /**
     * @brief Get the error in the optimized y2 parameter.
     * 
     * @return double Error in the optimized y2 parameter.
     */
    double GetY2Error() const;

private:
    TH1* sig;  ///< Pointer to the signal histogram.
    TH1* bkg;  ///< Pointer to the background histogram.
    TH1* data; ///< Pointer to the observed data histogram.
    double y1_opt, y2_opt; ///< Optimized parameters.
    double y1_err, y2_err; ///< Parameter errors.

    /**
     * @brief Calculate the Q value.
     * 
     * @param y1 Parameter y1.
     * @param y2 Parameter y2.
     * @return double Q value.
     */
    double calculateQ(double y1, double y2) const;

    /**
     * @brief Calculate the Cash statistic.
     * 
     * @param observed Observed value.
     * @param expected Expected value.
     * @return double Cash statistic.
     */
    double cash(double observed, double expected) const;
};

/**
 * @brief Run Q optimization.
 * 
 * This function runs the Q optimization algorithm on signal, background, and observed data histograms.
 * 
 * @param signal Pointer to the signal histogram.
 * @param background Pointer to the background histogram.
 * @param data Pointer to the observed data histogram.
 * @return std::vector<double> Vector of optimized parameters and their errors.
 */
std::vector<double> RunQOptimization(TH1* signal, TH1* background, TH1* data);

/**
 * @brief Perform TFractionFitter fit.
 * 
 * This function performs a fit using the TFractionFitter class on signal, background, and negative data histograms.
 * 
 * @param hsig Pointer to the signal histogram.
 * @param hbkg Pointer to the background histogram.
 * @param hneg Pointer to the negative data histogram.
 * @param nsig Reference to the number of signal events.
 * @param nsig_error Reference to the error in the number of signal events.
 */
void TFracFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig, double &nsig_error);
