#include "RooRealVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TApplication.h"
//ExpGausExp pdf in RooFit
class expgausexp : public RooAbsPdf {
    public:
        expgausexp() {}
        expgausexp(const char *name, const char *title,
                   RooAbsReal& _x,
                   RooAbsReal& _x0,
                   RooAbsReal& _sigmaL,
                   RooAbsReal& _alphaL,
                   RooAbsReal& _sigmaR,
                   RooAbsReal& _alphaR) :
            RooAbsPdf(name, title),
            x("x", "x", this, _x),
            x0("x0", "x0", this, _x0),
            sigmaL("sigmaL", "sigmaL", this, _sigmaL),
            alphaL("alphaL", "alphaL", this, _alphaL),
            sigmaR("sigmaR", "sigmaR", this, _sigmaR),
            alphaR("alphaR", "alphaR", this, _alphaR) {}

        expgausexp(const expgausexp& other, const char* name = 0) :
            RooAbsPdf(other, name),
            x("x", this, other.x),
            x0("x0", this, other.x0),
            sigmaL("sigmaL", this, other.sigmaL),
            alphaL("alphaL", this, other.alphaL),
            sigmaR("sigmaR", this, other.sigmaR),
            alphaR("alphaR", this, other.alphaR) {}

        virtual TObject* clone(const char* newname) const { return new expgausexp(*this, newname); }
        inline virtual ~expgausexp() {}

        // Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
        //     if (matchArgs(allVars, analVars, x)) return 1;
        //     return 0;
        // }

        // Double_t analyticalIntegral(Int_t code, const char* rangeName) const {
        //     assert(code == 1);
        //     constexpr double sqrtPiOver2 = 1.2533141373;
        //     constexpr double sqrt2 = 1.4142135624;
            double xmin = x.min(rangeName);
            double xmax = x.max(rangeName);
            double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
            double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
            double sum = 0;
            if (tmin < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double lv = tmin;
                double uv = std::min(tmax, -alphaL->getVal());
                sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
            }
            if (tmax > alphaR) {
                double a = 0.5 * alphaR * alphaR;
                double lv = std::max(tmin, alphaR->getVal());
                double uv = tmax;
                sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
            }
            if (tmin < alphaR && tmax > -alphaL) {
                double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
                double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
                sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR->getVal()) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
            }
            return sum;
        // }

    protected:
        RooRealProxy x;
        RooRealProxy x0;
        RooRealProxy sigmaL;
        RooRealProxy alphaL;
        RooRealProxy sigmaR;
        RooRealProxy alphaR;

        Double_t evaluate() const {
            constexpr double sqrtPiOver2 = 1.2533141373;
            constexpr double sqrt2 = 1.4142135624;

            double t = (x - x0) / (x < x0 ? sigmaL : sigmaR);
            double v = 0;
            if (t < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double b = alphaL * t;
                v = std::exp(a + b);
            } else if (t <= alphaR) {
                v = std::exp(-0.5 * t * t);
            } else {
                double a = 0.5 * alphaR * alphaR;
                double b = alphaR * (-t);
                v = std::exp(a + b);
            }

            return v;
        }
    };


void integral_test() {
    // 创建变量和参数
    RooRealVar x("x", "x", -10, 10);
    RooRealVar x0("x0", "x0", 0.0);
    RooRealVar sigmaL("sigmaL", "sigmaL", 1.0);
    RooRealVar alphaL("alphaL", "alphaL", 1.0);
    RooRealVar sigmaR("sigmaR", "sigmaR", 1.0);
    RooRealVar alphaR("alphaR", "alphaR", 1.0);
    
    // 创建 PDF 实例
    expgausexp pdf("pdf", "pdf", x, x0, sigmaL, alphaL, sigmaR, alphaR);
    
 
    
    // 2. 测试部分范围积分
    {
        // 设置部分积分范围
        x.setRange("partial", -2, 3);
        
        // 使用数值积分计算部分范围积分
        //pdf.declareAnalyticalIntegral(0);
        std::unique_ptr<RooAbsReal> numIntPartial(pdf.createIntegral(x, "partial"));
        double numIntPartialVal = numIntPartial->getVal();
        
        // 使用解析积分计算部分范围积分
       // double anaIntPartialVal = pdf.analyticalIntegral(1, "partial");
        
        // 打印结果
        std::cout << "\n=== Partial Range Test [-2, 3] ===" << std::endl;
        std::cout << "Numerical integral: " << numIntPartialVal << std::endl;
        //std::cout << "Analytical integral: " << anaIntPartialVal << std::endl;
        //std::cout << "Difference: " << std::abs(numIntPartialVal - anaIntPartialVal) << std::endl;
    }
    

}

