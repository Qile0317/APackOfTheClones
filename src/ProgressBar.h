#include <Rcpp.h>
#include <string>

class ProgressBar {
public:
    static void show(int x, int max) {
        double percent = 100.0 * double(x) / double(max);
        int filledWidth = static_cast<int>(percent * 0.5);
        int remainingWidth = 50 - filledWidth;

        std::string progressBar = "\r[";
        progressBar += std::string(filledWidth, '=');
        progressBar += std::string(remainingWidth, ' ') + "] ";
        progressBar += std::to_string(int(floor(percent))) + "%";

        Rcpp::Rcout << progressBar;
    }

    static void start() {
        show(0, 1);
    }

    static void finish() {
        show(1, 1);
    }
};
