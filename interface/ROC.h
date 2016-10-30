#ifndef ROC_H
#define ROC_H

#include <string> // std::string
#include <vector> // std::vector<>

namespace tthMEM
{
  /**
   * @brief Class which draws the ROC curve and calculates AUC (area under curve)
   *
   * @todo Extrapolate the class to a variable number of backgrounds
   */
  class ROC
  {
  public:
    ROC() = delete;
    ROC(const std::string & signalFileName,
        const std::vector<std::string> & bkgFileNames,
        const std::string & outFolderName,
        const std::string & treeName,
        const std::string & branchName,
        const std::vector<std::string> & labels);

    /**
     * @brief Sets legend position (in relative coordinates)
     * @param legPos (x1, y1, x2, y2) or lower left and upper right corner of the legend
     * @return EXIT_FAILURE (!= 0), if the provided coordinates are not valid,
     *         EXIT_SUCCESS (== 0) otherwise
     *
     * Note that the coordinates must be between 0 and 1, and x1 < x2, y1 < y2.
     */
    int
    setLegendPosition(const std::vector<double> & legPos);

    /**
     * @brief Saves the ROC curve to a file
     * @param idx Index of background sample
     * @return EXIT_FAILURE (!= 0), if there's enough input given to draw the cruve,
     *         EXIT_SUCCESS (== 0) otherwise
     */
    int
    plotROC();

    /**
     * @brief Calculates area under ROC curve (ROC AUC)
     * @param idx Index of background sample
     * @return ROC AUC (or -1, if not enough input given)
     *
     * Plots no title nor axis labels if the title_, xAxisTitle_ or yAxisTitle not set
     */
    double
    getAUC(unsigned idx);

    /**
     * @brief Calculates optimal cutoff point in the ROC curve
     * @param idx Index of background sample
     * @return The cutoff point (or -1, if not enough input given)
     *
     * There are various ways to calculate the optimal cutoff point, but the easiest
     * to understand is the method in which we find the point on the ROC curve closest
     * to the ideal working point (i.e. signal efficiency is 1 and background efficiency
     * is 0). Therefore we find the cutoff point (x, y) on ROC curve which minimizes
     * the distance function sqrt((x - 1)^2 + y^2).
     */
    double
    getOptimalCutoff(unsigned idx);

  private:
    std::string signalFileName_;            ///< signal file name
    std::vector<std::string> bkgFileNames_; ///< list of background file names
    std::string outFolderName_;             ///< output folder name
    std::string treeName_;                  ///< input tree name
    std::string branchName_;
    ///< branch name in which Neyman-Pearson scores are stored
    std::vector<std::string> labels_;       ///< labels describing each input file in one word

    std::vector<double> legPos_;           ///< legend position (in relative coordinates)
    const unsigned nofWPs = 101;           ///< number of working points (incl boundary)
    std::vector<double> signal;            ///< signal MEM scores
    std::vector<std::vector<double>> bkgs; ///< background MEM scores
    std::vector<double> wp;                ///< uniform vector of working points b/w 0 and 1

    std::vector<double> x; ///< ROC curve x-coordinates/signal efficiencies
    std::vector<std::vector<double>> ys;
    ///< ROC curve y-coordinates/background efficiencies

    /**
     * @brief Checks whether all input data is needed to find
     *        the coordinates of the ROC curve by checking:
     *   - signal file name
     *   - background file names
     *   - output folder name
     *   - input tree name
     *   - input branch name
     * @return false, if any of the above information is missing, true otherwise
     *
     * Note that if the output folder doesn't exist, it is created
     */
    bool
    checkIfReady() const;

    /**
     * @brief Reads input root files and stores MEM scores into a given vector
     * @param fileName Input root file name
     * @param values   Double vector of values/MEM scores (assumed empty)
     */
    void
    readInput(const std::string & fileName,
              std::vector<double> & values);

    /**
     * @brief Calculates the ROC curve coordinates by counting the number of events
     *        that pass any given working point, which is uniformly incremented
     *        between 0 and 1
     * @param values MEM scores
     * @param coords 1-dimensional coordinates corresponding to the MEM scores
     *
     * coords variable actually stores the efficiency for a range of working points
     */
    void
    getROCcoords(const std::vector<double> & values,
                 std::vector<double> & coords);
  };
}

#endif // ROC_H
