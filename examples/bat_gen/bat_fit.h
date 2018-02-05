#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "bat_yap_base.h"

#include <fwd/DataPartition.h>
#include <fwd/FreeAmplitude.h>
#include <fwd/IntegralElement.h>
#include <fwd/Model.h>
#include <fwd/Parameter.h>

#include <DataSet.h>
#include <FourVector.h>
#include <ModelIntegral.h>

#include <memory>
#include <functional>
#include <string>
#include <random>
#include <vector>

class TTree;

class bat_fit : public bat_yap_base
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {});

    /// add a complex parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex = "", std::string units = "");

    /// ada a real parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex = "", std::string units = "");

    /// set the priors for a FreeAmplitude's amplitude and phase
    void setPriors(std::shared_ptr<yap::FreeAmplitude> fa, BCPrior* amp_prior, BCPrior* arg_prior);
    
    /// set the range for a FreeAmplitude's real and imaginary parts
    void setRealImagRanges(std::shared_ptr<yap::FreeAmplitude> fa, double real_low, double real_high, double imag_low, double imag_high);
    
    void setRealImagRanges(double rangeFactor = 2.);

    void setAdmixtureRanges(double rangeFactor = 2.);

    /// set the range for a FreeAmplitude's abs and arg obvservables
    void setAbsArgRanges(std::shared_ptr<yap::FreeAmplitude> fa, double abs_low, double abs_high, double arg_low, double arg_high);

    void setUseJacobian(bool j = true)
    { useJacobian_ = j; }

    void setModelSelection(double s)
    { modelSelection_ = s; }

    /// fix a FreeAmplitude
    void fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase);

    // fix all free amplitudes
    void fixAmplitudes();

    // fix all free amplitudes
    void fixAdmixtures();

    /// log likelihood
    double LogLikelihood(const std::vector<double>& p) override;

    /// log prior
    double LogAPrioriProbability(const std::vector<double>& p) override;

    /// calculate  observables
    void CalculateObservables(const std::vector<double>& p) override;

    /// \return FitData_
    yap::DataSet& fitData()
    { return FitData_; }

    /// \return FitPartitions_
    yap::DataPartitionVector& fitPartitions()
    { return FitPartitions_; }

    /// \return IntegralData_
    yap::DataSet& integralData()
    { return IntegralData_; }

    /// \return IntegralPartitions_
    yap::DataPartitionVector& integralPartitions()
    { return IntegralPartitions_; }

    /// \return FitFractionData_
    yap::DataSet& fitFractionData()
    { return FitFractionData_; }

    /// \return FitFractionPartitions_
    yap::DataPartitionVector& fitFractionPartitions()
    { return FitFractionPartitions_; }

    /// \return FitFractionIntegral_
    yap::ModelIntegral& fitFractionIntegral()
    { return FitFractionIntegral_; }

    yap::ModelIntegral& modelIntegral()
    { return Integral_; }

    /// set parameters of integration
    /// \param N number of points to use for integration
    /// \param n batch size for integration
    /// \param t number of threads
    void setNIntegrationPoints(unsigned N, unsigned n, unsigned t = 1)
    { NIntegrationPoints_ = N; NIntegrationPointsBatchSize_ = n; NIntegrationThreads_ = t; }

    std::vector<double> getInitialPositions(bool centerZeroImag = false) const;
    std::vector<double> getRandomInitialPositions() const;

    /// \typedef Generator
    /// function for generating new points for integration
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;

    /// \return IntegrationPointGenerator_
    Generator& integrationPointGenerator()
    { return IntegrationPointGenerator_; }

    /// init size of CalculatedFitFractions_
    void MCMCUserInitialize() override;

    /// set parameters into model
    void setParameters(const std::vector<double>& p);

    /// perform the integration
    void integrate();

    /// find the position in the parameter list of the first element of a free amplitude
    size_t findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const;

    /// \return free amplitudes
    const yap::FreeAmplitudeVector& freeAmplitudes() const
    { return FreeAmplitudes_; }

    int firstParameter() const
    { return FirstParameter_; }

protected:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitPartitions_;

    /// DataSet to fit the model to
    yap::DataSet IntegralData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector IntegralPartitions_;

    /// DataSet to calculate fit fraction (phsp distributed)
    yap::DataSet FitFractionData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitFractionPartitions_;

    /// Number of points to integrate with
    unsigned NIntegrationPoints_;

    /// Batch size for generating integration points
    unsigned NIntegrationPointsBatchSize_;

    /// Number of threads for integration
    unsigned NIntegrationThreads_;

    /// generator for integration
    Generator IntegrationPointGenerator_;

    /// stores integral result
    yap::ModelIntegral Integral_;

    /// stores integral result
    yap::ModelIntegral FitFractionIntegral_;

    /// vector of parameters to set in model
    yap::ParameterVector Parameters_;

    /// offset of where first user-set parameter is
    int FirstParameter_;

    /// offset of where first user-set observable is
    int FirstObservable_;

    /// list of decay trees integrated over
    yap::DecayTreeVector DecayTrees_;

    /// Free amplitudes of model to set
    yap::FreeAmplitudeVector FreeAmplitudes_;

    yap::NonnegativeRealParameterVector Admixtures_;

    /// BCPrior on abs(FreeAmplitude)
    std::vector<std::unique_ptr<BCPrior> > AbsPriors_;

    /// BCPrior on arg(FreeAmplitude)
    std::vector<std::unique_ptr<BCPrior> > ArgPriors_;

    /// Calculated fit fractions (for observables)
    std::vector<yap::RealIntegralElementVector> CalculatedFitFractions_;

    // set to true when sampling
    bool useJacobian_;

    // use prior on fit fractions
    double modelSelection_;

};

/// load data from a TTree into a DataSet
/// \param data DataSet to load into
/// \param M Model to load with
/// \param A MassAxes to load with
/// \param t_mcmc TTree to load from
/// \param N max number of data points to (attempt to) load
/// \param lag Lag to apply to iterations when reading from TTree
///        per default, data points will be selected uniformly from t_mcmc
/// \param eps Amount to smear momenta by
size_t load_data(yap::DataSet& data, const yap::Model& M,
                 const yap::MassAxes& A, double initial_mass, TTree& t_mcmc,
                 int N = -1, int lag = -1);

/// find mass axes from TTree of parameters
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars);
#endif
