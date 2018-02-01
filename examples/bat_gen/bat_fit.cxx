#include "bat_fit.h"

#include "ConstantPrior.h"

#include <Attributes.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <SpinAmplitude.h>
#include <VariableStatus.h>

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCPrior.h>

#include <TTree.h>
#include <TRandom3.h>

#include <algorithm>
#include <assert.h>

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_yap_base(name, std::move(M)),
      FitData_(model()->createDataSet()),
      FitPartitions_(1, &FitData_),
      IntegralData_(model()->createDataSet()),
      IntegralPartitions_(1, &IntegralData_),
      FitFractionData_(model()->createDataSet()),
      FitFractionPartitions_(1, &FitFractionData_),
      NIntegrationPoints_(0),
      NIntegrationPointsBatchSize_(0),
      NIntegrationThreads_(1),
      Integral_(*model()),
      FitFractionIntegral_(*model()),
      FirstParameter_(-1),
      FirstObservable_(-1),
      useJacobian_(false),
      modelSelection_(0.)
{
    // create mass axes
    axes() = model()->massAxes(pcs);

    unsigned iVar(0);

    //
    // loop over all free amplitudes
    //
    // need to sort freeAmplitudes in a defined manner (default is by memory address)
    // sort alphabetically by name
    yap::FreeAmplitudeVector freeAmplitudes;
    for (auto& amp : free_amplitudes(*model())) {
        freeAmplitudes.push_back(amp);
    }
    std::sort(freeAmplitudes.begin(), freeAmplitudes.end(), [](const std::shared_ptr<yap::FreeAmplitude>& lhs, const std::shared_ptr<yap::FreeAmplitude>& rhs)
            {
                return (to_string(*lhs) < to_string(*rhs));
            });

    for (const auto& fa : freeAmplitudes) {
        // ignore fixed free amplitudes
        if (fa->variableStatus() == yap::VariableStatus::fixed)
            continue;

        // do not add shared free amplitudes
        bool cont(false);
        for (auto& famp : FreeAmplitudes_) {
            if (famp->freeAmplitude().get() == fa->freeAmplitude().get()) {
                LOG(INFO) << "free amplitude already in FreeAmplitudes_";
                cont = true;
                break;
            }
        }
        if (cont)
            continue;


        auto fa_name = std::to_string(iVar++) + " "
            + to_string(*fa->decayChannel())
            + " L = " + std::to_string(fa->spinAmplitude()->L())
            + " S = " + yap::spin_to_string(fa->spinAmplitude()->twoS());

        // add real parameter
        double range = 5. * abs(fa->value());
        AddParameter("real(" + fa_name + ")", -range, range);
        // add imag parameter
        AddParameter("imag(" + fa_name + ")", -range, range);

        AbsPriors_.push_back(std::make_unique<ConstantPrior>(0, range));
        ArgPriors_.push_back(std::make_unique<ConstantPrior>(-180, 180));

        // add amplitude observable
        AddObservable("amp(" + fa_name + ")", 0, range);
        // add phase observable
        AddObservable("arg(" + fa_name + ")", -180, 180);

        // add free amplitude to list
        FreeAmplitudes_.push_back(fa);
    }
    
    // loop over admixtures
    for (auto& comp : model()->components()) {
        // ignore fixed components
        if (comp.admixture()->variableStatus() == yap::VariableStatus::fixed)
            continue;
        auto adm_name = std::to_string(iVar++) + " "
                + "admixture_" + comp.particle()->name()
                + "_" + std::to_string(comp.decayTrees()[0]->initialTwoM());

        AddParameter(adm_name, 0, comp.admixture()->value() * 10.);

        Admixtures_.push_back(comp.admixture());
    }

    // // add observables for all fit fractions
    // int N = std::accumulate(Integral_.integrals().begin(), Integral_.integrals().end(), 0,
    //                         [](int n, const yap::IntegralMap::value_type& v)
    //                         {return n + v.second.decayTrees().size();});
    // if (N > 1) {
    /*for (const auto& mci : Integral_.integrals())
        for (const auto& dt : mci.Integral.decayTrees()) {
            DecayTrees_.push_back(dt);

            std::string str = std::to_string(iVar++) + " ""fit_frac(" + to_string(*dt) + ")";
            std::replace(str.begin(), str.end(), '-', 'm'); // - will be omitted by BATs "safe name", and when decay channels only differ in some spin projections (-1 vs 1), the safe name would be identical
            std::replace(str.begin(), str.end(), '\n', ';'); // make into one line
            std::replace(str.begin(), str.end(), '\t', ' ');

            AddObservable(str, 0, 1.1);
            //GetObservables().Back().SetNbins(1000);
        }*/
    // }

    FirstParameter_ = GetParameters().Size();

    assert(FirstParameter_ == int(2*FreeAmplitudes_.size() + Admixtures_.size()));
    FirstObservable_ = GetObservables().Size();
}

//-------------------------
std::vector<double> bat_fit::getInitialPositions(bool centerZeroImag) const
{
    std::vector<double> initialPositions;

    for (auto& fa : FreeAmplitudes_) {
        if (centerZeroImag and imag(fa->value()) == 0.) {
            initialPositions.push_back(0.);
            initialPositions.push_back(0.);
        }
        else {
            initialPositions.push_back(real(fa->value()));
            initialPositions.push_back(imag(fa->value()));
        }
    }

    for (auto& a : Admixtures_) {
        initialPositions.push_back(a->value());
    }

    for (auto& p : Parameters_) {
        // \todo might be nullptr
        initialPositions.push_back(dynamic_cast<yap::Parameter<double>*>(p.get())->value());
    }

    return initialPositions;
}

//-------------------------
std::vector<double> bat_fit::getRandomInitialPositions() const
{
    /*std::vector<double> initialPositions;

    for (auto& fa : FreeAmplitudes_) {
        while (true) {
            double range = std::max(1., 1.4 * abs(fa->value()));
            double real = gRandom->Uniform(-range, range);
            double imag = gRandom->Uniform(-range, range);

            if (real*real + imag*imag > range*range)
                continue;

            initialPositions.push_back(gRandom->Uniform(-range, range));
            initialPositions.push_back(gRandom->Uniform(-range, range));
            break;
        }

    }

    for (auto& a : Admixtures_) {
        initialPositions.push_back(2. * gRandom->Uniform(a->value()));
    }

    for (auto& p : Parameters_) {
        initialPositions.push_back(2. * gRandom->Uniform(dynamic_cast<yap::Parameter<double>*>(p.get())->value()));
    }

    if (!GetParameters().IsWithinLimits(initialPositions))
        return getRandomInitialPositions(); // try again

    return initialPositions;*/

    return GetParameters().GetUniformRandomValues(gRandom);
}

//-------------------------
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars)
{
    std::vector<std::vector<unsigned> > pcs;
    pcs.reserve(t_pars.GetEntries());

    bool parameter;
    char c_parname[10];
    t_pars.SetBranchAddress("parameter", &parameter);
    t_pars.SetBranchAddress("name", &c_parname);
    for (unsigned n = 0; n < t_pars.GetEntries(); ++n) {
        t_pars.GetEntry(n);
        if (!parameter)
            continue;
        std::string parname(c_parname);
        // parameter names should be m2_ij
        if (parname.find("m2_") != 0)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not match \"m2_ij\"",
                                             "find_mass_axes");
        if (parname.size() != 5)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" is not right length "
                                             + "(" + std::to_string(parname.size()) + " != 5)",
                                             "find_mass_axes");
        std::string indices_string = parname.substr(parname.rfind("_") + 1);
        if (indices_string.size() != 2)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not contain two one-digit indices", "find_mass_axes");
        if (!std::all_of(indices_string.begin(), indices_string.end(), isdigit))
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not contain two one-digit indices", "find_mass_axes");
        
        // read indices:
        std::vector<unsigned> indices;
        indices.reserve(2);
        std::transform(indices_string.begin(), indices_string.end(), std::back_inserter(indices), [](char c){return c - '0';});
        pcs.push_back(indices);
    }

    return pcs;
}

//-------------------------
void set_address(const yap::MassAxes::value_type& a,
                 std::vector<double>& m2,
                 TTree& t_mcmc)
{
    m2.push_back(0);
    t_mcmc.SetBranchAddress(indices_string(*a, "m2_", "").data(), &m2.back());
}

//-------------------------
size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A, double initial_mass, TTree& t_mcmc, int N, int lag)
{
    if (A.empty())
        throw yap::exceptions::Exception("mass axes empty", "load_data");

    LOG(INFO) << "Load with mass axes " << to_string(A);
    t_mcmc.Print();

    // set branch addresses
    std::vector<double> m2;
    m2.reserve(A.size());
    std::for_each(A.begin(), A.end(), std::bind(set_address, std::placeholders::_1, std::ref(m2), std::ref(t_mcmc)));
    if (m2.size() != A.size())
        throw yap::exceptions::Exception("not all mass axes loaded from TTree", "load_data");

    //
    // load data
    //
    int Phase = -1;
    t_mcmc.SetBranchAddress("Phase", &Phase);
    unsigned Iteration;
    t_mcmc.SetBranchAddress("Iteration", &Iteration);
    unsigned Chain;
    t_mcmc.SetBranchAddress("Chain", &Chain);

    unsigned long long n_entries = t_mcmc.GetEntries();


    if (N < 0)
        // attempt to load all data
        N = n_entries;

    if (lag < 0)
        // calculate lag
        lag = n_entries / N;
    lag = std::max(lag, 1);

    int n_attempted = 0;
    size_t old_size = data.size();

    for (unsigned long long n = 0; n < n_entries and n_attempted < N; ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        // if (fabs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
        //     continue;

        ++n_attempted;

        auto P = calculate_four_momenta(initial_mass, M, A, m2);
        if (P.empty())
            std::cout << "point is out of phase space!";
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)"
                << " from a tree of size " << n_entries << ", with a lag of " << lag;
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }

    if (int(data.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested. Reduce the lag (or set it to -1 to automatically determine the lag).";

    return data.size() - old_size;
}

//-------------------------
void bat_fit::setParameters(const std::vector<double>& p)
{
    DEBUG("p.size() = " << p.size());
    DEBUG("FreeAmplitudes_.size() = " << FreeAmplitudes_.size());
    DEBUG("Admixtures_.size()     = " << Admixtures_.size());
    DEBUG("Parameters_.size()     = " << Parameters_.size());
    DEBUG("FirstParameter_        = " << FirstParameter_);

    assert(p.size() >= 2 * FreeAmplitudes_.size() + Admixtures_.size() + Parameters_.size());

    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = std::complex<double>(p[i * 2], p[i * 2 + 1]);

    for (size_t i = 0; i < Admixtures_.size(); ++i)
        *Admixtures_[i] = p[FreeAmplitudes_.size()*2 + i];

    yap::set_values(Parameters_.begin(), Parameters_.end(),
            p.begin() + FirstParameter_, p.begin() + Parameters_.size());

    integrate();

    // calculate fit fractions
    // if (!CalculatedFitFractions_.empty()) {
    /*unsigned c = GetCurrentChain();
    size_t i = 0;
    for (const auto& mci : Integral_.integrals()) {
        auto ff = fit_fractions(mci.Integral);
        for (const auto& f : ff)
            CalculatedFitFractions_[c][i++] = f;
    }*/
    // }
}

//-------------------------
void bat_fit::integrate()
{

    if (IntegrationPointGenerator_)
        yap::ImportanceSampler::calculate(Integral_, IntegrationPointGenerator_, NIntegrationPoints_, NIntegrationPointsBatchSize_, NIntegrationThreads_);
    else {
        for (auto& p : IntegralPartitions_) {
            model()->updateCalculationStatus(*p);
        }

        if (yap::ImportanceSampler::calculate(Integral_, IntegralPartitions_) and modelSelection_ > 0.)
            yap::ImportanceSampler::calculate(FitFractionIntegral_, FitFractionPartitions_, true);
    }
}

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>& p)
{
    setParameters(p);

    auto integral_value = integral(Integral_).value();

    // print params
    std::cout << "\n";
    for (auto par : p)
        std::cout<<par<<"\t";
    std::cout << "\n";
    LOG(INFO) << "bat_fit::LogLikelihood; integral value = " << integral_value;

    double L = sum_of_log_intensity(*model(), FitPartitions_, log(integral_value));
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();

    // print L
    static unsigned nCalls(0);
    static double initL(L);
    static double maxL(L);
    maxL = std::max(maxL, L);
    LOG(INFO) << "\rLogLikelihood = " << L << "; \tmax = " << maxL
            << "; \tinitial = " << initL << "; \t calls = " << ++nCalls;

    return L;
}

//-------------------------
double bat_fit::LogAPrioriProbability(const std::vector<double>& p)
{
    //LOG(INFO) << "bat_fit::LogAPrioriProbability";
    double logP = 0;
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        if (GetParameter(i * 2).Fixed() or GetParameter(i * 2 + 1).Fixed())
            continue;
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        logP += AbsPriors_[i]->GetLogPrior(abs(A))
            + ArgPriors_[i]->GetLogPrior(yap::deg(arg(A)));
        if (useJacobian_)
            logP -= log(abs(A));      // jacobian
    }
    // parameters
    for (size_t i = FreeAmplitudes_.size() * 2; i < GetParameters().Size(); ++i)
        if (GetParameter(i).GetPrior())
            logP += GetParameter(i).GetLogPrior(p[i]);

    // model selection
    if (modelSelection_ > 0.) {

        // take normalized fit fractions of all components
        /*double sumIntegrals(0);
        // loop over admixtures
        for (const auto& mci : FitFractionIntegral_.integrals()) {
            sumIntegrals += mci.Admixture->value() * integral(mci.Integral).value();
        }

        double sum(0);
        for (const auto& mci : FitFractionIntegral_.integrals()) {
            //double sum_admixture(0);
            auto ff = fit_fractions(mci.Integral);
            for (size_t i = 0; i < ff.size(); ++i) {
                double fit_frac = mci.Admixture->value()  * integral(mci.Integral).value() / sumIntegrals * ff[i].value();
                //LOG(INFO) << "\t" << fit_frac*100. << " %; fit fraction in admixture = " << ff[i].value()*100. << " %";
                sum += fit_frac;
                //sum_admixture += ff[i].value();

                // LASSO
                logP -= fit_frac/modelSelection_;
            }
        }
        LOG(INFO) << "sum of fit fractions = " << 100.*sum << "%";*/

        // fit fraction priors
        // todo: do not hardcode
        static auto decayTrees = FitFractionIntegral_.integrals().at(0).Integral.decayTrees();
        static std::vector<yap::DecayTreeVector> groupedDecayTrees;
        {
            static BCGaussianPrior a1_rho_pi_S_prior(40., 5.);
            static BCGaussianPrior a1_sigma_pi_prior(9.6, 3.);
            static BCGaussianPrior rho_rho_prior(24., 5.);

            static std::vector<yap::DecayTreeVector> a1_rho_pi_S_groupedDecayTrees;
            static std::vector<yap::DecayTreeVector> a1_sigma_pi_groupedDecayTrees;
            static std::vector<yap::DecayTreeVector> rho_rho_groupedDecayTrees;
            static bool initialized(false);

            if (not initialized) {
                // a_1
                auto a_1_plus = particles(*model(), yap::is_named("a_1+")).empty() ? nullptr : std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), yap::is_named("a_1+")));
                auto rho = particles(*model(), yap::is_named("rho0")).empty() ? nullptr : std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), yap::is_named("rho0")));
                auto pipiS = particles(*model(), yap::is_named("pipiS")).empty() ? nullptr : std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), yap::is_named("pipiS")));

                if (a_1_plus) {

                    auto blub = yap::filter(decayTrees, yap::to(a_1_plus));
                    blub.erase(std::remove_if(blub.begin(), blub.end(),
                            [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), yap::to(rho), yap::l_equals(0))).empty();}),
                            blub.end());

                    LOG(INFO) << "a1_rho_pi_S_groupedDecayTrees";
                    a1_rho_pi_S_groupedDecayTrees.push_back(blub);
                    for (const auto& fa : a1_rho_pi_S_groupedDecayTrees[0])
                        LOG(INFO) << yap::to_string(*fa);

                    for(auto& dt : a1_rho_pi_S_groupedDecayTrees.back())  {
                        auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                        if(iter != decayTrees.end())
                            decayTrees.erase(iter);
                    }

                    blub = yap::filter(decayTrees, yap::to(a_1_plus));
                    blub.erase(std::remove_if(blub.begin(), blub.end(),
                            [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), yap::to(pipiS))).empty();}),
                            blub.end());

                    LOG(INFO) << "a1_sigma_pi_groupedDecayTrees";
                    a1_sigma_pi_groupedDecayTrees.push_back(blub);
                    for (const auto& fa : a1_sigma_pi_groupedDecayTrees[0])
                        LOG(INFO) << yap::to_string(*fa);

                    for(auto& dt : a1_sigma_pi_groupedDecayTrees.back())  {
                        auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                        if(iter != decayTrees.end())
                            decayTrees.erase(iter);
                    }
                }

                // rho rho
                if (rho) {
                    auto blub = yap::filter(decayTrees, yap::to(rho, rho));

                    LOG(INFO) << "rho_rho_groupedDecayTrees";
                    rho_rho_groupedDecayTrees.push_back(blub);
                    for (const auto& fa : rho_rho_groupedDecayTrees[0])
                        LOG(INFO) << yap::to_string(*fa);

                    for(auto& dt : rho_rho_groupedDecayTrees.back())  {
                        auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                        if(iter != decayTrees.end())
                            decayTrees.erase(iter);
                    }
                }

                LOG(INFO) << "------------------------------";
                groupedDecayTrees = group_by_free_amplitudes(decayTrees);

                initialized = true;
            }

            if (not a1_rho_pi_S_groupedDecayTrees.empty()) {
                double ff = 100. * fit_fractions(FitFractionIntegral_.integrals().at(0).Integral, a1_rho_pi_S_groupedDecayTrees).at(0).value();
                double prior = a1_rho_pi_S_prior.GetLogPrior(ff);
                LOG(INFO) << "a_1 -> (rho pi)_S fit fraction = " << ff << "%; " << prior;
                logP +=  prior;
            }

            if (not a1_sigma_pi_groupedDecayTrees.empty()) {
                double ff = 100. * fit_fractions(FitFractionIntegral_.integrals().at(0).Integral, a1_sigma_pi_groupedDecayTrees).at(0).value();
                double prior = a1_sigma_pi_prior.GetLogPrior(ff);
                LOG(INFO) << "a_1 -> (pipiS pi) fit fraction = " << ff << "%; " << prior;
                logP += prior;
            }

            if (not rho_rho_groupedDecayTrees.empty()) {
                double ff = 100. * fit_fractions(FitFractionIntegral_.integrals().at(0).Integral, rho_rho_groupedDecayTrees).at(0).value();
                double prior = rho_rho_prior.GetLogPrior(ff);
                LOG(INFO) << "D -> rho rho fit fraction = " << ff << "%; " << prior;
                logP += prior;
            }
        }

        auto logP_old = logP;

        // take fit fractions of signal component
        auto ff = fit_fractions(FitFractionIntegral_.integrals().at(0).Integral, groupedDecayTrees);
        double sum(0);
        for (size_t i = 0; i < ff.size(); ++i) {
            double fit_frac = ff[i].value();
            sum += fit_frac;

            // LASSO
            //logP -= fit_frac/modelSelection_;

            // BCM
            logP -= log(1. + fit_frac/(modelSelection_*modelSelection_));
        }
        LOG(INFO) << "sum of fit fractions = " << 100.*sum << "%; delta logP = " << logP - logP_old;

        //LOG(INFO) << to_string(mci.Integral);

    } // end if (modelSelection_ > 0.)

    LOG(INFO) << "bat_fit::LogAPrioriProbability = " << logP;

    return logP;
}

//-------------------------
void bat_fit::CalculateObservables(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        GetObservables()[i * 2 + 0] = abs(A);
        GetObservables()[i * 2 + 1] = yap::deg(arg(A));
    }
    /*unsigned c = GetCurrentChain();
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        GetObservables()[FreeAmplitudes_.size() * 2 + i] = CalculatedFitFractions_[c][i].value();
        */
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name + "_re", real(low), real(high), latex.empty() ? latex : "Re(" + latex + ")", units);
    AddParameter(name + "_im", imag(low), imag(high), latex.empty() ? latex : "Im(" + latex + ")", units);
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name, low, high, latex, units);
}

//-------------------------
size_t bat_fit::findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const
{
    auto it = std::find(FreeAmplitudes_.begin(), FreeAmplitudes_.end(), A);

    if (it == FreeAmplitudes_.end())
        throw yap::exceptions::Exception("FreeAmplitude not found", "setPrior");

    return (it - FreeAmplitudes_.begin()) * 2;
}

//-------------------------
void bat_fit::setPriors(std::shared_ptr<yap::FreeAmplitude> fa, BCPrior* amp_prior, BCPrior* arg_prior)
{
    //LOG(INFO) << "bat_fit::setPriors";
    if (!amp_prior)
        throw yap::exceptions::Exception("amp_prior is null", "bat_fit::setPrior");
    if (!arg_prior)
        throw yap::exceptions::Exception("phase_prior is null", "bat_fit::setPrior");
    
    auto i = findFreeAmplitude(fa);

    AbsPriors_[i / 2].reset(amp_prior);
    ArgPriors_[i / 2].reset(arg_prior);
}

//-------------------------
void bat_fit::setRealImagRanges(double rangeFactor)
{
    for (auto& fa : FreeAmplitudes_) {
        double range = rangeFactor * abs(fa->value());
        setRealImagRanges(fa, -range, range, -range, range);
    }
}

//-------------------------
void bat_fit::setAdmixtureRanges(double rangeFactor)
{
    unsigned i = FreeAmplitudes_.size() * 2;
    for (auto& adm : Admixtures_) {
        GetParameter(i++).SetLimits(0., rangeFactor * adm->value());
    }
}

//-------------------------
void bat_fit::setRealImagRanges(std::shared_ptr<yap::FreeAmplitude> fa, double real_low, double real_high, double imag_low, double imag_high)
{
    auto i = findFreeAmplitude(fa);
    GetParameter(i).SetLimits(real_low, real_high);
    GetParameter(i + 1).SetLimits(imag_low, imag_high);
}
    
//-------------------------
void bat_fit::setAbsArgRanges(std::shared_ptr<yap::FreeAmplitude> fa, double abs_low, double abs_high, double arg_low, double arg_high)
{
    auto i = findFreeAmplitude(fa);
    GetObservable(i).SetLimits(abs_low, abs_high);
    GetObservable(i + 1).SetLimits(arg_low, arg_high);
}

//-------------------------
void bat_fit::fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
{
    auto i = findFreeAmplitude(A);
    auto a = std::polar(amp, yap::rad(phase));
    GetParameter(i).Fix(real(a));
    GetParameter(i + 1).Fix(imag(a));
}

//-------------------------
void bat_fit::MCMCUserInitialize()
{
    bat_yap_base::MCMCUserInitialize();
    // if (!CalculatedFitFractions_.empty())
    //CalculatedFitFractions_.assign(GetNChains(), yap::RealIntegralElementVector(DecayTrees_.size()));
}

