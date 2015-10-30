#include "InitialStateParticle.h"

#include "Constants.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "logging.h"
#include "make_unique.h"

#include <TLorentzRotation.h>

#include <assert.h>
#include <future>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(q, mass, name, radialSize),
    Prepared_(false),
    FourMomenta_(),
    MeasuredBreakupMomenta_(),
    HelicityAngles_()
{
    FourMomenta_.setInitialStateParticle(this);
    MeasuredBreakupMomenta_.setInitialStateParticle(this);
    HelicityAngles_.setInitialStateParticle(this);

    setInitialStateParticle(this);
}

//-------------------------
InitialStateParticle::~InitialStateParticle()
{
    DataAccessors_.clear();
}

//-------------------------
std::complex<double> InitialStateParticle::amplitude(DataPoint& d, unsigned dataPartitionIndex) const
{
    std::complex<double> a = Complex_0;
    for (auto& pc : particleCombinations())
        a += amplitude(d, dataPartitionIndex);
    return a;
}

//-------------------------
double InitialStateParticle::partialSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D)
{
    if (!hasDataPartition(D))
        return 0;

    double L = 0;

    // loop over data points in partition
    for (DataIterator d = D->begin(); d != D->end(); ++d) {

        // reset calculation flags
        resetCalculationStatuses(D->index());

        L += logOfSquaredAmplitude(*d, D->index());
    }

    return L;
}

//-------------------------
double InitialStateParticle::sumOfLogsOfSquaredAmplitudes()
{
    // update global caclulationStatus's before looping over partitions
    updateGlobalCalculationStatuses();

    // create thread for calculation on each partition
    std::vector<std::future<double> > L;
    for (auto& D : DataPartitions_)
        L.push_back(std::async(&InitialStateParticle::partialSumOfLogsOfSquaredAmplitudes, this, D.get()));

    // wait for each partition to finish calculating
    double sum = 0;
    for (auto& l : L)
        sum += l.get();

    return sum;
}

//-------------------------
double InitialStateParticle::logLikelihood(DataPartitionBase* D)
{
    if (!hasDataPartition(D))
        return 0;

    /// \todo implement
    DEBUG("InitialStateParticle::logLikelihood()");

    // loop over DataPoints
    /// \todo investigate what needs to be done to use 'for (X : Y)'
    for (DataIterator d = D->begin(); d != D->end(); ++d) {
        DEBUG("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ");

        // reset calculation flags
        resetCalculationStatuses(D->index());

        // calculate amplitudes
        std::complex<double> a = Complex_0;
        for (auto& pc : particleCombinations())
            a += amplitude(*d, pc, D->index());

        DEBUG("InitialStateParticle amplitude = " << a);
    }


    setCachedDataValueFlagsToUnchanged(D->index());

    return 0;
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
bool InitialStateParticle::prepare()
{
    // check
    if (!consistent()) {
        LOG(ERROR) << "Cannot prepare InitialStateParticle, it is not consistent.";
        return false;
    }

    // particle combinations
    std::vector<std::shared_ptr<ParticleCombination> > PCs;
    for (auto& pc : particleCombinations()) {
        PCs.push_back(std::make_shared<ParticleCombination>(*pc));
    }
    ParticleCombination::makeParticleCombinationSetWithParents(PCs);
    clearSymmetrizationIndices();
    for (auto& pc : PCs)
        addSymmetrizationIndex(pc);

    // check
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (not pc->consistent()) {
            LOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombinationSet is not consistent.";
            return false;
        }
        if (pc->indices().size() < particleCombinations()[0]->indices().size() and not pc->parent()) {
            LOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombination is not consistent.";
            return false;
        }
    }

    //
    setSymmetrizationIndexParents();
    optimizeSpinAmplitudeSharing();

    // make sure that final state particles get the correct indices
    for (unsigned i = 0; i < particleCombinations()[0]->indices().size(); ++i) {
        std::shared_ptr<ParticleCombination> index = std::make_shared<ParticleCombination>(i);
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivByOrderlessContent(index, pc)) {
                FourMomenta_.addSymmetrizationIndex(pc);
                continue;
            }
        }
    }

    // add particle combinations to FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        FourMomenta_.addSymmetrizationIndex(pc);

        if (pc->indices().size() > 1) {
            HelicityAngles_.addSymmetrizationIndex(pc);
            MeasuredBreakupMomenta_.addSymmetrizationIndex(pc);
        }
    }

    FourMomenta_.findInitialStateParticle();

    setDataAcessorIndices();

    // fill DecayChannels_
    DecayChannels_.clear();
    for (DataAccessor* d : DataAccessors_) {
        DecayChannel* ch = dynamic_cast<DecayChannel*>(d);
        if (ch)
            DecayChannels_.push_back(ch);
    }

    // check
    if (!consistent()) {
        LOG(ERROR) << "Something went wrong while preparing InitialStateParticle, it is not consistent anymore.";
        return false;
    }

    Prepared_ = true;

    return true;
}

//-------------------------
std::vector<std::shared_ptr<ComplexParameter> > InitialStateParticle::freeAmplitudes() const
{
    std::vector<std::shared_ptr<ComplexParameter> > amps;
    amps.reserve(DecayChannels_.size());

    for (DecayChannel* ch : DecayChannels_)
        amps.push_back(ch->freeAmplitude());

    return amps;
}

//-------------------------
std::vector<DataPartitionBase*> InitialStateParticle::dataPartitions()
{
    std::vector<DataPartitionBase*> partitions;
    partitions.reserve(DataPartitions_.size());

    for (auto& d : DataPartitions_)
        partitions.push_back(d.get());

    return partitions;
}

//-------------------------
void InitialStateParticle::setDataPartitions(std::vector<DataPartitionBase*> partitions)
{
    DataPartitions_.clear();
    unsigned index(0);

    for (DataPartitionBase* p : partitions) {
        DataPartitions_.emplace_back(p->clone());
        DataPartitions_.back()->setIndex(index++);
    }

    setNumberOfDataPartitions(DataPartitions_.size());
}

//-------------------------
bool InitialStateParticle::addDataPoint(const std::vector<TLorentzVector>& fourMomenta)
{
    if (!Prepared_) {
        LOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        return false;
    }

    if (DataSet_.empty()) {
        return addDataPoint(DataPoint(fourMomenta));
    }

    DataSet_.push_back(DataPoint(DataSet_[0]));

    DataPoint& d = DataSet_.back();

    if (! d.setFourMomenta(fourMomenta))
        return false;

    FourMomenta_.calculate(d);
    MeasuredBreakupMomenta_.calculate(d);
    HelicityAngles_.calculate(d);

    if (!DataSet_.consistent(d))
        return false;

    return true;
}

//-------------------------
bool InitialStateParticle::addDataPoint(DataPoint&& d)
{
    if (!Prepared_) {
        LOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        return false;
    }

    d.allocateStorage(FourMomenta_, DataAccessors_);

    FourMomenta_.calculate(d);
    MeasuredBreakupMomenta_.calculate(d);
    HelicityAngles_.calculate(d);

    if (!DataSet_.consistent(d))
        return false;

    DataSet_.push_back(d);
    return true;
}

//-------------------------
bool InitialStateParticle::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}

//-------------------------
void InitialStateParticle::printDataAccessors(bool printParticleCombinations)
{
    // header
    std::cout << "DataAccessors of " << name() << "\n"
              << "index \tnSymIndices \taddress  \tname";
    if (printParticleCombinations)
        std::cout << "\t\tparticleCombinations";
    std::cout << "\n";

    for (DataAccessor* d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->maxSymmetrizationIndex() + 1
                  << "  \t\t" << d
                  << "  \t(" << typeid(*d).name() << ")  \t";
        if (dynamic_cast<Particle*>(d))
            std::cout << dynamic_cast<Particle*>(d)->name();
        else if (dynamic_cast<DecayChannel*>(d))
            std::cout << std::string(*dynamic_cast<DecayChannel*>(d));

        if (printParticleCombinations) {
            std::cout << " \t";

            for (auto& pc : d->particleCombinations())
                std::cout << std::string(*pc) << ";  ";
        }

        std::cout << "\n";
    }
    std::cout << std::endl;
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    unsigned size = particleCombinations()[0]->indices().size();
    assert(size > 1);

    clearSymmetrizationIndices();

    // get initial state PCs from set and add them
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->indices().size() < size)
            continue;

        if (pc->daughters()[0]->parent() == nullptr)
            continue;

        addSymmetrizationIndex(pc);
    }

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

}

//-------------------------
bool InitialStateParticle::hasDataPartition(DataPartitionBase* d)
{
    for (auto& dp : DataPartitions_) {
        if (d == dp.get())
            return true;
    }

    LOG(ERROR) << "InitialStateParticle::hasDataPartition - trying to calculate for a DataPartition, which is not stored in this InitialStateParticle.";
    return false;
}

//-------------------------
void InitialStateParticle::setNumberOfDataPartitions(unsigned n)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->setNumberOfDataPartitions(n);
        }
    }
}

//-------------------------
void InitialStateParticle::updateGlobalCalculationStatuses()
{
    for (DataAccessor* d : DataAccessors_) {
        /// \todo Think hard about a less stupid solution
        if (d == this)
            continue;
        d->updateGlobalCalculationStatuses();
    }
}

//-------------------------
void InitialStateParticle::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        //DEBUG("resetCalculationStatus for " << typeid(*d).name() << "  " << d);

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->resetCalculationStatus(dataPartitionIndex);
        }
    }
}

//-------------------------
void InitialStateParticle::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->setVariableStatus(kUnchanged, dataPartitionIndex);
        }
    }
}

//-------------------------
void InitialStateParticle::setParameterFlagsToUnchanged()
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            for (auto& p : c->ParametersItDependsOn_) {
                if (p->variableStatus() == kChanged) {
                    p->setVariableStatus(kUnchanged);
                }
            }
        }
    }
}

//-------------------------
void InitialStateParticle::addDataAccessor(DataAccessor* d)
{
    DataAccessors_.insert(d);
    if (prepared()) {
        LOG(ERROR) << "InitialStateParticle has already been prepared. "
                   << "Do NOT modify/add DecayChannels etc. after calling InitialStateParticle::prepare(), "
                   << "otherwise it will become inconsistent!";
    }
}


//-------------------------
void InitialStateParticle::removeDataAccessor(DataAccessor* d)
{
    if (! DataAccessors_.empty()) {

        DataAccessors_.erase(d);

        if (prepared()) {
            LOG(ERROR) << "InitialStateParticle has already been prepared. "
                       << "Do NOT modify/add DecayChannels etc. after calling InitialStateParticle::prepare(), "
                       << "otherwise it will become inconsistent!";
        }
    }
}

//-------------------------
void InitialStateParticle::setDataAcessorIndices()
{
    unsigned i(0);
    for (DataAccessor* d : DataAccessors_)
        d->setIndex(i++);
}


}
