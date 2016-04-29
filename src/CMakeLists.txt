set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	AmplitudePair.cxx
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedDataValue.cxx
	CachedValue.cxx
	ClebschGordan.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPointBase.cxx
	DataSetBase.cxx
	DecayChannel.cxx
	DecayingParticle.cxx
	FinalStateParticle.cxx
	Flatte.cxx
	FourMomenta.cxx
	HelicityAngles.cxx
	HelicityFormalism.cxx
	MassShape.cxx
	MassShapeWithNominalMass.cxx
	MeasuredBreakupMomenta.cxx
	Model.cxx
	ParticleCombination.cxx
	ParticleCombinationCache.cxx
	Particle.cxx
	ParticleFactory.cxx
	PoleMass.cxx
	QuantumNumbers.cxx
	Resonance.cxx
	SpinAmplitude.cxx
	StatusManager.cxx
	WignerD.cxx
	ZemachFormalism.cxx
)

make_static_library(YAP "${YAP_SOURCES}" "${ROOT_LIBS}")
#make_shared_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
