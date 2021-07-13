#ifndef SIGNALSMITH_DSP_ENVELOPES_H
#define SIGNALSMITH_DSP_ENVELOPES_H

#include "./common.h"

#include <cmath>
#include <random>

namespace signalsmith {
namespace envelopes {
	/**	@defgroup Envelopes Envelopes and LFOs
		@brief Envelopes, mostly fast and approximate
		
		@{
		@file
	*/
	
	/**	A cheap vaguely-sine-like LFO with random speed variation.
		
		Currently based on a cubic polynomial, but may change in future.
		
		It supports random variation in the rate, and the depth.  It will always remain bounded by the limits you specify (once it updates) so randomising the depth will produce smaller oscillations on average.
	*/
	class CheapLfo {
		float ratio = 0;
		float ratioStep = 0;
		
		float valueFrom = 0, valueTo = 1, valueRange = 1;
		float targetLow = 0, targetHigh = 1;
		float targetRate = 0;
		float rateRandom = 0.5, depthRandom = 0;
		bool freshReset = true;
		
		std::default_random_engine randomEngine;
		std::uniform_real_distribution<float> randomUnit;
		float random() {
			return randomUnit(randomEngine);
		}
		float randomRate() {
			return targetRate*exp(rateRandom*(random() - 0.5));
		}
		float randomTarget(float previous) {
			float randomOffset = depthRandom*random()*(targetLow - targetHigh);
			if (previous < (targetLow + targetHigh)*0.5f) {
				return targetHigh + randomOffset;
			} else {
				return targetLow - randomOffset;
			}
		}
	public:
		CheapLfo() : randomEngine(std::random_device()()), randomUnit(0, 1) {
			reset();
		}

		/// Resets the LFO state, with random phase.
		void reset() {
			ratio = random();
			ratioStep = randomRate();
			if (random() < 0.5) {
				valueFrom = targetLow;
				valueTo = targetHigh;
			} else {
				valueFrom = targetHigh;
				valueTo = targetLow;
			}
			valueRange = valueTo - valueFrom;
			freshReset = true;
		}
		/** Smoothly updates the LFO parameters.

		If called directly after `.reset()`, oscillation will immediately start within the specified range.  Otherwise, it will follow the new parameters after at most one cycle.

		The LFO will complete a full oscillation in (approximately) `1/rate` samples.  `rateVariation` can be any number, but 0-1 is a good range.  `depthVariation` must be in [0, 1], where 0.5 has random amplitude but still alternates up/down.
		*/
		void set(float low, float high, float rate, float rateVariation=0.5, float depthVariation=0) {
			rate *= 2; // We want to go up and down during this period
			targetRate = rate;
			targetLow = std::min(low, high);
			targetHigh = std::max(low, high);
			rateRandom = rateVariation;
			depthRandom = depthVariation;
			
			// If we haven't called .next() yet, don't bother being smooth.
			if (freshReset) return reset();

			// Only update the current rate if it's outside our new random-variation range
			float maxRandomRatio = exp((float)0.5*rateRandom);
			if (ratioStep > rate*maxRandomRatio || ratioStep < rate/maxRandomRatio) {
				ratioStep = randomRate();
			}
		}
		
		/// get the next output sample
		float next() {
			freshReset = false;
			float result = ratio*ratio*(3 - 2*ratio)*valueRange + valueFrom;

			ratio += ratioStep;
			if (ratio >= 1) {
				ratio = 0;
				ratioStep = randomRate();
				valueFrom = valueTo;
				valueTo = randomTarget(valueFrom);
				valueRange = valueTo - valueFrom;
			}
			return result;
		}
	};
	

/** @} */
}} // signalsmith::envelopes::
#endif // include guard
