#ifndef SIGNALSMITH_DSP_ENVELOPES_H
#define SIGNALSMITH_DSP_ENVELOPES_H

#include "./common.h"

#include <cmath>
#include <random>
#include <vector>

namespace signalsmith {
namespace envelopes {
	/**	@defgroup Envelopes Envelopes and LFOs
		@brief Envelopes, mostly fast and approximate
		
		@{
		@file
	*/
	
	/**	An LFO based on cubic segments.
		You can randomise the rate and/or the depth.  Randomising the depth past `0.5` means it no longer neatly alternates sides:
			\diagram{cubic-lfo-example.svg,Some example LFO curves.}
		Without randomisation, it is approximately sine-like:
			\diagram{cubic-lfo-spectrum-pure.svg}
	*/
	class CubicLfo {
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
		CubicLfo() : randomEngine(std::random_device()()), randomUnit(0, 1) {
			reset();
		}
		CubicLfo(long seed) : randomUnit(0, 1) {
			randomEngine.seed(seed);
			reset();
		}

		/// Resets the LFO state, starting with random phase.
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

		If called directly after `.reset()`, oscillation will immediately start within the specified range.  Otherwise, it will remain smooth and fit within the new range after at most one cycle:
			\diagram{cubic-lfo-changes.svg}

		The LFO will complete a full oscillation in (approximately) `1/rate` samples.  `rateVariation` can be any number, but 0-1 is a good range.
		
		`depthVariation` must be in the range [0, 1], where ≤ 0.5 produces random amplitude but still alternates up/down.
			\diagram{cubic-lfo-spectrum.svg,Spectra for the two types of randomisation - note the jump as depth variation goes past 50%}
		*/
		void set(float low, float high, float rate, float rateVariation=0, float depthVariation=0) {
			rate *= 2; // We want to go up and down during this period
			targetRate = rate;
			targetLow = std::min(low, high);
			targetHigh = std::max(low, high);
			rateRandom = rateVariation;
			depthRandom = std::min<float>(1, std::max<float>(0, depthVariation));
			
			// If we haven't called .next() yet, don't bother being smooth.
			if (freshReset) return reset();

			// Only update the current rate if it's outside our new random-variation range
			float maxRandomRatio = exp((float)0.5*rateRandom);
			if (ratioStep > rate*maxRandomRatio || ratioStep < rate/maxRandomRatio) {
				ratioStep = randomRate();
			}
		}
		
		/// Returns the next output sample
		float next() {
			freshReset = false;
			float result = ratio*ratio*(3 - 2*ratio)*valueRange + valueFrom;

			ratio += ratioStep;
			while (ratio >= 1) {
				ratio -= 1;
				ratioStep = randomRate();
				valueFrom = valueTo;
				valueTo = randomTarget(valueFrom);
				valueRange = valueTo - valueFrom;
			}
			return result;
		}
	};
	
	/** Variable-width rectangular sum */
	template<typename Sample=double>
	class BoxSum {
		int bufferLength, index;
		std::vector<Sample> buffer;
		Sample sum = 0, wrapJump = 0;
	public:
		BoxSum(int maxLength) {
			resize(maxLength);
		}

		/// Sets the maximum size (and reset contents)
		void resize(int maxLength) {
			bufferLength = maxLength + 1;
			buffer.resize(bufferLength);
			buffer.shrink_to_fit();
			reset();
		}
		
		/// Resets (with an optional "fill" value)
		void reset(Sample value=Sample()) {
			index = 0;
			sum = 0;
			for (size_t i = 0; i < buffer.size(); ++i) {
				buffer[i] = sum;
				sum += value;
			}
			wrapJump = sum;
			sum = 0;
		}
		
		Sample read(int width) {
			int readIndex = index - width;
			double result = sum;
			if (readIndex < 0) {
				result += wrapJump;
				readIndex += bufferLength;
			}
			return result - buffer[readIndex];
		}
		
		void write(Sample value) {
			++index;
			if (index == bufferLength) {
				index = 0;
				wrapJump = sum;
				sum = 0;
			}
			sum += value;
			buffer[index] = sum;
		}
		
		Sample readWrite(Sample value, int width) {
			write(value);
			return read(width);
		}
	};
	
	/** Variable-width rectangular moving average */
	template<typename Sample=double>
	class BoxAverage : private BoxSum<Sample> {
		using Super = BoxSum<Sample>;
		int _size, _maxSize;
		Sample multiplier;
	public:
		BoxAverage(int maxSize) : Super(maxSize) {
			resize(maxSize);
			set(maxSize);
		}
		/// Sets the maximum size (and reset contents)
		void resize(int maxSize) {
			_maxSize = maxSize;
			Super::resize(maxSize);
			if (maxSize < _size) set(maxSize);
		}
		/// Sets the current size (expanding size if needed)
		void set(int size) {
			_size = size;
			multiplier = Sample(1)/size;
			if (size > _maxSize) resize(size);
		}
		
		/// Resets (with an optional "fill" value)
		using Super::reset;
		using Super::write;
		Sample read() {
			return Super::read(_size)*multiplier;
		}
		Sample readWrite(Sample v) {
			return Super::readWrite(v, _size)*multiplier;
		}
	};

/** @} */
}} // signalsmith::envelopes::
#endif // include guard
