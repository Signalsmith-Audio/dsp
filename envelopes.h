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
	
	/** Rectangular moving average filter */
	template<typename Sample=double>
	class BoxFilter {
		BoxSum<Sample> boxSum;
		int _size, _maxSize;
		Sample multiplier;
	public:
		BoxFilter(int maxSize) : boxSum(maxSize) {
			resize(maxSize);
		}
		/// Sets the maximum size (and current size, and resets)
		void resize(int maxSize) {
			_maxSize = maxSize;
			boxSum.resize(maxSize);
			set(maxSize);
		}
		/// Sets the current size (expanding/allocating only if needed)
		void set(int size) {
			_size = size;
			multiplier = Sample(1)/size;
			if (size > _maxSize) resize(size);
		}
		
		/// Resets (with an optional "fill" value)
		void reset(Sample fill=Sample()) {
			boxSum.reset(fill);
		}
		
		Sample operator()(Sample v) {
			return boxSum.readWrite(v, _size)*multiplier;
		}
	};

	/** FIR filter made from a stack of `BoxFilter`s.
		This filter has a non-negative impulse (monotonic step response), making it useful for smoothing positive-only values.  The internal box-lengths are chosen to minimise peaks in the stop-band.
			\diagram{box-stack-long.svg,Impulse responses for various stack sizes at length N=1000}
		Since the box-averages must have integer width, the frequency response is less accurate for shorter lengths with higher numbers of layers:
			\diagram{box-stack-short-freq.svg,Frequency responses for various stack sizes at length N=30}
	*/
	template<typename Sample=double>
	class BoxStackFilter {
		struct Layer {
			double ratio = 0, lengthError = 0;
			int length = 0;
			BoxFilter<Sample> filter{0};
			Layer() {}
		};
		int _size;
		std::vector<Layer> layers;
		void setupLayers(int layerCount) {
			layers.resize(layerCount, Layer{});
			double invN = 1.0/layerCount, sqrtN = std::sqrt(layerCount);
			double p = 1 - invN;
			double k = 1 + 4.5/sqrtN + 0.08*sqrtN;

			double sum = 0;
			for (int i = 0; i < layerCount; ++i) {
				double x = i*invN;
				double power = -x*(1 - p*std::exp(-x*k));
				double length = std::pow(2, power);
				layers[i].ratio = length;
				sum += length;
			}
			double factor = 1/sum;
			for (auto &layer : layers) layer.ratio *= factor;
		}
	public:
		BoxStackFilter(int maxSize, int layers=4) {
			resize(maxSize, layers);
		}
		
		/** Approximate bandwidth for a given number of layers
		\diagram{box-stack-bandwidth.svg,Approximate main lobe width (bandwidth)}
		*/
		static constexpr double layersToBandwidth(int layers) {
			return 1.58*(layers + 0.1);
		}
		/** Approximate peak in the stop-band
		\diagram{box-stack-peak.svg,Heuristic stop-band peak}
		*/
		static constexpr double layersToPeakDb(int layers) {
			return 5 - layers*18;
		}
		
		/// Sets the maximum (and current) impulse response length and layer count, and resets the filter.
		void resize(int maxSize, int layerCount) {
			if (int(layers.size()) != layerCount) setupLayers(layerCount);
			for (auto &layer : layers) {
				layer.filter.resize(int(maxSize*layer.ratio + 2));
			}
			set(maxSize);
			reset();
		}
		
		/// Sets the impulse response length (does not reset if `size` ≤ `maxSize`)
		void set(int size) {
			if (_size == size) return;
			_size = size;
			int order = size - 1;
			int totalOrder  = 0;
			
			for (auto &layer : layers) {
				double layerOrderFractional = layer.ratio*order;
				int layerOrder = int(layerOrderFractional);
				layer.length = layerOrder + 1;
				layer.lengthError = layerOrder - layerOrderFractional;
				totalOrder += layerOrder;
			}
			// Round some of them up, so the total is correct - this is O(N²), but `layers.size()` is small
			while (totalOrder < order) {
				int minIndex = 0;
				double minError = layers[0].lengthError;
				for (size_t i = 1; i < layers.size(); ++i) {
					if (layers[i].lengthError < minError) {
						minError = layers[i].lengthError;
						minIndex = i;
					}
				}
				layers[minIndex].length++;
				layers[minIndex].lengthError += 1;
				totalOrder++;
			}
			for (auto &layer : layers) layer.filter.set(layer.length);
		}

		/// Resets the filter
		void reset(Sample fill=Sample()) {
			for (auto &layer : layers) layer.filter.reset(fill);
		}
		
		Sample operator()(Sample v) {
			for (auto &layer : layers) {
				v = layer.filter(v);
			}
			return v;
		}
	};
	
	/** Peak-hold filter.
		\diagram{peak-hold.svg}
		
		The size is variable, and can be changed instantly with `.set()`, or by using `.push()`/`.pop()` in an unbalanced way.

		To avoid allocations while running, it uses a fixed-size array (not a `std::deque`) which determines the maximum length.
	*/
	template<typename Sample=double>
	class PeakHold {
		int frontIndex = 0, backIndex = 0;
		int indexMask;
		struct Segment {
			Sample value;
			int length;
		};
		std::vector<Segment> segments;
	public:
		PeakHold(int maxLength) {
			resize(maxLength);
		}
		/// Sets the maximum guaranteed size.
		void resize(int maxLength) {
			int length = 1;
			while (length < maxLength) length *= 2;
			indexMask = length - 1;
			
			frontIndex = backIndex = 0;
			// Sets the size
			segments.assign(length, Segment{0, maxLength});
		}
		/// Calculates the current size
		int size() const {
			int result = 0;
			for (int i = frontIndex; i <= backIndex; ++i) {
				result += segments[i&indexMask].length;
			}
			return result;
		}
		/// Resets the filter, preserving the size
		void reset(Sample fill=Sample()) {
			int s = size();
			frontIndex = backIndex = 0;
			segments[0] = {fill, s};
		}
		
		/// Sets a new size, extending the oldest value if needed
		void set(int newSize) {
			int oldSize = size();
			if (newSize > oldSize) {
				auto &front = segments[frontIndex&indexMask];
				front.length += newSize - oldSize;
			} else {
				for (int i = 0; i < oldSize - newSize; ++i) {
					pop();
				}
			}
		}
		/// Adds a new value and drops an old one.
		Sample operator()(Sample v) {
			pop();
			push(v);
			return read();
		}

		/// Drops the oldest value from the peak tracker.
		void pop() {
			auto &front = segments[frontIndex&indexMask];
			if (--front.length == 0) {
				++frontIndex;
			}
		}
		/// Adds a value to the peak tracker
		void push(Sample v) {
			int length = 1;
			while (backIndex >= frontIndex && segments[backIndex&indexMask].value <= v) {
				// Consume the segment
				length += segments[backIndex&indexMask].length;
				--backIndex;
			}
			++backIndex;
			segments[backIndex&indexMask] = {v, length};
		}
		/// Returns the current value
		Sample read() const {
			return segments[frontIndex&indexMask].value;
		}
	};
	
	/** Peak-decay filter.
		\diagram{peak-hold.svg}
		This always returns a value greater than the input sample, and (given a constant input) will reach that value in a fixed amount of time.
	*/
	template<typename Sample=double>
	class PeakDecayLinear {
		PeakHold<Sample> peakHold;
		Sample value = 0;
		Sample stepMultiplier = 1;
	public:
		PeakDecayLinear(int maxLength) : peakHold(maxLength) {
			set(maxLength);
		}
		void resize(int maxLength) {
			peakHold(maxLength);
			set(maxLength);
		}
		void set(double length) {
			peakHold.set(length);
			stepMultiplier = Sample(1)/length;
		}
		void reset(Sample start=Sample()) {
			peakHold.reset(start);
			value = start;
		}
		
		Sample operator ()(Sample v) {
			Sample peak = peakHold.read();
			peakHold(v);
			return value = std::max<Sample>(v, value + (v - peak)*stepMultiplier);
		}
	};

/** @} */
}} // signalsmith::envelopes::
#endif // include guard
