#ifndef SIGNALSMITH_DSP_ENVELOPES_H
#define SIGNALSMITH_DSP_ENVELOPES_H

#include "./common.h"

#include <cmath>
#include <random>
#include <vector>
#include <iterator>
#include <algorithm> // std::stable_sort

namespace signalsmith {
namespace envelopes {
	/**	@defgroup Envelopes Envelopes and LFOs
		@brief LFOs, envelopes and filters for manipulating them
		
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
			if (maxLength != 0) buffer.shrink_to_fit();
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
	
	/** Rectangular moving average filter (FIR).
		\diagram{box-filter-example.svg}
		A filter of length 1 has order 0 (i.e. does nothing). */
	template<typename Sample=double>
	class BoxFilter {
		BoxSum<Sample> boxSum;
		int _length, _maxLength;
		Sample multiplier;
	public:
		BoxFilter(int maxLength) : boxSum(maxLength) {
			resize(maxLength);
		}
		/// Sets the maximum size (and current size, and resets)
		void resize(int maxLength) {
			_maxLength = maxLength;
			boxSum.resize(maxLength);
			set(maxLength);
		}
		/// Sets the current size (expanding/allocating only if needed)
		void set(int length) {
			_length = length;
			multiplier = Sample(1)/length;
			if (length > _maxLength) resize(length);
		}
		
		/// Resets (with an optional "fill" value)
		void reset(Sample fill=Sample()) {
			boxSum.reset(fill);
		}
		
		Sample operator()(Sample v) {
			return boxSum.readWrite(v, _length)*multiplier;
		}
	};

	/** FIR filter made from a stack of `BoxFilter`s.
		This filter has a non-negative impulse (monotonic step response), making it useful for smoothing positive-only values.  It provides an optimal set of box-lengths, chosen to minimise peaks in the stop-band:
			\diagram{box-stack-long.svg,Impulse responses for various stack sizes at length N=1000}
		Since the underlying box-averages must have integer width, the peaks are slightly higher for shorter lengths with higher numbers of layers:
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

		template<class Iterable>
		void setupLayers(const Iterable &ratios) {
			layers.resize(0);
			double sum = 0;
			for (auto ratio : ratios) {
				Layer layer;
				layer.ratio = ratio;
				layers.push_back(layer);
				sum += ratio;
			}
			double factor = 1/sum;
			for (auto &l : layers) {
				l.ratio *= factor;
			}
		}
	public:
		BoxStackFilter(int maxSize, int layers=4) {
			resize(maxSize, layers);
		}

		/// Returns an optimal set of length ratios (heuristic for larger depths)
		static std::vector<double> optimalRatios(int layerCount) {
			// Coefficients up to 6, found through numerical search
			static double hardcoded[] = {1, 0.58224186169, 0.41775813831, 0.404078562416, 0.334851475794, 0.261069961789, 0.307944914938, 0.27369945234, 0.22913263601, 0.189222996712, 0.248329349789, 0.229253789144, 0.201191468123, 0.173033035122, 0.148192357821, 0.205275202874, 0.198413552119, 0.178256637764, 0.157821404506, 0.138663023387, 0.121570179349 /*, 0.178479592135, 0.171760666359, 0.158434068954, 0.143107825806, 0.125907148711, 0.11853946895, 0.103771229086, 0.155427880834, 0.153063152848, 0.142803459422, 0.131358358458, 0.104157805178, 0.119338029601, 0.0901675284678, 0.103683785192, 0.143949349747, 0.139813248378, 0.132051305252, 0.122216776152, 0.112888320989, 0.102534988632, 0.0928386714364, 0.0719750997699, 0.0817322396428, 0.130587011572, 0.127244563184, 0.121228748787, 0.113509941974, 0.105000272288, 0.0961938290157, 0.0880639725438, 0.0738389766046, 0.0746781936619, 0.0696544903682 */};
			if (layerCount <= 0) {
				return {};
			} else if (layerCount <= 6) {
				double *start = &hardcoded[layerCount*(layerCount - 1)/2];
				return std::vector<double>(start, start + layerCount);
			}
			std::vector<double> result(layerCount);

			double invN = 1.0/layerCount, sqrtN = std::sqrt(layerCount);
			double p = 1 - invN;
			double k = 1 + 4.5/sqrtN + 0.08*sqrtN;

			double sum = 0;
			for (int i = 0; i < layerCount; ++i) {
				double x = i*invN;
				double power = -x*(1 - p*std::exp(-x*k));
				double length = std::pow(2, power);
				result[i] = length;
				sum += length;
			}
			double factor = 1/sum;
			for (auto &r : result) r *= factor;
			return result;
		}
		/** Approximate (optimal) bandwidth for a given number of layers
		\diagram{box-stack-bandwidth.svg,Approximate main lobe width (bandwidth)}
		*/
		static constexpr double layersToBandwidth(int layers) {
			return 1.58*(layers + 0.1);
		}
		/** Approximate (optimal) peak in the stop-band
		\diagram{box-stack-peak.svg,Heuristic stop-band peak}
		*/
		static constexpr double layersToPeakDb(int layers) {
			return 5 - layers*18;
		}
		
		/// Sets size using an optimal (heuristic at larger sizes) set of length ratios
		void resize(int maxSize, int layerCount) {
			resize(maxSize, optimalRatios(layerCount));
		}
		/// Sets the maximum (and current) impulse response length and explicit length ratios
		template<class List>
		auto resize(int maxSize, List ratios) -> decltype(void(std::begin(ratios)), void(std::end(ratios))) {
			setupLayers(ratios);
			for (auto &layer : layers) layer.filter.resize(0); // .set() will expand it later
			_size = -1;
			set(maxSize);
			reset();
		}
		void resize(int maxSize, std::initializer_list<double> ratios) {
			resize<const std::initializer_list<double> &>(maxSize, ratios);
		}
		
		/// Sets the impulse response length (does not reset if `size` ≤ `maxSize`)
		void set(int size) {
			if (layers.size() == 0) return; // meaningless

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

		This has complexity O(1) every sample when the length remains constant (balanced `.push()`/`.pop()`, or using `filter(v)`), and amortised O(1) complexity otherwise.  To avoid allocations while running, it pre-allocates a vector (not a `std::deque`) which determines the maximum length.
	*/
	template<typename Sample>
	class PeakHold {
		static constexpr Sample lowest = std::numeric_limits<Sample>::lowest();
		int bufferMask;
		std::vector<Sample> buffer;
		int backIndex = 0, middleStart = 0, workingIndex = 0, middleEnd = 0, frontIndex = 0;
		Sample frontMax = lowest, workingMax = lowest, middleMax = lowest;
		
	public:
		PeakHold(int maxLength) {
			resize(maxLength);
		}
		int size() {
			return frontIndex - backIndex;
		}
		void resize(int maxLength) {
			int bufferLength = 1;
			while (bufferLength < maxLength) bufferLength *= 2;
			buffer.resize(bufferLength);
			bufferMask = bufferLength - 1;
			
			frontIndex = backIndex + maxLength;
			reset();
		}
		void reset(Sample fill=lowest) {
			int prevSize = size();
			buffer.assign(buffer.size(), fill);
			frontMax = workingMax = middleMax = lowest;
			middleEnd = workingIndex = frontIndex = 0;
			middleStart = middleEnd - (prevSize/2);
			backIndex = frontIndex - prevSize;
		}
		/** Sets the size immediately.
		Must be `0 <= newSize <= maxLength` (see constructor and `.resize()`).
		
		Shrinking doesn't destroy information, and if you expand again (with `preserveCurrentPeak=false`), you will get the same output as before shrinking.  Expanding when `preserveCurrentPeak` is enabled is destructive, re-writing its history such that the current output value is unchanged.*/
		void set(int newSize, bool preserveCurrentPeak=false) {
			while (size() < newSize) {
				Sample &backPrev = buffer[backIndex&bufferMask];
				--backIndex;
				Sample &back = buffer[backIndex&bufferMask];
				back = preserveCurrentPeak ? backPrev : std::max(back, backPrev);
			}
			while (size() > newSize) {
				pop();
			}
		}
		
		void push(Sample v) {
			buffer[frontIndex&bufferMask] = v;
			++frontIndex;
			frontMax = std::max(frontMax, v);
		}
		void pop() {
			if (backIndex == middleStart) {
				// Move along the maximums
				workingMax = lowest;
				middleMax = frontMax;
				frontMax = lowest;

				int prevFrontLength = frontIndex - middleEnd;
				int prevMiddleLength = middleEnd - middleStart;
				if (prevFrontLength <= prevMiddleLength + 1) {
					// Swap over simply
					middleStart = middleEnd;
					middleEnd = frontIndex;
					workingIndex = middleEnd;
				} else {
					// The front is longer than the middle - only happens if unbalanced
					// We don't move *all* of the front over, keeping half the surplus in the front
					int middleLength = (frontIndex - middleStart)/2;
					middleStart = middleEnd;
					middleEnd += middleLength;

					// Working index is close enough that it will be finished by the time the back is empty
					int backLength = middleStart - backIndex;
					int workingLength = std::min(backLength, middleEnd - middleStart);
					workingIndex = middleStart + workingLength;

					// Since the front was not completely consumed, we re-calculate the front's maximum
					for (int i = middleEnd; i != frontIndex; ++i) {
						frontMax = std::max(frontMax, buffer[i&bufferMask]);
					}
					// The index might not start at the end of the working block - compute the last bit immediately
					for (int i = middleEnd - 1; i != workingIndex - 1; --i) {
						buffer[i&bufferMask] = workingMax = std::max(workingMax, buffer[i&bufferMask]);
					}
				}

				// Is the new back (previous middle) empty? Only happens if unbalanced
				if (backIndex == middleStart) {
					 // swap over again (front's empty, no change)
					workingMax = lowest;
					middleMax = frontMax;
					frontMax = lowest;
					middleStart = workingIndex = middleEnd;
	
					if (backIndex == middleStart) {
						--backIndex; // Only happens if you pop from an empty list - fail nicely
					}
				}
				
				buffer[frontIndex&bufferMask] = lowest; // In case of length 0, when everything points at this value
			}

			++backIndex;
			if (workingIndex != middleStart) {
				--workingIndex;
				buffer[workingIndex&bufferMask] = workingMax = std::max(workingMax, buffer[workingIndex&bufferMask]);
			}
		}
		Sample read() {
			Sample backMax = buffer[backIndex&bufferMask];
			return std::max(backMax, std::max(middleMax, frontMax));
		}
		
		// For simple use as a constant-length filter
		Sample operator ()(Sample v) {
			push(v);
			pop();
			return read();
		}
	};
	
	/** Peak-decay filter with a linear shape and fixed-time return to constant value.
		\diagram{peak-decay-linear.svg}
		This is equivalent to a `BoxFilter` which resets itself whenever the output would be less than the input.
	*/
	template<typename Sample=double>
	class PeakDecayLinear {
		static constexpr Sample lowest = std::numeric_limits<Sample>::lowest();
		PeakHold<Sample> peakHold;
		Sample value = lowest;
		Sample stepMultiplier = 1;
	public:
		PeakDecayLinear(int maxLength) : peakHold(maxLength) {
			set(maxLength);
		}
		void resize(int maxLength) {
			peakHold.resize(maxLength);
			reset();
		}
		void set(double length) {
			peakHold.set(std::ceil(length));
			// Overshoot slightly but don't exceed 1
			stepMultiplier = Sample(1.0001)/std::max(1.0001, length);
		}
		void reset(Sample start=lowest) {
			peakHold.reset(start);
			set(peakHold.size());
			value = start;
		}
		
		Sample operator ()(Sample v) {
			Sample peak = peakHold.read();
			peakHold(v);
			return value = std::max<Sample>(v, value + (v - peak)*stepMultiplier);
		}
	};

	/// A real-valued cubic curve
	template<typename Sample=double>
	class CubicSegment {
		Sample xStart, a0, a1, a2, a3;
		
		// Only use with y0 != y1
		static inline double gradient(Sample x0, Sample x1, Sample y0, Sample y1) {
			return (y1 - y0)/(x1 - x0);
		}
		// Ensure a gradient produces monotonic segments
		static inline void ensureMonotonic(Sample &curveGrad, Sample gradA, Sample gradB) {
			if ((gradA <= 0 && gradB >= 0) || (gradA >= 0 && gradB <= 0)) {
				curveGrad = 0; // point is a local minimum/maximum
			} else {
				if (std::abs(curveGrad) > std::abs(gradA*3)) {
					curveGrad = gradA*3;
				}
				if (std::abs(curveGrad) > std::abs(gradB*3)) {
					curveGrad = gradB*3;
				}
			}
		}
		// When we have duplicate x-values (either side) make up a gradient
		static inline void chooseGradient(Sample &curveGrad, Sample grad1, Sample curveGradOther, Sample y0, Sample y1, bool monotonic) {
			curveGrad = 2*grad1 - curveGradOther;
			if (y0 != y1 && (y1 > y0) != (grad1 >= 0)) { // not duplicate y, but a local min/max
				curveGrad = 0;
			} else if (monotonic) {
				if (grad1 >= 0) {
					curveGrad = std::max<Sample>(0, curveGrad);
				} else {
					curveGrad = std::min<Sample>(0, curveGrad);
				}
			}
		}
	public:
		CubicSegment() : CubicSegment(0, 0, 0, 0, 0) {}
		CubicSegment(Sample xStart, Sample a0, Sample a1, Sample a2, Sample a3) : xStart(xStart), a0(a0), a1(a1), a2(a2), a3(a3) {}
		
		Sample operator ()(Sample x) const {
			x -= xStart;
			return a0 + x*(a1 + x*(a2 + x*a3));
		}
		/// Differentiate
		CubicSegment dx() const {
			return {xStart, a1, 2*a2, 3*a3, 0};
		}
		/// The reference x-value, used as the centre of the cubic expansion
		Sample start() {
			return xStart;
		}
		
		/// Cubic segment based on start/end values and gradients
		static CubicSegment hermite(Sample x0, Sample x1, Sample y0, Sample y1, Sample g0, Sample g1) {
			Sample xScale = 1/(x1 - x0);
			return {
				x0, y0, g0,
				(3*(y1 - y0)*xScale - 2*g0 - g1)*xScale,
				(2*(y0 - y1)*xScale + g0 + g1)*(xScale*xScale)
			};
		}
		
		/** Cubic segment (valid between `x1` and `x2`), which is smooth when applied to an adjacent set of points.
		If `x0 == x1` or `x2 == x3` it will choose a gradient which continues in a quadratic curve, or 0 if the point is a local minimum/maximum.
		*/
		static CubicSegment smooth(Sample x0, Sample x1, Sample x2, Sample x3, Sample y0, Sample y1, Sample y2, Sample y3, bool monotonic=false) {
			if (x1 == x2) return {0, y1, 0, 0, 0}; // zero-width segment, just return constant

			Sample grad1 = gradient(x1, x2, y1, y2);
			Sample curveGrad1 = grad1;
			if (x0 != x1) { // we have a defined x0-x1 gradient
				Sample grad0 = gradient(x0, x1, y0, y1);
				curveGrad1 = (grad0 + grad1)*Sample(0.5);
				if (monotonic) ensureMonotonic(curveGrad1, grad0, grad1);
			} else if (y0 != y1) {
				if ((y1 > y0) != (grad1 >= 0)) curveGrad1 = 0; // set to 0 if it's a min/max
			}
			Sample curveGrad2;
			if (x2 != x3) { // we have a defined x1-x2 gradient
				Sample grad2 = gradient(x2, x3, y2, y3);
				curveGrad2 = (grad1 + grad2)*Sample(0.5);
				if (monotonic) ensureMonotonic(curveGrad2, grad1, grad2);

				if (x0 == x1) { // If the other gradient isn't defined, make one up
					chooseGradient(curveGrad1, grad1, curveGrad2, y0, y1, monotonic);
				}
			} else {
				chooseGradient(curveGrad2, grad1, curveGrad1, y2, y3, monotonic);
			}
			return hermite(x1, x2, y1, y2, curveGrad1, curveGrad2);
		}
	};
	
	/** Smooth interpolation (optionally monotonic) between points, using cubic segments.
	\diagram{cubic-segments-example.svg,Example curve including a repeated point at (2.5\, 3.5) and an instantaneous jump at x=4.  The curve is flat beyond the first/last points.}
	To produce a sharp corner, use a repeated point. The gradient is flat at the edges, unless you use repeated points at the start/end.*/
	template<typename Sample=double>
	class CubicSegmentCurve {
		struct Point {
			Sample x, y;
			bool operator <(const Point &other) const {
				return x < other.x;
			}
		};
		std::vector<Point> points;
		Point first{0, 0}, last{0, 0};
		
		using Segment = CubicSegment<Sample>;
		std::vector<Segment> segments{1};
	public:
		/// Clear existing points and segments
		void clear() {
			points.resize(0);
			segments.resize(0);
			first = last = {0, 0};
		}
	
		/// Add a new point, but does not recalculate the segments.  `corner` just writes the point twice, for convenience.
		CubicSegmentCurve & add(Sample x, Sample y, bool corner=false) {
			points.push_back({x, y});
			if (corner) points.push_back({x, y});
			return *this;
		}
		
		/// Recalculates the segments.
		void update(bool monotonic=false) {
			if (points.empty()) add(0, 0);
			std::stable_sort(points.begin(), points.end()); // Ensure ascending order
			segments.resize(0);
			first = points[0];
			last = points.back();
			for (size_t i = 1; i < points.size(); ++i) {
				Point p1 = points[i - 1];
				Point p2 = points[i];
				if (p1.x != p2.x) {
					Point p0 = (i > 1) ? points[i - 2] : Point{p1.x, p2.y};
					Point p3 = (i + 1 < points.size()) ? points[i + 1] : Point{p2.x, p1.y};
					segments.push_back(Segment::smooth(p0.x, p1.x, p2.x, p3.x, p0.y, p1.y, p2.y, p3.y, monotonic));
				}
			}
		}
		
		/// Reads a value out from the curve.
		Sample operator ()(Sample x) {
			if (x <= first.x) return first.y;
			if (x >= last.x) return last.y;
			size_t index = 1;
			while (index < segments.size() && segments[index].start() <= x) {
				++index;
			}
			return segments[index - 1](x);
		}
	};

/** @} */
}} // signalsmith::envelopes::
#endif // include guard
