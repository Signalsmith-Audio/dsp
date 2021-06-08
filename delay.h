#ifndef SIGNALSMITH_DSP_DELAY_H
#define SIGNALSMITH_DSP_DELAY_H

#include "./common.h"

#include <vector>
#include <array>
#include <cmath> // for std::ceil()
#include <type_traits>

#include <complex>
#include "./fft.h"
#include "./windows.h"

namespace signalsmith {
namespace delay {
	/**	@defgroup Delay Delay utilities
		@brief Standalone templated classes for delays
		
		You can set up a `Buffer` or `MultiBuffer`, and get interpolated samples using a `Reader` (separately on each channel in the multi-channel case) - or you can use `Delay`/`MultiDelay` which include their own buffers.

		Interpolation quality is chosen using a template class, from @ref Interpolators.

		@{
		@file
	*/

	/** @brief Single-channel delay buffer
 
		Access is used with `buffer[]`, relative to the internal read/write position ("head").  This head is moved using `++buffer` (or `buffer += n`), such that `buffer[1] == (buffer + 1)[0]` in a similar way iterators/pointers.
		
		Operations like `buffer - 10` or `buffer++` return a View, which holds a fixed position in the buffer (based on the read/write position at the time).
		
		The capacity includes both positive and negative indices.  For example, a capacity of 100 would support using any of the ranges:
		
		* `buffer[-99]` to buffer[0]`
		* `buffer[-50]` to buffer[49]`
		* `buffer[0]` to buffer[99]`

		Although buffers are usually used with historical samples accessed using negative indices e.g. `buffer[-10]`, you could equally use it flipped around (moving the head backwards through the buffer using `--buffer`).
	*/
	template<typename Sample>
	class Buffer {
		int bufferIndex;
		int bufferMask;
		std::vector<Sample> buffer;
	public:
		Buffer(int minCapacity=0) {
			resize(minCapacity);
		}
		// We shouldn't accidentally copy a delay buffer
		Buffer(const Buffer &other) = delete;
		Buffer & operator =(const Buffer &other) = delete;

		void resize(int minCapacity, Sample value=Sample()) {
			int bufferLength = 1;
			while (bufferLength < minCapacity) bufferLength *= 2;
			buffer.assign(bufferLength, value);
			bufferMask = bufferLength - 1;
			bufferIndex = 0;
		}
		void reset(Sample value=Sample()) {
			buffer.assign(buffer.size(), value);
		}

		/// Holds a view for a particular position in the buffer
		template<bool isConst>
		class View {
			using CBuffer = typename std::conditional<isConst, const Buffer, Buffer>::type;
			using CSample = typename std::conditional<isConst, const Sample, Sample>::type;
			CBuffer *buffer = nullptr;
			int bufferIndex = 0;
		public:
			View(CBuffer &buffer, int offset=0) : buffer(&buffer), bufferIndex(buffer.bufferIndex + offset) {}
			View(const View &other, int offset=0) : buffer(other.buffer), bufferIndex(other.bufferIndex + offset) {}
			
			CSample & operator[](int offset) {
				return buffer->buffer[(bufferIndex + offset)&buffer->bufferMask];
			}
			const Sample & operator[](int offset) const {
				return buffer->buffer[(bufferIndex + offset)&buffer->bufferMask];
			}

			/// Write data into the buffer
			template<typename Data>
			void write(Data &&data, int length) {
				for (int i = 0; i < length; ++i) {
					(*this)[i] = data[i];
				}
			}
			/// Read data out from the buffer
			template<typename Data>
			void read(int length, Data &&data) const {
				for (int i = 0; i < length; ++i) {
					data[i] = (*this)[i];
				}
			}

			View operator +(int offset) const {
				return View(*this, offset);
			}
			View operator -(int offset) const {
				return View(*this, -offset);
			}
		};
		using MutableView = View<false>;
		using ConstView = View<true>;
		
		MutableView view(int offset=0) {
			return MutableView(*this, offset);
		}
		ConstView view(int offset=0) const {
			return ConstView(*this, offset);
		}
		ConstView constView(int offset=0) const {
			return ConstView(*this, offset);
		}

		Sample & operator[](int offset) {
			return buffer[(bufferIndex + offset)&bufferMask];
		}
		const Sample & operator[](int offset) const {
			return buffer[(bufferIndex + offset)&bufferMask];
		}

		/// Write data into the buffer
		template<typename Data>
		void write(Data &&data, int length) {
			for (int i = 0; i < length; ++i) {
				(*this)[i] = data[i];
			}
		}
		/// Read data out from the buffer
		template<typename Data>
		void read(int length, Data &&data) const {
			for (int i = 0; i < length; ++i) {
				data[i] = (*this)[i];
			}
		}
		
		Buffer & operator ++() {
			++bufferIndex;
			return *this;
		}
		Buffer & operator +=(int i) {
			bufferIndex += i;
			return *this;
		}
		Buffer & operator --() {
			--bufferIndex;
			return *this;
		}
		Buffer & operator -=(int i) {
			bufferIndex -= i;
			return *this;
		}

		MutableView operator ++(int) {
			MutableView view(*this);
			++bufferIndex;
			return view;
		}
		MutableView operator +(int i) {
			return MutableView(*this, i);
		}
		ConstView operator +(int i) const {
			return ConstView(*this, i);
		}
		MutableView operator --(int) {
			MutableView view(*this);
			--bufferIndex;
			return view;
		}
		MutableView operator -(int i) {
			return MutableView(*this, -i);
		}
		ConstView operator -(int i) const {
			return ConstView(*this, -i);
		}
	};

	/** @brief Multi-channel delay buffer

		This behaves similarly to the single-channel `Buffer`, with the following differences:
		
		* `buffer[c]` returns a view for a single channel, which behaves like the single-channel `Buffer::View`.
		* The constructor and `.resize()` take an additional first `channel` argument.
	*/
	template<typename Sample>
	class MultiBuffer {
		int channels, stride;
		Buffer<Sample> buffer;
	public:
		using ConstChannel = typename Buffer<Sample>::ConstView;
		using MutableChannel = typename Buffer<Sample>::MutableView;

		MultiBuffer(int channels=0, int capacity=0) : channels(channels), stride(capacity), buffer(channels*capacity) {}

		void resize(int channels, int capacity, Sample value=Sample()) {
			this->channels = channels;
			stride = capacity;
			buffer.resize(channels*capacity, value);
		}
		void reset(Sample value=Sample()) {
			buffer.reset(value);
		}

		/// A reference-like multi-channel result for a particular sample index
		template<bool isConst>
		class Stride {
			using CChannel = typename std::conditional<isConst, ConstChannel, MutableChannel>::type;
			using CSample = typename std::conditional<isConst, const Sample, Sample>::type;
			CChannel view;
			int channels, stride;
		public:
			Stride(CChannel view, int channels, int stride) : view(view), channels(channels), stride(stride) {}
			
			CSample & operator[](int channel) {
				return view[channel*stride];
			}
			const Sample & operator[](int channel) const {
				return view[channel*stride];
			}

			/// Reads from the buffer into a multi-channel result
			template<class Data>
			void get(Data &&result) const {
				for (int c = 0; c < channels; ++c) {
					result[c] = view[c*stride];
				}
			}
			/// Writes from multi-channel data into the buffer
			template<class Data>
			void set(Data &&data) {
				for (int c = 0; c < channels; ++c) {
					view[c*stride] = data[c];
				}
			}
			template<class Data>
			Stride & operator =(const Data &data) {
				set(data);
				return *this;
			}
			Stride & operator =(const Stride &data) {
				set(data);
				return *this;
			}
		};
		
		Stride<false> at(int offset) {
			return {buffer.view(offset), channels, stride};
		}
		Stride<true> at(int offset) const {
			return {buffer.view(offset), channels, stride};
		}

		/// Holds a particular position in the buffer
		template<bool isConst>
		class View {
			using CChannel = typename std::conditional<isConst, ConstChannel, MutableChannel>::type;
			CChannel view;
			int channels, stride;
		public:
			View(CChannel view, int channels, int stride) : view(view), channels(channels), stride(stride) {}
			
			CChannel operator[](int channel) {
				return view + channel*stride;
			}
			ConstChannel operator[](int channel) const {
				return view + channel*stride;
			}

			Stride<isConst> at(int offset) {
				return {view + offset, channels, stride};
			}
			Stride<true> at(int offset) const {
				return {view + offset, channels, stride};
			}
		};
		using MutableView = View<false>;
		using ConstView = View<true>;

		MutableView view(int offset=0) {
			return MutableView(buffer.view(offset), channels, stride);
		}
		ConstView view(int offset=0) const {
			return ConstView(buffer.view(offset), channels, stride);
		}
		ConstView constView(int offset=0) const {
			return ConstView(buffer.view(offset), channels, stride);
		}

		MutableChannel operator[](int channel) {
			return buffer + channel*stride;
		}
		ConstChannel operator[](int channel) const {
			return buffer + channel*stride;
		}
		
		MultiBuffer & operator ++() {
			++buffer;
			return *this;
		}
		MultiBuffer & operator +=(int i) {
			buffer += i;
			return *this;
		}
		MutableView operator ++(int) {
			return MutableView(buffer++, channels, stride);
		}
		MutableView operator +(int i) {
			return MutableView(buffer + i, channels, stride);
		}
		ConstView operator +(int i) const {
			return ConstView(buffer + i, channels, stride);
		}
		MultiBuffer & operator --() {
			--buffer;
			return *this;
		}
		MultiBuffer & operator -=(int i) {
			buffer -= i;
			return *this;
		}
		MutableView operator --(int) {
			return MutableView(buffer--, channels, stride);
		}
		MutableView operator -(int i) {
			return MutableView(buffer - i, channels, stride);
		}
		ConstView operator -(int i) const {
			return ConstView(buffer - i, channels, stride);
		}
	};
	
	/** \defgroup Interpolators Interpolators
		\ingroup Delay
		@{ */
	template<typename Sample>
	struct InterpolatorNearest {
		static constexpr int inputLength = 1;
		static constexpr Sample latency = -0.5; // Because we're truncating, which rounds down too often
	
	protected:
		template<class Data>
		static Sample fractional(const Data &data, Sample) {
			return data[0];
		}
	};
	template<typename Sample>
	struct InterpolatorLinear {
		static constexpr int inputLength = 2;
		static constexpr int latency = 0;
	
	protected:
		template<class Data>
		static Sample fractional(const Data &data, Sample fractional) {
			Sample a = data[0], b = data[1];
			return a + fractional*(b - a);
		}
	};
	// Spline cubic
	template<typename Sample>
	struct InterpolatorCubic {
		static constexpr int inputLength = 4;
		static constexpr int latency = 1;
	
	protected:
		template<class Data>
		static Sample fractional(const Data &data, Sample fractional) {
			// Cubic interpolation
			Sample a = data[0], b = data[1], c = data[2], d = data[3];
			Sample cbDiff = c - b;
			Sample k1 = (c - a)*0.5;
			Sample k3 = k1 + (d - b)*0.5 - cbDiff*2;
			Sample k2 = cbDiff - k3 - k1;
			return b + fractional*(k1 + fractional*(k2 + fractional*k3)); // 16 ops total, not including the indexing
		}
	};
	template<typename Sample, int n>
	struct InterpolatorLagrangeN {
		static constexpr int inputLength = n + 1;
		static constexpr int latency = (n - 1)/2;

		std::array<Sample, (n + 1)> invDivisors;
		
		InterpolatorLagrangeN() {
			for (int j = 0; j <= n; ++j) {
				double divisor = 1;
				for (int k = 0; k < j; ++k) divisor *= (j - k);
				for (int k = j + 1; k <= n; ++k) divisor *= (j - k);
				invDivisors[j] = 1/divisor;
			}
		}
		
		template<class Data>
		Sample fractional(const Data &data, Sample fractional) const {
			std::array<Sample, (n + 1)> sums;
			
			Sample x = fractional + latency;

			Sample forwardFactor = 1;
			sums[0] = data[0];
			for (int i = 1; i <= n; ++i) {
				forwardFactor *= x - (i - 1);
				sums[i] = forwardFactor*data[i];
			}
			
			Sample backwardsFactor = 1;
			Sample result = sums[n]*invDivisors[n];
			for (int i = n - 1; i >= 0; --i) {
				backwardsFactor *= x - (i + 1);
				result += sums[i]*invDivisors[i]*backwardsFactor;
			}
			return result;
		}
	};
	
	template<typename Sample>
	using InterpolatorLagrange3 = InterpolatorLagrangeN<Sample, 3>;
	template<typename Sample>
	using InterpolatorLagrange7 = InterpolatorLagrangeN<Sample, 7>;
	template<typename Sample>
	using InterpolatorLagrange19 = InterpolatorLagrangeN<Sample, 19>;

	template<typename Sample, int n, bool minimumPhase=false>
	struct InterpolatorKaiserSincN {
		static constexpr int inputLength = n;
		static constexpr Sample latency = minimumPhase ? 0 : (n*Sample(0.5) - 1);

		int subSampleSteps;
		std::vector<Sample> coefficients;

		InterpolatorKaiserSincN(double bandwidthTarget=0) {
			subSampleSteps = 2*n; // Heuristic again.  Really it depends on the bandwidth as well.
			if (bandwidthTarget == 0) {
				bandwidthTarget = 1 - 0.9/std::sqrt(n);
			}
			double kaiserBandwidth = (1 - bandwidthTarget)*(n + 1.0/subSampleSteps);
			kaiserBandwidth += 1.25/kaiserBandwidth; // We want to place the first zero, but (because using this to window a sinc essentially integrates it in the freq-domain), our ripples (and therefore zeroes) are out of phase.  This is a heuristic fix.
			
			double centreIndex = n*subSampleSteps*0.5, scaleFactor = 1.0/subSampleSteps;
			std::vector<Sample> windowedSinc(subSampleSteps*n + 1);
			
			::signalsmith::windows::Kaiser::withBandwidth(kaiserBandwidth, false).fill(windowedSinc, windowedSinc.size());

			for (size_t i = 0; i < windowedSinc.size(); ++i) {
				double x = (i - centreIndex)*scaleFactor;
				int intX = std::round(x);
				if (intX != 0 && std::abs(x - intX) < 1e-6) {
					// Exact 0s
					windowedSinc[i] = 0;
				} else if (std::abs(x) > 1e-6) {
					double p = x*M_PI;
					windowedSinc[i] *= std::sin(p)/p;
				}
			}
			
			if (minimumPhase) {
				signalsmith::FFT<Sample> fft(windowedSinc.size()*2, 1);
				windowedSinc.resize(fft.size(), 0);
				std::vector<std::complex<Sample>> spectrum(fft.size());
				std::vector<std::complex<Sample>> cepstrum(fft.size());
				fft.fft(windowedSinc, spectrum);
				for (size_t i = 0; i < fft.size(); ++i) {
					spectrum[i] = std::log(std::abs(spectrum[i]) + 1e-30);
				}
				fft.fft(spectrum, cepstrum);
				for (size_t i = 1; i < fft.size()/2; ++i) {
					cepstrum[i] *= 0;
				}
				for (size_t i = fft.size()/2 + 1; i < fft.size(); ++i) {
					cepstrum[i] *= 2;
				}
				Sample scaling = Sample(1)/fft.size();
				fft.ifft(cepstrum, spectrum);

				for (size_t i = 0; i < fft.size(); ++i) {
					Sample phase = spectrum[i].imag()*scaling;
					Sample mag = std::exp(spectrum[i].real()*scaling);
					spectrum[i] = {mag*std::cos(phase), mag*std::sin(phase)};
				}
				fft.ifft(spectrum, cepstrum);
				windowedSinc.resize(subSampleSteps*n + 1);
				for (size_t i = 0; i < windowedSinc.size(); ++i) {
					windowedSinc[i] = cepstrum[i].real()*scaling;
				}
			}
			
			// Re-order into FIR fractional-delay blocks
			coefficients.resize(n*(subSampleSteps + 1));
			for (int k = 0; k <= subSampleSteps; ++k) {
				for (int i = 0; i < n; ++i) {
					coefficients[k*n + i] = windowedSinc[(subSampleSteps - k) + i*subSampleSteps];
				}
			}
		}
		
		template<class Data>
		Sample fractional(const Data &data, Sample fractional) const {
			Sample subSampleDelay = fractional*subSampleSteps;
			int lowIndex = subSampleDelay;
			if (lowIndex >= subSampleSteps) lowIndex = subSampleSteps - 1;
			Sample subSampleFractional = subSampleDelay - lowIndex;
			int highIndex = lowIndex + 1;
			
			Sample sumLow = 0, sumHigh = 0;
			const Sample *coeffLow = coefficients.data() + lowIndex*n;
			const Sample *coeffHigh = coefficients.data() + highIndex*n;
			for (int i = 0; i < n; ++i) {
				sumLow += data[i]*coeffLow[i];
				sumHigh += data[i]*coeffHigh[i];
			}
			return sumLow + (sumHigh - sumLow)*subSampleFractional;
		}
	};

	template<typename Sample>
	using InterpolatorKaiserSinc20 = InterpolatorKaiserSincN<Sample, 20>;
	template<typename Sample>
	using InterpolatorKaiserSinc8 = InterpolatorKaiserSincN<Sample, 8>;
	template<typename Sample>
	using InterpolatorKaiserSinc4 = InterpolatorKaiserSincN<Sample, 4>;

	template<typename Sample>
	using InterpolatorKaiserSinc20Min = InterpolatorKaiserSincN<Sample, 20, true>;
	template<typename Sample>
	using InterpolatorKaiserSinc8Min = InterpolatorKaiserSincN<Sample, 8, true>;
	template<typename Sample>
	using InterpolatorKaiserSinc4Min = InterpolatorKaiserSincN<Sample, 4, true>;
	///  @}
	
	/** @brief A delay-line reader which uses an external buffer
 
		This is useful if you have multiple delay-lines reading from the same buffer.
	*/
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class Reader : public Interpolator<Sample> /* so we can get the empty-base-class optimisation */ {
		using Super = Interpolator<Sample>;
	public:
		Reader () {}
		/// Pass in a configured interpolator
		Reader (const Interpolator<Sample> &interpolator) : Super(interpolator) {}
	
		template<typename Buffer>
		Sample read(const Buffer &buffer, Sample delaySamples) const {
			int startIndex = delaySamples;
			Sample remainder = delaySamples - startIndex;
			
			// Delay buffers use negative indices, but interpolators use positive ones
			using View = decltype(buffer - startIndex);
			struct Flipped {
				 View view;
				 Sample operator [](int i) const {
					return view[-i];
				 }
			};
			return Super::fractional(Flipped{buffer - startIndex}, remainder);
		}
	};

	/**	@brief A single-channel delay-line containing its own buffer. */
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class Delay : private Reader<Sample, Interpolator> {
		using Super = Reader<Sample, Interpolator>;
		Buffer<Sample> buffer;
	public:
		static constexpr Sample latency = Super::latency;

		Delay(int capacity=0) : buffer(capacity + Super::inputLength) {}
		/// Pass in a configured interpolator
		Delay(const Interpolator<Sample> &interp, int capacity=0) : Super(interp), buffer(capacity + Super::inputLength) {}
		
		void reset(Sample value=Sample()) {
			buffer.reset(value);
		}
		void resize(int minCapacity, Sample value=Sample()) {
			buffer.resize(minCapacity + Super::inputLength, value);
		}
		
		Sample read(Sample delaySamples) const {
			return Super::read(buffer, delaySamples);
		}
		void write(Sample value) {
			buffer[0] = value;
			++buffer;
		}
		Sample readWrite(Sample value, Sample delaySamples) {
			buffer[0] = value;
			Sample result = Super::read(buffer, delaySamples);
			++buffer;
			return result;
		}
	};

	/**	@brief A multi-channel delay-line with its own buffer. */
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class MultiDelay : private Reader<Sample, Interpolator> {
		using Super = Reader<Sample, Interpolator>;
		int channels;
		MultiBuffer<Sample> multiBuffer;
	public:
		static constexpr Sample latency = Super::latency;

		MultiDelay(int channels=0, int capacity=0) : channels(channels), multiBuffer(channels, capacity + Super::inputLength) {}

		void reset(Sample value=Sample()) {
			multiBuffer.reset(value);
		}
		void resize(int channels, int capacity, Sample value=Sample()) {
			this->channels = channels;
			multiBuffer.resize(channels, capacity + Super::inputLength, value);
		}
		
		/// A single-channel delay-line view, similar to a `const Delay`
		struct ChannelView {
			static constexpr Sample latency = Super::latency;

			Super &reader;
			typename MultiBuffer<Sample>::ConstChannel &channel;
			
			Sample read(Sample delaySamples) const {
				return reader.read(channel, delaySamples);
			}
		};
		ChannelView operator [](int channel) {
			return ChannelView{*this, multiBuffer[channel]};
		}

		/// A multi-channel result, lazily calculating samples
		struct DelayView {
			Super &reader;
			typename MultiBuffer<Sample>::ConstView view;
			Sample delaySamples;
			
			// Calculate samples on-the-fly
			Sample operator [](int c) const {
				return reader.read(view[c], delaySamples);
			}
		};
		DelayView read(Sample delaySamples) {
			return DelayView{*this, multiBuffer.constView(), delaySamples};
		}
		template<class Data>
		void write(const Data &data) {
			for (int c = 0; c < channels; ++c) {
				multiBuffer[c][0] = data[c];
			}
			++multiBuffer;
		}
		template<class Data>
		DelayView readWrite(const Data &data, Sample delaySamples) {
			for (int c = 0; c < channels; ++c) {
				multiBuffer[c][0] = data[c];
			}
			DelayView result = read(delaySamples);
			write(data);
			return result;
		}
	};

/** @} */
}} // signalsmith::delay::
#endif // SIGNALSMITH_DSP_DELAY_H
