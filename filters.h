#ifndef SIGNALSMITH_DSP_FILTERS_H
#define SIGNALSMITH_DSP_FILTERS_H

#include "./common.h"
#include "./perf.h"

#include <cmath>

namespace signalsmith {
namespace filters {
	/**	@defgroup Filters Basic filters
		@brief Classes for some common filter types
		
		@{
		@file
	*/
	
	/** A standard biquad.

		This is not guaranteed to be stable if modulated at audio rate.
		
		The default highpass/lowpass bandwidth (1.9) produces a Butterworth filter when bandwidth-compensation is disabled.
		
		Bandwidth compensation defaults to `true` for all filter types aside from highpass/lowpass.  When in "cookbook mode", it roughly preserves the ratio between upper/lower boundaries (-3dB point or half-gain).  Otherwise, it exactly preserves the lower boundary.*/
	template<typename Sample, bool cookbookBandwidth=false>
	class BiquadStatic {
		Sample a1 = 0, a2 = 0, b0 = 1, b1 = 0, b2 = 0;
		Sample x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		
		// Straight from the cookbook: https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
		enum class Type {highpass, lowpass, highShelf, lowShelf, bandpass, bandStop};
		SIGNALSMITH_INLINE void configure(Type type, double scaledFreq, double octaves, double sqrtGain, bool correctBandwidth) {
			scaledFreq = std::max(0.0001, std::min(0.4999, scaledFreq));
			double w0 = 2*M_PI*scaledFreq;
			double cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
			if (correctBandwidth) {
				if (cookbookBandwidth) {
					// Approximately preserves bandwidth between halfway points
					octaves *= w0/sin_w0;
				} else {
					// Bilinear warps frequencies like tan(pi*x)/pi, so this places the lower boundary in the correct place
					double lowerRatio = std::pow(2, octaves*-0.5);
					octaves = 2*std::log2(std::tan(M_PI*scaledFreq)/std::tan(M_PI*scaledFreq*lowerRatio));
				}
			}
			double alpha = sin_w0*std::sinh(std::log(2)*0.5*octaves);
			double A = sqrtGain, sqrtA2alpha = std::sqrt(A)*alpha;

			double a0;
			if (type == Type::highpass) {
				b1 = -1 - cos_w0;
				b0 = b2 = (1 + cos_w0)*0.5;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::lowpass) {
				b1 = 1 - cos_w0;
				b0 = b2 = b1*0.5;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::highShelf) {
				b0 = A*((A+1)+(A-1)*cos_w0+sqrtA2alpha);
				b2 = A*((A+1)+(A-1)*cos_w0-sqrtA2alpha);
				b1 = -2*A*((A-1)+(A+1)*cos_w0);
				a0 = (A+1)-(A-1)*cos_w0+sqrtA2alpha;
				a2 = (A+1)-(A-1)*cos_w0-sqrtA2alpha;
				a1 = 2*((A-1)-(A+1)*cos_w0);
			} else if (type == Type::lowShelf) {
				b0 = A*((A+1)-(A-1)*cos_w0+sqrtA2alpha);
				b2 = A*((A+1)-(A-1)*cos_w0-sqrtA2alpha);
				b1 = 2*A*((A-1)-(A+1)*cos_w0);
				a0 = (A+1)+(A-1)*cos_w0+sqrtA2alpha;
				a2 = (A+1)+(A-1)*cos_w0-sqrtA2alpha;
				a1 = -2*((A-1)+(A+1)*cos_w0);
			} else if (type == Type::bandpass) {
				b0 = alpha;
				b1 = 0;
				b2 = -alpha;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::bandStop) {
				b0 = 1;
				b1 = -2*cos_w0;
				b2 = 1;
				a0 = 1 + alpha;
				a1 = b1;
				a2 = 1 - alpha;
			} else {
				// reset to neutral
				a1 = a2 = b1 = b2 = 0;
				a0 = b0 = 1;
			}
			double invA0 = 1/a0;
			b0 *= invA0;
			b1 *= invA0;
			b2 *= invA0;
			a1 *= invA0;
			a2 *= invA0;
		}
	public:

		Sample operator ()(Sample x0) {
			Sample y0 = x0*b0 + x1*b1 + x2*b2 - y1*a1 - y2*a2;
			y2 = y1;
			y1 = y0;
			x2 = x1;
			x1 = x0;
			return y0;
		}
		
		void reset() {
			x1 = x2 = y1 = y2 = 0;
		}

		void highpass(double scaledFreq, double octaves=1.9, bool correctBandwidth=false) {
			configure(Type::highpass, scaledFreq, octaves, 0, correctBandwidth);
		}
		void lowpass(double scaledFreq, double octaves=1.9, bool correctBandwidth=false) {
			configure(Type::lowpass, scaledFreq, octaves, 0, correctBandwidth);
		}
		void bandpass(double scaledFreq, double octaves=1, bool correctBandwidth=true) {
			configure(Type::bandpass, scaledFreq, octaves, 0, correctBandwidth);
		}
		void highShelf(double scaledFreq, double gain, double octaves=2, bool correctBandwidth=true) {
			configure(Type::highShelf, scaledFreq, octaves, std::sqrt(gain), correctBandwidth);
		}
		void lowShelf(double scaledFreq, double gain, double octaves=2, bool correctBandwidth=true) {
			configure(Type::lowShelf, scaledFreq, octaves, std::sqrt(gain), correctBandwidth);
		}
		void highShelfDb(double scaledFreq, double db, double octaves=2, bool correctBandwidth=true) {
			double sqrtGain = std::pow(10, db*0.025);
			configure(Type::highShelf, scaledFreq, octaves, sqrtGain, correctBandwidth);
		}
		void lowShelfDb(double scaledFreq, double db, double octaves=2, bool correctBandwidth=true) {
			double sqrtGain = std::pow(10, db*0.025);
			configure(Type::lowShelf, scaledFreq, octaves, sqrtGain, correctBandwidth);
		}
		void bandStop(double scaledFreq, double octaves=1, bool correctBandwidth=true) {
			configure(Type::bandStop, scaledFreq, octaves, 0, correctBandwidth);
		}
	};

/** @} */
}} // signalsmith::filters::
#endif // include guard
